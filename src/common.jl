using Glob
using Dates
using CategoricalArrays
using NetCDF
using StatsBase
using FLoops

using
    CSV,
    Dates,
    DataFrames,
    YAXArrays,
    DimensionalData,
    DataStructures

using
    ColorSchemes,
    GeoMakie,
    # GLMakie,
    GraphMakie

using
    Statistics,
    Bootstrap,
    LinearAlgebra

using Graphs, SimpleWeightedGraphs
import Graphs.Parallel
import GeoDataFrames as GDF
import GeoFormatTypes as GFT
import ArchGDAL as AG
import GeoInterface as GI
using ADRIA

include("plotting_functions.jl")
include("ADRIA_DataCube_functions.jl")

gbr_domain_path = "c:/Users/bgrier/Documents/Projects/ADRIA_Domains/rme_ml_2024_01_08/"
conn_path = joinpath(gbr_domain_path, "data_files/con_bin")

"""
    find_intersections(
        x::DataFrame,
        y::DataFrame,
        x_id::Symbol,
        y_id::Symbol,
        y_geom_col::Symbol=:geometry;
        proportion::Bool=false,
    )::DataFrame

Find the areas of `y` that intersect with each polygon in `x`.
`rel_areas` contains corresponding `y_id` for each intersecting polygon in x (can then be
joined to `x`).

If `proportion = true`: polygons of `y` are only chosen if the intersection with `x` is >
50% the area of `x`.

# Arguments
- `x` : The target GeoDataFrame to compare with
- `y` : GeoDataFrame containing polygons to match against `x`
- `xid` : Column name holding unique IDs for x geometries (referred to as GBRMPA_ID in rel_areas)
- `yid` : Column name holding variable of interest for y geometries
- `y_geom_col` : Column name holding geometries in y
- `proportion` : Only select y polygons if the intersection with x polygon is > 50% of x polygon area
                 (default: `false`).
"""
function find_intersections(
    x::DataFrame,
    y::DataFrame,
    x_id::Symbol,
    y_id::Symbol,
    y_geom_col::Symbol=:geometry;
    proportion::Bool=false
)::DataFrame
    rel_areas = DataFrame(
        [Vector{Any}(missing, size(x, 1)) for _ in 1:2],
        [:GBRMPA_ID, :area_ID]
    )

    for (x_i, reef_poly) in enumerate(eachrow(x))
        intersecting = DataFrame(
            [Vector{Any}(missing, size(y, 1)) for _ in 1:3],
            [:GBRMPA_ID, :area_ID, :inter_area]
        )

        for (y_i, interest_area) in enumerate(eachrow(y))
            if AG.intersects(reef_poly.geometry, interest_area[y_geom_col])
                inter_area = AG.intersection(
                    reef_poly.geometry, interest_area[y_geom_col]
                )

                inter_area = AG.geomarea(inter_area)
                if proportion
                    prop_area = inter_area / AG.geomarea(reef_poly.geometry)

                    if prop_area >= 0.5
                        data = [reef_poly[x_id], interest_area[y_id], inter_area]

                    else
                        data = [missing, missing, missing]
                    end
                else
                    data = [reef_poly[x_id], interest_area[y_id], inter_area]
                end
            else
                data = [reef_poly[x_id], missing, missing]
            end

            intersecting[y_i, :] = data
        end

        if all(ismissing, intersecting.area_ID)
            x_data = [intersecting[1, x_id], intersecting[1, :area_ID]]
        else
            dropmissing!(intersecting)
            max_inter_area = argmax(intersecting.inter_area)
            x_data = [intersecting[max_inter_area, x_id], intersecting[max_inter_area, :area_ID]]
        end

        rel_areas[x_i, :] = x_data
    end

    return rel_areas
end

DATE_FORMAT = "YYYY-mm-dd-THH-MM-SS"

"""
    _get_file_timestamp(file_path, dt_length)::DateTime

Extract the timestamp from a given file name.
"""
function _get_file_timestamp(file_path, dt_length)::DateTime
    # Get name of file without extension
    filename = splitext(basename(file_path))[1]

    local fn_timestamp
    try
        fn_timestamp = Dates.DateTime(filename[end-(dt_length-1):end], DATE_FORMAT)
    catch err
        if !(err isa ArgumentError)
            rethrow(err)
        end

        # Otherwise, some unexpected date format was encountered so we assign an
        # very early date/time.
        fn_timestamp = Dates.DateTime("1900-01-01-T00-00-00", DATE_FORMAT)
    end

    # Return datetime stamp
    return fn_timestamp
end

"""
    find_latest_file(
        target_dir::String;
        prefix::String="rrap_canonical",
        ext::String="gpkg"
    )::String

Identify the latest output file in a directory based on the timestamp included in the
file name (default: `YYYY-mm-dd-THH-MM-SS`). Intended to find the latest output file for
input into the next script.

# Arguments
- `target_dir` : Target directory
- `prefix` : prefix of target file
- `ext` : the file extension

# Returns
Path to latest output file.
"""
function find_latest_file(
    target_dir::String;
    prefix::String="rrap_canonical",
    ext::String="gpkg",
    DATE_FORMAT=DATE_FORMAT
)::String
    # Get list of files matching the given pattern
    candidate_files = glob("$(prefix)*.$(ext)", target_dir)

    timestamps = map(f -> _get_file_timestamp(f, length(DATE_FORMAT)), candidate_files)
    latest = candidate_files[argmax(timestamps)]

    return latest
end

"""
    relative_site_cover(x)

Calculate the cover of a site at each timestep of x relative to the site's initial cover (x[1]) (standardise all series to start at 1).
"""
function relative_site_cover(x)
    init = x[1]
    for (index, step) in enumerate(x)
        x[index] = x[index]/init
    end

    return x
end

canonical_reefs = find_latest_file("../../canonical-reefs/output/")
canonical_reefs = GDF.read(canonical_reefs)

# Get location indices
REGION_REEFS = DataStructures.OrderedDict(
    "Far Northern Management Area"=>Int64[],
    "Cairns/Cooktown Management Area"=>Int64[],
    "Townsville/Whitsunday Management Area"=>Int64[],
    "Mackay/Capricorn Management Area"=>Int64[],
)

for mgmt_area in collect(keys(REGION_REEFS))
    if mgmt_area == "NA"
        continue
    end

    REGION_REEFS[mgmt_area] = findall(canonical_reefs.management_area .== mgmt_area)
end

"""
    connectivity_scoring(
        conn_matrix::YAXArray;
        gdf::DataFrame=nothing,
        context_layer::Symbol=nothing,
        conn_col_name::Symbol=nothing,
        by_layer::Bool=false
    )::DataFrame

Calculate eigenvector_centrality connectivity scores for GBR reefs at different spatial scales.
When by_layer=true and a DataFrame is given for gdf and a layer symbol is given for context_layer,
then reefs are subset into categories via context_layer groups and then eigenvectors are calculated.

# Arguments
- `conn_matrix` : YAXArray containing the connectivity values. Diagonal should be set to 0 prior to use.
- `gdf` : DataFrame containing context_layer.
- `context_layer` : Categorical column with levels for grouping.
- `by_layer` : Calculate by grouping reefs by context_layer (vs whole GBR connectivity).
"""
function connectivity_scoring(
    conn_matrix::YAXArray;
    gdf=nothing,
    context_layer=nothing,
    conn_col_name=nothing,
    by_layer::Bool=false
)::DataFrame
    RME_UNIQUE_ID = collect(getAxis("Source", conn_matrix).val)
    connectivity_scores = DataFrame(
        [RME_UNIQUE_ID, Vector{Union{Missing, Float64}}(missing, 3806)],
        [:RME_UNIQUE_ID, :conn_score]
    )

    if by_layer
        for level in unique(gdf[:, context_layer])
            level_reefs = gdf[gdf[:, context_layer] .== level, :RME_UNIQUE_ID]
            level_matrix = conn_matrix[conn_matrix.Source .∈ [level_reefs], conn_matrix.Sink .∈ [level_reefs]]

            g = SimpleWeightedDiGraph(level_matrix)
            conn_scores = Dict(zip(collect(level_matrix.Source), eigenvector_centrality(g)))

            for (k, v) in conn_scores
                connectivity_scores[connectivity_scores.RME_UNIQUE_ID .== k, :conn_score] .= v
            end
        end

        rename!(connectivity_scores, :conn_score => conn_col_name)
    else
        g = SimpleWeightedDiGraph(conn_matrix)
        conn_score = eigenvector_centrality(g)
        connectivity_scores.conn_score = conn_score
    end

    return connectivity_scores
end

"""
    weight_by_context(
        gdf::DataFrame,
        target_col::Symbol,
        context_layer::Symbol,
        new_col_name::Symbol
    )::DataFrame

Rank the values in the target_col by their numerical order within each level of context_layer category.

# Arguments
- `gdf` : DataFrame containing all columns.
- `target_col` : Column containing values for ranking.
- `context_layer` : Categorical column with levels for ranking.
- `new_col_name` : Column name for new ranked values.
"""
function weight_by_context(
    gdf::DataFrame,
    target_col::Symbol,
    context_layer::Symbol,
    new_col_name::Symbol
    )::DataFrame

    gdf[:, new_col_name] .= 0.0
    for level in unique(gdf[:, context_layer])
        gdf_level = gdf[gdf[:, context_layer] .== level, :]

        max_target = maximum(gdf_level[:, target_col])
        gdf[gdf[:, context_layer] .== level, new_col_name] = gdf_level[:, target_col] ./ max_target
    end

    return gdf
end

"""
    normalise(x, (a, b))

Normalise the vector `x` so that it's minimum value is `a` and its maximum value is `b`.

# Examples
- `normalise([1, 5, 10, 78] (0,1))` to return a vector with min=0 and max=1.
- `normalise([1, 5, 10, 78] (-1,1))` to return a vector with min=-1 and max=1.
"""
function normalise(x, (a, b))
    x_norm = (b - a) .* ((x .- minimum(filter(!isnan,x))) ./ (maximum(filter(!isnan,x)) .- minimum(filter(!isnan,x)))) .+ a
    return x_norm
end

function _coral_evenness(r_taxa_cover::AbstractArray{T}; method="scaled_evenness_additive", evenness_weight=1, cover_weight=1)::AbstractArray{T} where {T<:Real}
    # Evenness as a functional diversity metric
    n_steps, n_locs, n_grps = size(r_taxa_cover)

    # Sum across groups represents functional diversity
    # Group evenness (Hill 1973, Ecology 54:427-432)
    loc_cover = dropdims(sum(r_taxa_cover, dims=3), dims=3)
    simpsons_diversity::YAXArray = ZeroDataCube((:timesteps, :locations), (n_steps, n_locs))

    if method == "scaled_evenness_multiplicative"
        for loc in axes(loc_cover, 2)
            norm_evenness = normalise((1.0 ./ sum((r_taxa_cover[:, loc, :] ./ loc_cover[:, loc]) .^ 2, dims=2)), (0,1))
            norm_loc_cover = normalise(loc_cover[:, loc], (0,1))
            simpsons_diversity[:, loc] = norm_evenness .* norm_loc_cover
        end
    elseif method == "scaled_evenness_additive"
        for loc in axes(loc_cover, 2)
            norm_evenness = normalise((1.0 ./ sum((r_taxa_cover[:, loc, :] ./ loc_cover[:, loc]) .^ 2, dims=2)), (0,1))
            norm_loc_cover = normalise(loc_cover[:, loc], (0,1))
            simpsons_diversity[:, loc] = (evenness_weight .* norm_evenness) .+ (cover_weight .* norm_loc_cover)
        end
    elseif method == "normalised_evenness"
        for loc in axes(loc_cover, 2)
            simpsons_diversity[:, loc] = normalise((1.0 ./ sum((r_taxa_cover[:, loc, :] ./ loc_cover[:, loc]) .^ 2, dims=2)), (0,1))
        end
    elseif method == "raw_evenness"
        for loc in axes(loc_cover, 2)
            simpsons_diversity[:, loc] = (1.0 ./ sum((r_taxa_cover[:, loc, :] ./ loc_cover[:, loc]) .^ 2, dims=2))
        end
    elseif method == "shannon_index"
        for loc in axes(loc_cover, 2)
            simpsons_diversity[:, loc] = 1.0 ./ sum((r_taxa_cover[:, loc, :] ./ loc_cover[:, loc]) .* log.(r_taxa_cover[:, loc, :] ./ loc_cover[:, loc]), dims=2)
        end
    end

    return replace!(simpsons_diversity, NaN => 0.0, Inf => 0.0) ./ n_grps
end

"""
    extract_timeseries(rs_YAXArray, reefs, context_cols)

Extract the timeseries data for each reef in `reefs` dataframe and attach `context_cols` from
`reefs` to the output dataframe.

# Arguments
- `rs_YAXArray` : YAXArray containing sites that are RME_UNIQUE_IDs and timeseries data for 78 year timeseries.
- `reefs` : Context dataframe containing the target reefs to keep and their context columns.
- `context_cols` : Names of desired context columns for attaching to output timeseries dataframe
"""
function extract_timeseries(rs_YAXArray, reefs, context_cols)
    local reef_ids = []
    try
        reef_ids = collect(getAxis("locations", rs_YAXArray).val)
    catch
        @warn "locations dimension not found in YAXArray, attaching from `reefs`"
        reef_ids = reefs.UNIQUE_ID
    end

    df = DataFrame(rs_YAXArray.data, reef_ids)
    df.year = [string(i) for i in 1:size(df,1)]
    select!(df, :year, Not(:year))

    data = permutedims(df, 1, "UNIQUE_ID")
    data = data[data.UNIQUE_ID .∈ [reefs.UNIQUE_ID],:]
    data = leftjoin(data, reefs[:, vcat(["UNIQUE_ID", context_cols]...)], on=:UNIQUE_ID)
    data = dropmissing(data)

    return data
end

"""
    stat_range(x::AbstractArray)

Return the range of values in `x`.
"""
function stat_range(x::AbstractArray)
    isempty(x) && error("Does not support empty vectors")
    min = typemax(eltype(x))
    max = typemin(eltype(x))
    for xi in x
        isnan(xi) && error("Does not support NaNs in vectors")
        min = ifelse(min > xi, xi, min)
        max = ifelse(max < xi, xi, max)
    end
    return max - min
end

"""
    grouping_counts(
        grouping_col::Union{Symbol, String}, 
        dataset::DataFrame, 
        clustering_col::Union{Symbol, String}, 
        n_clusters::Int64,
        n_reefs::Int64
    )::Vector{String}

Calculate the number of clusters within each grouping. If the number of clusters is less
than expected (`n_clusters`) then a group is returned in `removed_groups`. Calculate the number of reefs
within each cluster of each group. If the number of reefs is less than desired (`n_reefs`) then
a group is returned in `removed_groups`.

# Arguments
- `grouping_col` : Name of column containing group assignments.
- `dataset` : Dataframe containing `grouping_col` and `clustering_col` information.
- `clustering_col` : Name of column containing clustered assignments.
- `n_clusters` : Desired minimum number of unique clusters within each group.
- `n_reefs` : Desired minimum number of reefs within each cluster of each group.

# Examples
```
groups_too_few_clusters = grouping_counts(
    :bioregion,
    context_layers,
    :bioregion_timeseries_clusters,
    3,
    5
)
```
"""
function grouping_counts(
    grouping_col::Union{Symbol, String}, 
    dataset::DataFrame, 
    clustering_col::Union{Symbol, String}, 
    n_clusters::Int64,
    n_reefs::Int64
)::Vector{String}
    reef_counts = combine(groupby(dataset, [String(grouping_col), clustering_col]), nrow => :nrow)
    groups_few_reefs = unique(reef_counts[reef_counts.nrow .< n_reefs, grouping_col])

    cluster_counts = combine(groupby(dataset, grouping_col)) do sdf
        return (clusters = length(unique(sdf[:, clustering_col])), total_reefs=nrow(sdf))
    end
    groups_few_clusters = unique(cluster_counts[cluster_counts.clusters .< n_clusters, grouping_col])
    
    removed_groups = unique(vcat(groups_few_reefs, groups_few_clusters))

    return removed_groups
end

"""
    GCM_analysis_results(results_set::ADRIA.ADRIAResultSet)::Dataset

Use `results_set` to calculate metrics for further analyses. These metrics include
total coral cover for each reef and scenario (`absolute_scenario_cover`), median total
coral cover timeseries (`absolute_median_cover`) and coral cover relative to reef initial
cover (`relative_cover`).
"""
function GCM_analysis_results(results_set::ADRIA.ADRIAResultSet)::Dataset
    
    # Extract metric from scenarios
    tac = ADRIA.metrics.total_absolute_cover(results_set)

    # Get a timeseries summarizing the scenarios for each site
    absolute_cover = ADRIA.metrics.loc_trajectory(median, tac)

    # # Calculate scenario taxa evenness and then summarise to median location level
    # taxa_evenness = ADRIA.metrics.coral_evenness(rs)
    # taxa_evenness = ADRIA.metrics.loc_trajectory(median, taxa_evenness)

    # Calculate cover relative to the reef's initial coral cover
    rel_cover = Float64.(mapslices(relative_site_cover, absolute_cover[:, :], dims=[:timesteps]))

    return Dataset(;
        Dict(
            :absolute_scenario_cover => tac,
            :absolute_median_cover => absolute_cover,
            :relative_cover => rel_cover
        )...
    )
end

"""
    grouped_timeseries_clustering(
        timeseries::AbstractMatrix, 
        grouping::Vector; 
        n_clusters::Int64=3, 
        length_t=1:size(timeseries, 1)
    )::Tuple{Vector{Int64}, Vector{String}}

Perform ADRIA timeseries clustering within each group of `grouping`. Return two vectors
matching the location-order of `timeseries` and `grouping` containing the integer cluster
assignments and the matching coral cover levels.
Function is designed to assign 3 clusters in each group and label these clusters `low`, 
`medium` or `high` depending on their relative coral cover levels between timesteps 6:26 of
the timeseries (2030:2050).
Cover-labelling is performed based on these years of timeseries because there are commonly
large declines in cover in years 1-5 and low levels of cover in years post-2050.
"""
function grouped_timeseries_clustering(
    timeseries::AbstractMatrix, 
    grouping::Vector;
    n_clusters::Int64=3, 
    length_t=1:size(timeseries, 1)
)::Tuple{Vector{Int64}, Vector{String}}
    # clustering = Dict{String, Any}(group => missing for group in grouping)
    cluster_assignemnts = zeros(Int64, length(grouping))
    cluster_cats = Vector{Union{String, Missing}}(missing, length(grouping))

    for group in unique(grouping)
        indices = findall(grouping .== group)
        group_ts = timeseries[:, indices]

        clusters = ADRIA.analysis.cluster_series(group_ts[length_t, :], n_clusters)
        unique_clusters = unique(clusters)

        # Find which numerical cluster belongs to each 'high, medium or low' group. These can differ between bioregions
        # so we need to find which cluster number belongs to the high, med or low median group. This is just based on 
        # years 2030:2050 (time indices 6:26)
        cluster_medians = [median(group_ts[6:26, findall(clusters .== cluster)]) for cluster in unique_clusters]
        cluster_categories = Dict(
            "high" => unique_clusters[argmax(cluster_medians)],
            "low" => unique_clusters[argmin(cluster_medians)],
        )
        merge!(
            cluster_categories, 
            Dict("medium" => 
            unique_clusters[
                findfirst(unique_clusters .∉ [[cluster_categories["high"], cluster_categories["low"]]])
                ]
            )
        )
        cluster_categories = [find_dict_val(value, cluster_categories) for value in clusters]
        cluster_cats[indices] = cluster_categories

        # Once we know which reefs are 'high, medium or low' we can then number them accordingly so the labels
        # are consistent across regions. Low = 1, Medium = 2, High = 3.
        cluster_numbers = Dict("high" => 3, "medium" => 2, "low" => 1)
        cluster_assignemnts[indices] = [cluster_numbers[key] for key in cluster_categories]
        
    end

    return (cluster_assignemnts, cluster_cats)
end

"""
    find_dict_val(val, dict::Dict)

Identify the key in `dict` that holds the target `val`. If there are multiple matches for
`val` the first matching key is returned. If no matching values are contained in `dict` then
nothing is returned.
"""
function find_dict_val(val, dict::Dict)
    for (k,v) in dict
        if v == val
            return k
        end
    end

    return nothing
end

"""
    change_ADRIA_debug(to::Bool; config_fn::String="config.toml")::Nothing

Edit ADRIA debug mode specified in the `config.toml` to become `to`.
ADRIA must be reloaded after after this change is made to take effect.

# Arguments
- `to` : Desired debug option to set in file. (possible options = true, false)
- `config_fn` : path to config.toml file. Default path is "./config.toml"
"""
function change_ADRIA_debug(to::Bool; config_fn::String="config.toml")::Nothing
    config_contents = TOML.parsefile(config_fn)

    if (config_contents["operation"]["debug"] != to)
        config_contents["operation"]["debug"] = to
        
        open(config_fn, "w") do io
            TOML.print(io, config_contents)
        end

        @info "config debug set to $(to)"
    else
        @info "config debug already set to $(to)"
    end
    return nothing
end

"""
    cross_correlation(
        x::AbstractVector{<:Real},
        y::AbstractVector{<:Real},
        lags::AbstractVector{<:Integer},
        demean::Bool
    )

Calculate the normalised cross correlation of two vectors x and y with time series
lags. If `x` is ahead of `y` then a positive lag will result in positive correlation. If `y`
is ahead of `x`, then a negative lag will result in positive correlation.
E.g. If testing for x reef to be ahead of y reef, test for correlation at positive lag.

Based on StatsBase https://github.com/JuliaStats/StatsBase.jl/blob/60fb5cd400c31d75efd5cdb7e4edd5088d4b1229/src/signalcorr.jl#L400
and https://paulbourke.net/miscellaneous/correlate/

# Arguments
- `x` : Vector of interest to test for being ahead or behind `y`
- `y` : Vector to test lags of `x` against
- `lags` : Vector of lags to apply to vector. Positive lags test for `x` leading `y`, negative lags test for `y` leading `x`.
- `demean` : Subtract the mean of each vector from each element of `x` and `y`. If demean is intended include it as true, otherwise do not include `demean` argument.

# Returns
Vector of correlation values for each lag in `lags`.
"""
function cross_correlation(
    x::AbstractVector{<:Real},
    y::AbstractVector{<:Real},
    lag::Int64,
    correlation_function::Function,
    demean::Bool
    )

    #r = Vector{Float64}()
    lx = length(x)
    #m = length(lags)

    if demean
        zx::Vector{Float64} = x .- mean(x)
    else
        throw("`demean` must be true if included. Intended use for applying mean subtraction to `x` and `y`.")
    end

    if demean
        zy::Vector{Float64} = y .- mean(y)
    end

    #for k = 1 : m  # foreach lag value
        #l = lags[k]
        l=lag

        if l >= 0
           sub_x = zx[1:lx-l]
           sub_y = zy[1+l:lx]
        else
           sub_x = zx[1-l:lx]
           sub_y = zy[1:lx+l]
        end

        #push!(r, correlation_function(sub_x, sub_y))
        r = correlation_function(sub_x, sub_y)
    #end

   return r
end

function pearsons_cor(x, y)
    sc = sqrt(dot(x, x) * dot(y, y))
    r = dot(x, y) / sc

    return r
end

function CE(x)
    return sqrt(sum(diff(x).^2))
end

function CF(x, y)
    return max(CE(x), CE(y)) / min(CE(x), CE(y))
end

function CID(x, y)
    return (sqrt(sum((x - y) .^2)) * CF(x, y))
end

function carbonate_threshold_cover(area, cover; threshold=0.2)
    carbonate_cover = area * threshold

    return (cover .- carbonate_cover) ./ carbonate_cover .* 100
end

function threshold_cover_timeseries(areas, cover_timeseries, threshold)
    threshold_cover = YAXArray(
        (cover_timeseries.timesteps, cover_timeseries.locations),
        zeros(size(cover_timeseries))
    )

    for loc in eachindex(areas)
        threshold_cover[:, loc] = carbonate_threshold_cover(
            areas[loc], 
            cover_timeseries[:, loc]; 
            threshold=threshold
        )
    end

    return threshold_cover
end
