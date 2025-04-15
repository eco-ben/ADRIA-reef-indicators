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

# """
#     mapslices_toFloat64(func, data::YAXArray, dim::Symbol)

# Apply a function across grouped slices of a YAXArray.
# E.g. To reduce from 3 to 2 dimensions, use a summarising `func` and the dimension to be
# reduced as `dim`.
# E.g. To apply function over timesteps for each reef, use a non-summarising `func` and
# timesteps as `dim`.

# # Arguments
# - `func` : Function to apply to array slices
# - `data` : YAXArray for slicing (2D or 3D)
# - `dim` : Target dimension symbol (if reducing use unwanted dim, if transforming across each series use timesteps)

# # Returns
# YAXArray after applying function for use in further ADRIA analysis/viz functions.
# """
# function mapslices_toFloat64(func, data::YAXArray, dim::Symbol)
#     init_array = mapslices(func, data, dims=[dim])
#     dat = convert.(Float64, init_array.data)
#     new_array::YAXArray{Float64} = YAXArray(dims(init_array), dat)

#     return new_array
# end

"""
    cross_correlation(
        x::AbstractVector{<:Real},
        y::AbstractVector{<:Real},
        lags::AbstractVector{<:Integer}
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

# Returns
Vector of correlation values for each lag in `lags`.
"""
function cross_correlation(
    x::AbstractVector{<:Real},
    y::AbstractVector{<:Real},
    lags::AbstractVector{<:Integer}
    )

    r = Vector{Float64}()
    lx = length(x)
    m = length(lags)

    zx::Vector{Float64} = x
    zy::Vector{Float64} = y

    for k = 1 : m  # foreach lag value
        l = lags[k]

        if l >= 0
           sub_x = zx[1:lx-l]
           sub_y = zy[1+l:lx]
        else
           sub_x = zx[1-l:lx]
           sub_y = zy[1:lx+l]
        end

        sc = sqrt(dot(sub_x, sub_x) * dot(sub_y, sub_y))

        push!(r, dot(sub_x, sub_y) / sc)
    end

   return r
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

# """
#     lagged_cluster_analysis(
#         region_rel_cover::YAXArray,
#         region_clusters::Vector{Int64},
#         lags::AbstractVector{<:Integer}
#     )::DataFrame

# Perform lagged cross correlation analysis across a number of clusters for a target
# region. Uses `mapslices_toFloat64()` and `cross_correlation()` functions.

# # Arguments
# - `region_rel_cover` : YAXArray of cover trajectories for each reef in the target region.
# - `region_clusters` : Vector of cluster categories assigned to each reef.
# - `lags` : Vector of lags to apply with `cross_correlation()`. Positive lags test for `x` leading `y`, negative lags test for `y` leading `x`.

# # Returns
# DataFrame of correlation values for each reef in each cluster with columns for the reef RME_UNIQUE_ID, cluster category and correlation values for each lag step.
# """
# function lagged_cluster_analysis(
#         region_rel_cover::YAXArray,
#         region_clusters::Vector{Int64},
#         lags::AbstractVector{<:Integer}
#     )::DataFrame

#     lags_symbols = [Symbol("lag" * string(lags[a])) for a in eachindex(lags)]
#     cross_cor = DataFrame(
#         [Vector{Any}() for _ in 1:(length(lags)+2)],
#         [:RME_UNIQUE_ID; :t_cluster ; lags_symbols]
#     )

#     for cluster in 1:n_clusters
#         target_cluster = region_rel_cover[:, (region_clusters .== cluster)]

#         cluster_median = Float64.(mapslices(median, target_cluster, dims=[:sites])) # Have to remove the target reef before

#         for (ind, reef) in enumerate(eachcol(target_cluster))
#             reef_name = target_cluster.sites[ind]
#             correlation = cross_correlation(reef, cluster_median, lags, true)

#             push!(cross_cor, [reef_name; cluster; correlation])
#         end
#     end

#     return cross_cor
# end

# """
#     lagged_region_analysis(
#         region_rel_cover::YAXArray,
#         reg::String,
#         lags::AbstractVector{<:Integer}
#     )::DataFrame

# Perform lagged cross correlation analysis across a a target region.
# Uses `mapslices_toFloat64()` and `cross_correlation()` functions.

# # Arguments
# - `region_rel_cover` : YAXArray of cover trajectories for each reef in the target region.
# - `reg` : String of region name.
# - `lags` : Vector of lags to apply with `cross_correlation()`. Positive lags test for `x` leading `y`, negative lags test for `y` leading `x`.

# # Returns
# DataFrame of correlation values for each reef in target_region with columns for the reef RME_UNIQUE_ID, region and correlation values for each lag step.
# """
# function lagged_region_analysis(
#         region_rel_cover::YAXArray,
#         reg::String,
#         lags::AbstractVector{<:Integer}
#     )::DataFrame

#     lags_symbols = [Symbol("lag" * string(lags[a])) for a in eachindex(lags)]
#     cross_cor = DataFrame(
#         [Vector{Any}() for _ in 1:(length(lags)+2)],
#         [:RME_UNIQUE_ID; :region; lags_symbols]
#     )

#     reg_median = Float64.(mapslices(median, region_rel_cover, dims=[:sites]))

#     for (ind, reef) in enumerate(eachcol(region_rel_cover))
#         reef_name = region_rel_cover.sites[ind]
#         correlation = cross_correlation(reef, reg_median, lags, true)

#         push!(cross_cor, [reef_name; reg; correlation])
#     end

#     return cross_cor
# end

# """
#     subregion_analysis(
#         subregions::Vector{String},
#         rel_cover::YAXArray,
#         context_layers::DataFrame,
#         category::Symbol
#     )::DataFrame

# Apply lagged_region_analysis across a list of subregions such as closest_ports or bioregions.
# Uses `lagged_region_analysis()` function, lags must be specified.

# # Arguments
# - `subregions` : Vector of subregion names to apply analysis to. Must be found in `context_layers.category` column.
# - `rel_cover` : YAXArray of cover trajectories for each reef.
# - `context_layers` : DataFrame containing columns `RME_UNIQUE_ID` and category column matching `subregions`.
# - `category` : Column name in `context_layers` that contains `subregions` for each reef.
# - `lags` : Vector of lags to apply with `cross_correlation()`. Positive lags test for `x` leading `y`, negative lags test for `y` leading `x`.

# # Returns
# DataFrame of correlation values for each reef and surrounding reefs in 'subregion'. Contains RME_UNIQUE_ID, region and correlation values for each lag step.
# """
# function subregion_analysis(
#     subregions::Vector{String},
#     rel_cover::YAXArray,
#     context_layers::DataFrame,
#     category::Symbol,
#     lags::AbstractVector{<:Integer}
#     )::DataFrame

#     lagged_analysis_sub = DataFrame()
#     for subregion in subregions
#         subregion_reefs = context_layers[(context_layers[:, category] .== subregion), :RME_UNIQUE_ID]
#         subregion_cover = rel_cover[:, (findall(rel_cover.sites .∈ [subregion_reefs]))]

#         subregion_lagged_analysis = lagged_region_analysis(subregion_cover, subregion, lags)
#         lagged_analysis_sub = vcat(lagged_analysis_sub, subregion_lagged_analysis)
#     end

#     return lagged_analysis_sub
# end

"""
    cluster_correlation(reef_clusters, timeseries, lag)

Calculate the correlation between each reef timeseries and the median timeseries for other reefs
in the same cluster. Clusters are defined with `reef`
"""
function cluster_correlation(reef_clusters, timeseries, lag, func)
    correlation_values = Vector{Union{Missing, Float64}}(missing, length(reef_clusters))

    for (ind, cluster) in enumerate(reef_clusters)
        indices_not_target = findall(reef_clusters .== cluster)
        indices_not_target = indices_not_target[indices_not_target .!= ind]

        target_timeseries = timeseries[locations = ind]
        non_target_median = timeseries[locations = indices_not_target]
        non_target_median = vec(mapslices(median, non_target_median.data, dims=2))
        correlation_values[ind] = cross_correlation(target_timeseries, non_target_median, lag, func, true)
    end

    return correlation_values
end


canonical_reefs = find_latest_file("../../canonical-reefs/output/")
canonical_reefs = GDF.read(canonical_reefs)

#MANAGEMENT_AREAS = unique(canonical_reefs.management_area)

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

function prepare_ReefMod_results(
    result_store_dir::String,
    location_ids,
    start_year,
    end_year,
    fn,
    evenness_method;
    evenness_weight=1,
    cover_weight=1
)
    results = open_dataset("$(result_store_dir)/results_no_duplicates.nc")
    total_cover = results.total_cover
    taxa_cover = results.total_taxa_cover

    cover_median = Float64.(mapslices(median, total_cover, dims=[:scenarios]))
    taxa_median = Float64.(mapslices(median, taxa_cover, dims=[:scenarios]))
    taxa_evenness = _coral_evenness(taxa_median; method=evenness_method, evenness_weight=evenness_weight, cover_weight=cover_weight)

    axlist = (
                Dim{:timesteps}(start_year:end_year),
                Dim{:sites}(location_ids)
    )
    cover_median = rebuild(cover_median, dims=axlist)
    taxa_evenness = rebuild(taxa_evenness, dims=axlist)
    arrays = Dict(
        :total_relative_cover_median => cover_median,
        :scaled_taxa_evenness => taxa_evenness
    )
    ds = Dataset(; arrays...)
    savedataset(ds, path="$(result_store_dir)/$(fn)", driver=:netcdf, overwrite=true)

    return nothing
end

function prepare_ReefMod_results(
    results::YAXArrays.Dataset,
    result_store_dir,
    location_ids,
    start_year,
    end_year,
    fn,
    evenness_method;
    evenness_weight=1,
    cover_weight=1
)
    total_cover = results.total_cover
    taxa_cover = results.total_taxa_cover

    cover_median = Float64.(mapslices(median, total_cover, dims=[:scenarios]))
    taxa_median = Float64.(mapslices(median, taxa_cover, dims=[:scenarios]))
    cyc_median = Float64.(mapslices(median, results.cyc_mortality; dims=[:scenarios]))
    cots_median = Float64.(mapslices(median, results.cots_mortality; dims=[:scenarios]))
    dhw_mortality_median = Float64.(mapslices(median, results.dhw_mortality; dims=[:scenarios]))
    taxa_evenness = _coral_evenness(taxa_median; method=evenness_method, evenness_weight=evenness_weight, cover_weight=cover_weight)

    dhw_median = Float64.(mapslices(median, results.dhw, dims=[:scenarios]))
    negative_normal_dhw = -normalise(dhw_median.data, (0,1))

    axlist = (
                Dim{:timesteps}(start_year:end_year),
                Dim{:sites}(location_ids)
    )
    cover_median = rebuild(cover_median, dims=axlist)
    taxa_cover = rebuild(
        taxa_median,
        dims=(
            Dim{:timesteps}(start_year:end_year),
            Dim{:sites}(location_ids),
            Dim{:species}(1:6)
        )
    )
    taxa_evenness = rebuild(taxa_evenness, dims=axlist)
    cyc_mortality = rebuild(cyc_median, dims=axlist)
    cots_mortality = rebuild(cots_median, dims=axlist)
    dhw_mortality = rebuild(dhw_mortality_median, dims=axlist)
    dhw = rebuild(dhw_median, dims=axlist)
    negative_normal_dhw = YAXArray(axlist, negative_normal_dhw)
    arrays = Dict(
        :total_relative_cover_median => cover_median,
        :taxa_cover => taxa_cover,
        :scaled_taxa_evenness => taxa_evenness,
        :dhw_mortality => dhw_mortality,
        :cots_mortality => cots_mortality,
        :cyc_mortality => cyc_mortality,
        :dhw => dhw,
        :negative_normal_dhw => negative_normal_dhw
    )
    ds = Dataset(; arrays...)
    savedataset(ds, path="$(result_store_dir)/$(fn)", driver=:netcdf, overwrite=true)

    return nothing
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
    remove_duplicate_reps(rs_dataset, start_year, end_year, location_ids, n_reps)

Find the indices of unique scenarios when there are duplicated scenarios and rebuild
the scenarios axis in `rebuild_RME_dataset()` to contain only a single copy of unique scenarios.
"""
function remove_duplicate_reps(rs_dataset, start_year, end_year, location_ids, n_reps)
    cover = rs_dataset.total_cover

    for year_reef1 in cover.timesteps
        cover_scen = cover[At(year_reef1),1,:]
        if size(unique(cover_scen.data), 1) == n_reps
            global unique_indices = unique(i -> cover_scen.data[i], 1:length(cover_scen.data))
            @info "200 unique reps found."
            break
        end
    end

    rs_dataset = rebuild_RME_dataset(
        rs_dataset,
        start_year, end_year,
        location_ids,
        n_reps,
        unique_indices
    )

    return rs_dataset
end

"""
    rebuild_RME_dataset(
        rs_dataset,
        start_year,
        end_year,
        location_ids,
        n_reps,
        unique_indices
    )

Rebuild a RME dataset that has duplicated scenarios. For example, when RME outputs counterfactual runs with duplicate scenario data.

# Arguments
- `rs_dataset` : The RME dataset with duplicated scenarios.
- `start_year` : Start year of timesteps dimension.
- `end_year` : End year of timesteps dimension.
- `location_ids` : Location IDs to be held in sites dimension.
- `n_reps` : The intended number of scenarios that should be in the returned dataset (after removing duplicate scenarios).
- `unique_indices` : The first index of each unique scenario to keep (excludes indices of duplicate scenarios).
"""
function rebuild_RME_dataset(
    rs_dataset,
    start_year,
    end_year,
    location_ids,
    n_reps,
    unique_indices
)
    variable_keys = keys(rs_dataset.cubes)

    arrays = Dict()
    for variable in variable_keys
        if variable == :total_taxa_cover
            axlist = (
                Dim{:timesteps}(start_year:end_year),
                Dim{:sites}(location_ids),
                Dim{:taxa}(1:6),
                Dim{:scenarios}(1:n_reps)
            )
        else
            axlist = (
                Dim{:timesteps}(start_year:end_year),
                Dim{:sites}(location_ids),
                Dim{:scenarios}(1:n_reps)
            )
        end

        # Remove duplicated scenarios
        yarray = rs_dataset[variable][scenarios = unique_indices]
        # Rebuild to ensure correct scenario lookup axis.
        yarray = rebuild(yarray, dims=axlist)
        push!(arrays, variable => yarray)
    end

    return Dataset(; arrays...)
end

"""
    concat_RME_netcdfs(dataset_1, dataset_s...)

Combine RME result netcdf datasets along the `scenarios` dimension to
combine scenarios that have been run separately into a single dataset.

# Example
results_dataset_300scens = concat_RME_netcdfs(
    results_dataset_200scens,
    results_dataset_50scens,
    results_dataset_50scens
)
"""
function concat_RME_netcdfs(dataset_1, dataset_s...)
    variable_keys = keys(dataset_1.cubes)
    Xin = [dataset_1, dataset_s...]
    arrays = Dict()

    for variable in variable_keys
        if variable == :total_taxa_cover
            yarrays = [x[variable] for x in Xin]
            yarray = YAXArrays.cat(yarrays...; dims=4) # In RME YAXArrays with taxa the 4th dimension is scenarios

            # For some reason after concattenating you need to rebuild the scenario axis
            axlist = (
                yarray.axes[1],
                yarray.axes[2],
                yarray.axes[3],
                Dim{:scenarios}(1:size(yarray,4))
            )
            yarray = rebuild(yarray, dims=axlist)
        else
            yarrays = [x[variable] for x in Xin]
            yarray = YAXArrays.cat(yarrays...; dims=3) # In RME YAXArrays without taxa the 3rd dimension is scenarios

            # For some reason after concattenating you need to rebuild the scenario axis
            axlist = (
                yarray.axes[1],
                yarray.axes[2],
                Dim{:scenarios}(1:size(yarray,3))
            )
            yarray = rebuild(yarray, dims=axlist)
        end

        push!(arrays, variable => yarray)
    end

    return Dataset(; arrays...)
end

"""
    load_result_store(dir_name::String, n_reps::Int64)::ResultStore

Save ResultStore from saved results.nc and scenarios.csv files to allow modification.
"""
function load_result_store(dir_name::String, n_reps::Int64)::ResultStore
    result_path = joinpath(dir_name, "results.nc")
    results = open_dataset(result_path, driver=:netcdf)
    start_year = first(results.timesteps)
    end_year = last(results.timesteps)
    n_reefs = length(results.locations)

    scenario_path = joinpath(dir_name, "scenarios.csv")
    scenario = CSV.read(scenario_path, DataFrame)

    return ResultStore(
        results,
        scenario,
        start_year,
        end_year,
        (end_year - start_year) + 1,
        n_reefs,
        n_reps
    )
end

function stat_range(x)
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

function bellwether_glmm_fit(dataset)
    y_term_cover = term("$(GCM)_target_reefs_bior")
    mean_dhw_term = term("$(GCM)_mean_dhw")
    icc_term = term("$(GCM)_initial_coral_cover")
    #dhw_cover_cor_term = term("$(GCM)_dhw_cover_cor")
    cots_mortality_term = term("$(GCM)_mean_cots_mortality")

    glmm_form =
        y_term_cover ~
        mean_dhw_term +
        term(:so_to_si) +
        term(:conn_score) +
        term(:total_strength) +
        icc_term +
        #dhw_cover_cor_term +
        cots_mortality_term +
        (term(1)|term(:bioregion))
    glmm_fit = fit(MixedModel, glmm_form, dataset, Bernoulli())
    coef_table = coeftable(glmm_fit)
    results = DataFrame(
        Names=coef_table.rownms,
        Coef=coef_table.cols[1],
        std_error=coef_table.cols[2],
        z=coef_table.cols[3],
        p=coef_table.cols[4]
    )

    return results
end

function grouping_counts(grouping_col::Symbol, dataset, bellwether_reefs_col, fn, GCM)
    reef_counts = combine(groupby(dataset, grouping_col)) do sdf
        return (bellwether_reefs = sum(sdf[:, bellwether_reefs_col]), total_reefs=nrow(sdf))
    end

    removed_groups = reef_counts[reef_counts[:, :bellwether_reefs] .< 5, grouping_col]

    CSV.write(fn, reef_counts)

    return removed_groups
end

function grouping_counts(grouping_col::Symbol, dataset, clustering_col, n_clusters::Int64)
    reef_counts = combine(groupby(dataset, grouping_col)) do sdf
        return (clusters = length(unique(sdf[:, clustering_col])), total_reefs=nrow(sdf))
    end

    removed_groups = reef_counts[reef_counts[:, :clusters] .< n_clusters, grouping_col]


    return removed_groups
end

function naive_split_metric(obs::AbstractVector, n_members::Int, metric::Function=mean)
    obs_chunks = Iterators.partition(obs, n_members)
    scores = Vector{Float64}(undef, ceil(Int64, length(obs) / n_members))

    for (idx, chunk) in enumerate(obs_chunks)
        scores[idx] = metric(chunk)
    end

    return scores
end

# """
#     temporal_variability(x::AbstractVector{<:Real})
#     temporal_variability(x::AbstractArray{<:Real, 2})
#     temporal_variability(x::AbstractArray{<:Real}, func_or_data...)

# The V meta-metric.

# As a meta-metric, it can be applied to any combination of
# metrics (including itself), assuming \$x\$ is bound between 0 and 1.
# If this is not the case, consider normalizing values first.

# # Examples
# ```julia-repl
# # Apply V to a time series
# julia> temporal_variability(rand(50))

# # Apply V to an ensemble of time series
# julia> x = rand(50, 200)
# julia> temporal_variability(x)

# # Create and apply a modified V metric to an ensemble of time series.
# # Where the argument is an array and not a function, the data is used directly
# # and so it is assumed all matrices are of the same size and shape.
# julia> temporal_variability(x, temporal_variabilty, temporal_variability(P(x)))
# julia> temporal_variability(x, temporal_variabilty, P(x), D(x), E(x))
# ```
# """
# function temporal_variability(x::AbstractVector{<:Real}; w=[0.9, 0.1])
#     return mean([mean(x)*w[1], std(x)*w[2]])
# end


# """

# """
# function cluster_temporal_variability(reef_clusters, temporal_variability_results, threshold)
#     bellwether_reefs = zeros(Bool, length(reef_clusters))

#     for (ind, cluster) in enumerate(reef_clusters)
#         indices_not_target = findall(reef_clusters .== cluster)
#         indices_not_target = indices_not_target[indices_not_target .!= ind]

#         target_temporal_variability = temporal_variability_results[ind]
#         non_target_temp_variability = median(temporal_variability_results[indices_not_target])

#         bellwether_reefs[ind] = (non_target_temp_variability - target_temporal_variability) > threshold
#     end

#     return bellwether_reefs
# end

# function cluster_temporal_variability(reef_clusters, temporal_variability_results; quantile_threshold=0.25)
#     bellwether_reefs = zeros(Bool, length(reef_clusters))

#     for (ind, cluster) in enumerate(reef_clusters)
#         indices_not_target = findall(reef_clusters .== cluster)
#         indices_not_target = indices_not_target[indices_not_target .!= ind]

#         target_temporal_variability = temporal_variability_results[ind]
#         non_target_quantile = quantile(temporal_variability_results[indices_not_target], quantile_threshold)

#         bellwether_reefs[ind] = target_temporal_variability < non_target_quantile
#     end

#     return bellwether_reefs
# end

function GCM_analysis_results(results_set)
    
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

function grouped_timeseries_clustering(timeseries, grouping, n_clusters; length_t=1:size(timeseries, 1))
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

function find_dict_val(val, dict)
    for (k,v) in dict
        if v == val
            return k
        end
    end
end

# function find_substantial_difference(values, groups)
#     unique_groups = unique(groups)

#     medians = [median(values[groups .== group]) for group in unique_groups]
#     quantile_75 = [quantile(values[groups .== group], 0.75) for group in unique_groups]
#     quantile_25 = [quantile(values[groups .== group], 0.25) for group in unique_groups]

