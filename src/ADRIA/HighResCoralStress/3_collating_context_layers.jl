"""
Access the domain connectivity and calculate the strength and number of outgoing/incoming
connections for each reef. Add context information such as reef DHW levels, years above carbonate
budget thresholds, DHW tolerance levels and number of reefs in each bioregion.
"""

using ArchGDAL

include("../../common.jl")

# Change ADRIA debug mode to true to extract DHW tolerance data from runs
# (removes parallel processing)
change_ADRIA_debug(true)


using ADRIA

# Load context layers with bioregions and target reefs
context_layers = GDF.read(joinpath(output_path, "clustered_reefs_carbonate.gpkg"))

# GBR wide domain
gbr_dom = ADRIA.load_domain(gbr_domain_path, "45")
dhw_scenarios = open_dataset(joinpath(gbr_domain_path, "DHWs/dhwRCP45.nc"))
GCMs = dhw_scenarios.dhw.properties["members"]
gbr_dom.loc_data.geometry = Vector{ArchGDAL.IGeometry}(gbr_dom.loc_data.geometry) # Need to recast gbr_dom geometry col for ADRIA.run_scenarios(). Possibly issue with GeoDataFrames version.
gbr_dom_filtered = gbr_dom.loc_data[gbr_dom.loc_data.UNIQUE_ID.∈[context_layers.UNIQUE_ID], :]
filtered_indices = indexin(gbr_dom_filtered.UNIQUE_ID, gbr_dom.loc_data.UNIQUE_ID)

# 1. Attach connectivity data to context_layers
connectivity_matrix = gbr_dom.conn

# Calculate connectivity metrics for context analysis
reefs = collect(getAxis("Source", connectivity_matrix).val)

col_types = [
    String,
    Float64,
    Int64,
    Float64,
    Float64,
    Int64,
    Float64,
    Float64,
    Int64,
    Float64,
    Float64
]

source_to_sink = DataFrame(
    [T[] for T in col_types],
    [
        :UNIQUE_ID,
        :income_strength,
        :income_count,
        :income_comb,
        :out_strength,
        :out_count,
        :out_comb,
        :total_strength,
        :total_count,
        :total_comb,
        :so_to_si
    ]
)

for GCM in GCMs
    context_layers[!, "$(GCM)_weighted_incoming_conn"] .= 1*10^-5
end

for reef in eachindex(reefs)
    reef_id = reefs[reef]

    if !(reef_id ∈ context_layers.UNIQUE_ID)
        continue
    end

    outgoing = connectivity_matrix[connectivity_matrix.Source.==reef_id, :]
    incoming = connectivity_matrix[:, connectivity_matrix.Sink.==reef_id]

    income_strength = sum(incoming)
    income_count = count(incoming .> 0)
    income_comb = income_strength / income_count
    income_comb = isnan(income_comb) ? 0.0 : income_comb

    out_strength = sum(outgoing)
    out_count = count(outgoing .> 0)
    out_comb = out_strength / out_count
    out_comb = isnan(out_comb) ? 0.0 : out_comb

    total_strength = income_strength + out_strength
    total_count = income_count + out_count
    total_comb = total_strength / total_count

    if (income_comb == 0) & (out_comb > 0)
        so_to_si = out_comb
    else
        so_to_si = out_comb / income_comb
    end

    # Get all incoming connectivity indices
    incoming_positive = findall(incoming[:, 1] .> 0)

    if !isempty(incoming_positive)
        for GCM in GCMs
            # Load all source reefs for each GCM
            absolute_cover = readcubedata(open_dataset(joinpath(output_path, "processed_model_outputs/median_cover_$(GCM).nc")).layer)
            absolute_cover = absolute_cover[locations=At(String.(incoming.Source[incoming_positive]))]
            
            weighted_conn_sources = zeros(length(incoming_positive))
            for (i, in_conn) in enumerate(incoming_positive)
                weighted_conn = incoming[in_conn, 1] * mean(absolute_cover[:, i].data)
                weighted_conn_sources[i] = weighted_conn
            end

            context_layers[context_layers.UNIQUE_ID .== reef_id, "$(GCM)_weighted_incoming_conn"] .= sum(weighted_conn_sources)
        end
    end

    push!(
        source_to_sink,
        [
            reef_id,
            income_strength,
            income_count,
            income_comb,
            out_strength,
            out_count,
            out_comb,
            total_strength,
            total_count,
            total_comb,
            so_to_si
        ]
    )
end

context_layers = leftjoin(context_layers, source_to_sink; on=:UNIQUE_ID, order=:left)

context_layers.initial_coral_cover = vec(
    sum(gbr_dom.init_coral_cover[locations=At(context_layers.UNIQUE_ID)]; dims=1).data' .*
    gbr_dom_filtered.area .*
    gbr_dom_filtered.k
)
context_layers.initial_proportion = (
    context_layers.initial_coral_cover ./
    (context_layers.area .* context_layers.k)
)

areas = gbr_dom_filtered.area
thresholds = 10:1:20

# Attaching GCM-dependent context layers
for (i_gcm, GCM) in enumerate(GCMs)

    # 3. Calculate number of years each reef is above a carbonate budget threshold
    absolute_cover = readcubedata(open_dataset(joinpath(output_path, "processed_model_outputs/median_cover_$(GCM).nc")).layer)
    absolute_cover = absolute_cover[locations=At(context_layers.UNIQUE_ID)]
    for threshold in thresholds
        t_threshold = threshold / 100
        reef_thresholds = context_layers.area .* t_threshold

        n_years_above = [sum(absolute_cover[:, loc] .> reef_thresholds[loc]) for loc in eachindex(reef_thresholds)]
        context_layers[!, "$(GCM)_years_above_$(threshold)"] = n_years_above
    end

    # 3. Attach DHW data to context_layers
    # mean DHW data for sites

    dhw_time = gbr_dom.dhw_scens[:, filtered_indices, i_gcm]
    dhw_locs = Float64.(mapslices(mean, dhw_time[1:30, :], dims=[:timesteps]))

    context_layers[:, "$(GCM)_mean_dhw"] = dhw_locs

    threshold_cover = percentage_cover_timeseries(context_layers.area, absolute_cover)[1:30, :]
    negative_normal_dhw = -normalise(dhw_time[:, :], (0, 1))

    dhw_cover_cor = zeros(Float64, length(context_layers.RME_UNIQUE_ID))

    for ind in 1:length(context_layers.UNIQUE_ID)
        reef_total_cover = threshold_cover[locations=ind]
        reef_dhw = negative_normal_dhw[sites=ind]

        dhw_cover_cor[ind] = cross_correlation(reef_total_cover.data, reef_dhw.data, 0, pearsons_cor, true)
    end

    context_layers[!, "$(GCM)_dhw_cover_cor"] = dhw_cover_cor

    # 4. Re-run 5 scenarios for each GCM and extract mean DHW tolerance from timeseries
    scens = CSV.read(joinpath(output_path, "HighResCoralStress_ADRIA_scens_$(GCM).csv"), DataFrame)
    scens = scens[1:4, :]
    scens.wave_scenario .= 1.0

    rs_dhw = ADRIA.run_scenarios(gbr_dom, scens, "45")
    dhw_tol_log = Float64.(mapslices(median, rs_dhw.coral_dhw_tol_log, dims=[:scenarios, :species]))

    dhw_tol_mean = dropdims(mean(dhw_tol_log[10:end, :], dims=1), dims=1).data
    context_layers[!, "$(GCM)_mean_DHW_tol"] = dhw_tol_mean[filtered_indices]

    context_layers[:, "$(GCM)_DHW_tol_rate"] .= 0.0
    for (r, reef) in enumerate(eachrow(context_layers))
        tol_rate = (maximum(dhw_tol_log[:, r].data) - minimum(dhw_tol_log[:, r].data)) / (size(dhw_tol_log, 1) / 10)
        reef["$(GCM)_DHW_tol_rate"] = tol_rate
    end
end

# 5. Calculate the number of reefs in each bioregion
context_layers.bioregion_count .= 1
for reef in eachrow(context_layers)
    reef.bioregion_count = count(context_layers.bioregion .== [reef.bioregion])
end

# 6. Add average latitude for each bioregion and management area
bioregion_average_latitude = combine(groupby(context_layers, :bioregion)) do sdf
    return average_latitude = mean([ArchGDAL.gety(ArchGDAL.centroid(geom), 0) for geom in sdf.geometry])
end
rename!(bioregion_average_latitude, :x1 => :bioregion_average_latitude)
bioregion_average_latitude = Dict(tuple.(bioregion_average_latitude.bioregion, bioregion_average_latitude.bioregion_average_latitude))
context_layers.bioregion_average_latitude = [bioregion_average_latitude[bior] for bior in context_layers.bioregion]

management_area_average_latitude = combine(groupby(context_layers, :management_area)) do sdf
    return average_latitude = mean([ArchGDAL.gety(ArchGDAL.centroid(geom), 0) for geom in sdf.geometry])
end
rename!(management_area_average_latitude, :x1 => :management_area_average_latitude)
management_area_average_latitude = Dict(tuple.(management_area_average_latitude.management_area, management_area_average_latitude.management_area_average_latitude))
context_layers.management_area_average_latitude = [management_area_average_latitude[bior] for bior in context_layers.management_area]

# 7. Check for reef cluster consistency across GCMs at management area scale
man_area_gcm_cols = ["$(GCM)_management_area_cluster_cats" for GCM in GCMs]
context_layers.man_area_consistent_reefs = Vector{Any}([allequal(row) for row in eachrow(context_layers[:, man_area_gcm_cols])])

for reef in eachrow(context_layers)
    if reef.man_area_consistent_reefs
        reef.man_area_consistent_reefs = unique(reef[man_area_gcm_cols])[1]
    else
        reef.man_area_consistent_reefs = "variable"
    end
end

# Format columns for writing to geopackage
context_layers.income_strength = Float64.(context_layers.income_strength)
context_layers.income_count .= convert.(Int64, context_layers.income_count)
context_layers.income_comb .= convert.(Float64, context_layers.income_comb)
context_layers.out_strength .= convert.(Float64, context_layers.out_strength)
context_layers.out_count .= convert.(Int64, context_layers.out_count)
context_layers.out_comb .= convert.(Float64, context_layers.out_comb)
context_layers.total_strength .= convert.(Float64, context_layers.total_strength)
context_layers.total_count .= convert.(Int64, context_layers.total_count)
context_layers.total_comb .= convert.(Float64, context_layers.total_comb)
context_layers.so_to_si .= convert.(Float64, context_layers.so_to_si)
context_layers.man_area_consistent_reefs .= convert.(String, context_layers.man_area_consistent_reefs)

GDF.write(joinpath(output_path, "analysis_context_layers_carbonate.gpkg"), context_layers; crs=GFT.EPSG(7844), overwrite=true)

change_ADRIA_debug(false) # reset ADRIA debug mode to false as required by other scripts
