"""
Access the domain connectivity and calculate the strength and number of outgoing/incoming
connections for each reef. Add context information such as reef DHW levels, years above carbonate
budget thresholds, DHW tolerance levels and number of reefs in each bioregion.
"""

using ArchGDAL

change_ADRIA_debug(true) # Change ADRIA debug mode to true to extract DHW tolerance data from runs

include("../../common.jl")

# Load context layers with bioregions and target reefs
context_layers = GDF.read("../outputs/ADRIA_results/HighResCoralStress/bellwether_reefs_carbonate.gpkg")

# GBR wide domain
gbr_dom = ADRIA.load_domain("../../ADRIA Domains/GBR_2024_10_15_HighResCoralStress/", "45")
dhw_scenarios = open_dataset("../../ADRIA Domains/GBR_2024_10_15_HighResCoralStress/DHWs/dhwRCP45.nc")

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
    Float64,
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
        :so_to_si,
    ]
)

for reef in eachindex(reefs)
    reef_id = reefs[reef]
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
    sum(gbr_dom.init_coral_cover; dims=1).data' .*
    gbr_dom.loc_data.area .*
    gbr_dom.loc_data.k
)
context_layers.initial_proportion = (
    context_layers.initial_coral_cover ./
    (context_layers.area .* context_layers.k)
)

thresholds = 10:1:20
# Attaching GCM-dependent context layers
for (i_gcm, GCM) in enumerate(dhw_scenarios.dhw.properties["members"])

    rs = ADRIA.load_results("../outputs/ADRIA_results/HighResCoralStress/GBR_2024_10_15_HighResCoralStress__RCPs_45_$(GCM)")
    GCM_results = GCM_analysis_results(rs)

    # 3. Calculate number of years each reef is above a carbonate budget threshold
    @floop for threshold in thresholds
        t_threshold = threshold / 100
        reef_thresholds = context_layers.area .* t_threshold

        absolute_cover = GCM_results.absolute_median_cover
        n_years_above = [sum(absolute_cover[:, loc] .> reef_thresholds[loc]) for loc in eachindex(reef_thresholds)]
        context_layers[!, "$(GCM)_years_above_$(threshold)"] = n_years_above
    end

    # 3. Attach DHW data to context_layers
    # mean DHW data for sites

    dhw_time = gbr_dom.dhw_scens[:, :, i_gcm]
    dhw_locs = Float64.(mapslices(mean, dhw_time[1:30, :], dims=[:timesteps]))

    context_layers[:, "$(GCM)_mean_dhw"] = dhw_locs

    rel_cover = GCM_results.relative_cover[1:30, :]
    negative_normal_dhw = -normalise(dhw_time[:, :], (0, 1))

    dhw_cover_cor = zeros(Float64, length(context_layers.RME_UNIQUE_ID))

    for ind in 1:length(context_layers.UNIQUE_ID)
        reef_total_cover = rel_cover[locations=ind]
        reef_dhw = negative_normal_dhw[sites=ind]

        dhw_cover_cor[ind] = cross_correlation(reef_total_cover.data, reef_dhw.data, 0, pearsons_cor, true)

    end

    context_layers[!, "$(GCM)_dhw_cover_cor"] = dhw_cover_cor

    # 4. Re-run 5 scenarios for each GCM and extract mean DHW tolerance from timeseries
    # scens = rs.inputs[1:4, Not(318)] # Only use 1st 5 scenarios for speed and remove RCP variable in inputs
    # scens.wave_scenario .= 1.0
    # rs_dhw = ADRIA.run_scenarios(gbr_dom, scens, "45")
    # dhw_tol_log = Float64.(mapslices(median, rs_dhw.coral_dhw_tol_log, dims = [:scenarios, :species]))
    # dhw_tol_mean = dropdims(mean(dhw_tol_log[10:end, :], dims = 1), dims = 1).data

    # context_layers[!, "$(GCM)_mean_DHW_tol"] = dhw_tol_mean
end

# 5. Calculate the number of reefs in each bioregion
context_layers.bioregion_count .= 1
for reef in eachrow(context_layers)
    reef.bioregion_count = count(context_layers.bioregion .== [reef.bioregion])
end

# 6. Add average latitude for each bioregion
bioregion_average_latitude = combine(groupby(context_layers, :bioregion)) do sdf
    return average_latitude = mean([ArchGDAL.gety(ArchGDAL.centroid(geom), 0) for geom in sdf.geometry])
end
rename!(bioregion_average_latitude, :x1 => :bioregion_average_latitude)
bioregion_average_latitude = Dict(tuple.(bioregion_average_latitude.bioregion, bioregion_average_latitude.bioregion_average_latitude))
context_layers.bioregion_average_latitude = [bioregion_average_latitude[bior] for bior in context_layers.bioregion]

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

GDF.write("../outputs/ADRIA_results/HighResCoralStress/analysis_context_layers_carbonate.gpkg", context_layers; crs=GFT.EPSG(7844), overwrite=true)

change_ADRIA_debug(false) # reset ADRIA debug mode to false as required by other scripts