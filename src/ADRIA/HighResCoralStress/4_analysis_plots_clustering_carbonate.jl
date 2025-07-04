"""
Create analysis plots for clustering and carbonate budget analyses. Performed for each GCM separately,
and multi-GCM comparison plot created at end of the script.
"""

using Revise, Infiltrator

using GLMakie, GeoMakie, GraphMakie
using Statistics
using YAXArrays
using CairoMakie

using ADRIA

include("../../common.jl")
include("../../plotting_functions.jl")

CairoMakie.activate!()

# # Select the target GCM from
# dhw_scenarios = open_dataset("../../ADRIA Domains/GBR_2024_10_15_HighResCoralStress/DHWs/dhwRCP45.nc")
# GCMs = dhw_scenarios.dhw.properties["members"]
GCMs = [
    "EC-Earth3-Veg",
    "ACCESS-ESM1-5",
    "ACCESS-CM2",
    "NorESM2-MM",
    "GFDL-CM4"
]

context_layers = GDF.read(joinpath(output_path, "analysis_context_layers_carbonate.gpkg"))
context_layers.gbr .= "Great Barrier Reef"
context_layers.log_so_to_si = log10.(context_layers.so_to_si)
context_layers.log_total_strength = log10.(context_layers.total_strength)

gbr_dom = ADRIA.load_domain(gbr_domain_path, "45")
areas = gbr_dom.loc_data.area

thresholds = 10:1:20
year_cols = Vector{String}()
for GCM in GCMs
    for threshold in thresholds
        push!(year_cols, "$(GCM)_years_above_$(threshold)")
    end
end

reefs_long = stack(
    context_layers[:, ["UNIQUE_ID", "depth_med", "log_total_strength", year_cols...]],
    year_cols
)

for (i_gcm, GCM) in enumerate(GCMs)
    # Select GCM and load relevant results
    @info "analysing reef clustering for $(GCM)"

    fig_out_dir = joinpath(output_path, "figs/$(GCM)")

    absolute_cover = open_dataset(joinpath(output_path, "processed_model_outputs/median_cover_$(GCM).nc")).layer
    # threshold_cover = threshold_cover_timeseries(areas, absolute_cover, 0.17)
    threshold_cover = percentage_cover_timeseries(areas, absolute_cover)

    dhw_ts = gbr_dom.dhw_scens[:,:,i_gcm]
    dhw_ts = rebuild(dhw_ts, dims=threshold_cover.axes)

    # Subset layers to remove reefs that start with very low coral cover
    analysis_layers = context_layers[(context_layers.management_area .!= "NA") .& (context_layers.bioregion .!= "NA"), :]

    # Analysing clusters identified at bioregion level
    cluster_analysis_plots(GCM, analysis_layers, threshold_cover, dhw_ts, :bioregion, fig_out_dir)
    cluster_analysis_plots(GCM, analysis_layers, threshold_cover, dhw_ts, :management_area, fig_out_dir)
    cluster_analysis_plots(GCM, analysis_layers, threshold_cover, dhw_ts, :gbr, fig_out_dir)

    gcm_reefs_long = reefs_long[contains.(reefs_long.variable, GCM), :]
    gcm_reefs_long.variable = last.(split.(gcm_reefs_long.variable, "_"))
    gcm_reefs_long_depth = gcm_reefs_long[
        gcm_reefs_long.UNIQUE_ID .∈ [context_layers[context_layers.depth_qc .== 0, :UNIQUE_ID]],
    :]

    depth_carbonate_scatter = carbonate_budget_variable_scatter(gcm_reefs_long_depth, :depth_med, :value, :variable; color_label = "Median reef depth (m)")
    save(
        joinpath(fig_out_dir, "depth_carbonate_scatter.png"), 
        depth_carbonate_scatter, 
        px_per_unit = dpi
    )
    total_connectivity_carbonate_scatter = carbonate_budget_variable_scatter(gcm_reefs_long, :log_total_strength, :value, :variable; color_label = "Log total connectivity strength")
    save(
        joinpath(fig_out_dir, "log_total_strength_carbonate_scatter.png"), 
        total_connectivity_carbonate_scatter, 
        px_per_unit = dpi
    )

end

man_area_gcm_cluster_cols = [Symbol("$(GCM)_management_area_clusters") for GCM in GCMs]
bioregion_gcm_cluster_cols = [Symbol("$(GCM)_bioregion_clusters") for GCM in GCMs]
gbr_gcm_cluster_cols = [Symbol("$(GCM)_gbr_clusters") for GCM in GCMs]
analysis_layers_long = stack(
    context_layers[:, [:UNIQUE_ID, :management_area, man_area_gcm_cluster_cols...]], 
    man_area_gcm_cluster_cols
)
analysis_layers_long.GCM = [first(split(name, "_")) for name in analysis_layers_long.variable]

rel_cover_arrays = [
    percentage_cover_timeseries(
        areas, 
        open_dataset(joinpath(output_path, "processed_model_outputs/median_cover_$(GCM).nc")).layer
        ) for GCM in GCMs
]

rel_cover_arrays = concatenatecubes(rel_cover_arrays, Dim{:GCM}(GCMs))
analysis_layers_long = analysis_layers_long[analysis_layers_long.management_area .!= "NA", :]
rel_cover_arrays = rel_cover_arrays[locations = (rel_cover_arrays.locations .∈ [unique(analysis_layers_long.UNIQUE_ID)])]

GCM_comparison_plot = grouped_GCM_cluster_timeseries_plots(
    rel_cover_arrays,
    analysis_layers_long,
    :value,
    [:management_area, :GCM],
    1:50
)
save(
    joinpath(output_path, "figs/GCM_timeseries_plot.png"), 
    GCM_comparison_plot,
    px_per_unit = dpi
)

analysis_layers = context_layers[context_layers.UNIQUE_ID .∈ [unique(analysis_layers_long.UNIQUE_ID)], :]
reefs_with_clusters = (
    (analysis_layers[:, gbr_gcm_cluster_cols[1]] .!= 0) .&
    (analysis_layers[:, gbr_gcm_cluster_cols[2]] .!= 0) .&
    (analysis_layers[:, gbr_gcm_cluster_cols[3]] .!= 0) .&
    (analysis_layers[:, gbr_gcm_cluster_cols[4]] .!= 0) .&
    (analysis_layers[:, gbr_gcm_cluster_cols[5]] .!= 0)
)
analysis_layers = analysis_layers[reefs_with_clusters, :]
consistent_reefs = (
    (analysis_layers[:, gbr_gcm_cluster_cols[1]] .== analysis_layers[:, gbr_gcm_cluster_cols[2]]) .&
    (analysis_layers[:, gbr_gcm_cluster_cols[1]] .== analysis_layers[:, gbr_gcm_cluster_cols[3]]) .&
    (analysis_layers[:, gbr_gcm_cluster_cols[1]] .== analysis_layers[:, gbr_gcm_cluster_cols[4]]) .&
    (analysis_layers[:, gbr_gcm_cluster_cols[1]] .== analysis_layers[:, gbr_gcm_cluster_cols[5]])
)
sum(consistent_reefs) / length(consistent_reefs)

# Identify reefs that have a median depth of 1 to 10m and the proportion of these reefs that are in low-medium cover clusters
reefs_1_to_10 = context_layers[
    (context_layers.depth_qc .== 0) .& 
    (context_layers.depth_med .>= 1) .& 
    (context_layers.depth_med .<= 10), :]
reefs_1_to_10.low_medium_clusters .= 0
reefs_1_to_10_clusters = analysis_layers_long[analysis_layers_long.UNIQUE_ID .∈ [reefs_1_to_10.UNIQUE_ID], :]
reefs_1_to_10_clusters = reefs_1_to_10_clusters[reefs_1_to_10_clusters.variable .!= 0.0, :]

for reef in eachrow(reefs_1_to_10)
    reef_id = reef.UNIQUE_ID
    reef.low_medium_clusters = all(reefs_1_to_10_clusters[reefs_1_to_10_clusters.UNIQUE_ID .== reef_id, :value] .<= 2.0)
end

sum(reefs_1_to_10.low_medium_clusters) / nrow(reefs_1_to_10)