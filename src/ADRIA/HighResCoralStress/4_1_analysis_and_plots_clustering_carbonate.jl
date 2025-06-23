using Revise, Infiltrator

using GLMakie, GeoMakie, GraphMakie
using Statistics
using YAXArrays
using CairoMakie

using ADRIA

include("../../common.jl")
include("../../plotting_functions.jl")

CairoMakie.activate!()

# Select the target GCM from
dhw_scenarios = open_dataset("../../ADRIA Domains/GBR_2024_10_15_HighResCoralStress/DHWs/dhwRCP45.nc")
GCMs = dhw_scenarios.dhw.properties["members"]

context_layers = GDF.read("../outputs/ADRIA_results/HighResCoralStress/analysis_context_layers_carbonate.gpkg")
context_layers.gbr .= "Great Barrier Reef"
context_layers.log_so_to_si = log10.(context_layers.so_to_si)
context_layers.log_total_strength = log10.(context_layers.total_strength)

gbr_dom = ADRIA.load_domain("../../ADRIA Domains/GBR_2024_10_15_HighResCoralStress/", "45")
areas = gbr_dom.loc_data.area

for (i_gcm, GCM) in enumerate(GCMs)
    # Select GCM and load relevant results
    @info "analysing reef clustering for $(GCM)"

    fig_out_dir = "../outputs/ADRIA_results/HighResCoralStress/figs/$(GCM)"

    rs = ADRIA.load_results("../outputs/ADRIA_results/HighResCoralStress/GBR_2024_10_15_HighResCoralStress__RCPs_45_$(GCM)")
    GCM_results = GCM_analysis_results(rs)
    absolute_cover = GCM_results.absolute_median_cover
    # threshold_cover = threshold_cover_timeseries(areas, absolute_cover, 0.17)
    threshold_cover = percentage_cover_timeseries(areas, absolute_cover)

    dhw_ts = gbr_dom.dhw_scens[:,:,i_gcm]
    ax = (
        Dim{:timesteps}(2025:2099),
        Dim{:locations}(collect(GCM_results.relative_cover.locations))
    )
    dhw_ts = rebuild(dhw_ts, dims=ax)

    # Subset layers to remove reefs that start with very low coral cover
    analysis_layers = context_layers[(context_layers.management_area .!= "NA") .& (context_layers.bioregion .!= "NA"), :]

    # Analysing clusters identified at bioregion level
    cluster_analysis_plots(GCM, analysis_layers, threshold_cover, dhw_ts, :bioregion, fig_out_dir)
    cluster_analysis_plots(GCM, analysis_layers, threshold_cover, dhw_ts, :management_area, fig_out_dir)
    cluster_analysis_plots(GCM, analysis_layers, threshold_cover, dhw_ts, :gbr, fig_out_dir)

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
        GCM_analysis_results(
            ADRIA.load_results("../outputs/ADRIA_results/HighResCoralStress/GBR_2024_10_15_HighResCoralStress__RCPs_45_$(GCM)")
            ).absolute_median_cover
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
    "../outputs/ADRIA_results/HighResCoralStress/figs/GCM_timeseries_plot.png", 
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

# fig = Figure()
# ax = Axis(
#     fig[1,1],
#     xlabel = "GCM",
#     ylabel = "Cluster",
#     xticks = (1:length(GCMs), GCMs)
# )
# man_cluster_changes = [unique(collect(row[gcm_cluster_cols])) for row in eachrow(analysis_layers)]
# lines!.(ax, cluster_changes; color = (:gray, 0.1))

# function transition_count_matrix(v, unique_vals)
#     mat = zeros(length(unique_vals), length(unique_vals))

#     for i in 1:length(v)-1
#         current_val = Int64(v[i])
#         next_val = Int64(v[i + 1])
#         if current_val != next_val
#             mat[current_val, next_val] += 1
#         end
#     end

#     return mat
# end
# cluster_transitions = [transition_count_matrix(vals, 1:3) for vals in cluster_changes]
# cluster_transitions_mat = zeros(3,3)
# for matrix in cluster_transitions
#     cluster_transitions_mat = cluster_transitions_mat + matrix
# end
# cluster_trans_props = cluster_transitions_mat ./ sum(cluster_transitions_mat)

# unique_cluster_changes = [(maximum(vals) - minimum(vals)) for vals in cluster_changes]
# sum(unique_cluster_changes .== 0) / length(unique_cluster_changes)
# sum(unique_cluster_changes .== 1) / length(unique_cluster_changes)
# sum(unique_cluster_changes .== 2) / length(unique_cluster_changes)


year_cols = Vector{String}()
for GCM in GCMs
    for threshold in thresholds
        push!(year_cols, "$(GCM)_years_above_$(threshold)")
    end
end

long_depth_reefs = stack(
    context_layers[:, ["UNIQUE_ID", "depth_med", year_cols...]],
    year_cols
)
ACCESS_ESM1_5_long = long_depth_reefs[contains.(long_depth_reefs.variable, "ACCESS-ESM1-5"), :]
ACCESS_ESM1_5_long.variable = last.(split.(ACCESS_ESM1_5_long.variable, "_"))

fig = Figure()
ax = Axis(
    fig[1,1],
    ylabel = "Years above carbonate budget threshold",
    xlabel = "Carbonate budget threshold (%)",
    xticks = (1:1:11, unique(ACCESS_ESM1_5_long.variable))
)
rain = rainclouds!(
    ACCESS_ESM1_5_long.variable, 
    ACCESS_ESM1_5_long.value;
    color=ACCESS_ESM1_5_long.depth_med, 
    plot_boxplots=false, 
    clouds=nothing,
    show_median=false, 
    markersize=6, 
    jitter_width=0.9
)
# scat = scatter!(
#     ACCESS_ESM1_5_long.variable,
#     ACCESS_ESM1_5_long.value,
#     color=ACCESS_ESM1_5_long.depth_med,
#     alpha = 0.5
# )
Colorbar(fig[1,2], rain, label = "Median reef depth (m)")