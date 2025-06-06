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

context_layers = GDF.read("../outputs/ADRIA_results/HighResCoralStress/analysis_context_layers.gpkg")
context_layers.gbr .= "Great Barrier Reef"
context_layers.log_so_to_si = log10.(context_layers.so_to_si)
context_layers.log_total_strength = log10.(context_layers.total_strength)

gbr_dom = ADRIA.load_domain("../../ADRIA Domains/GBR_2024_10_15_HighResCoralStress/", "45")

for (i_gcm, GCM) in enumerate(GCMs)
    # Select GCM and load relevant results
    @info "analysing reef clustering for $(GCM)"

    fig_out_dir = "../outputs/ADRIA_results/HighResCoralStress/figs/$(GCM)"

    rs = ADRIA.load_results("../outputs/ADRIA_results/HighResCoralStress/GBR_2024_10_15_HighResCoralStress__RCPs_45_$(GCM)")
    GCM_results = GCM_analysis_results(rs)
    rel_cover = GCM_results.relative_cover

    rel_cover_less_than_5 = [all(rel_cover[:, i].data .< 5) for i in 1:size(rel_cover, 2)]
    rel_cover_less_than_5_ind = findall(rel_cover_less_than_5)

    dhw_ts = gbr_dom.dhw_scens[:,:,i_gcm]
    ax = (
        Dim{:timesteps}(2025:2099),
        Dim{:locations}(collect(GCM_results.relative_cover.locations))
    )
    dhw_ts = rebuild(dhw_ts, dims=ax)

    # Subset layers to remove reefs that start with very low coral cover
    analysis_layers = context_layers[rel_cover_less_than_5, :]
    analysis_layers = analysis_layers[(analysis_layers.management_area .!= "NA") .& (analysis_layers.bioregion .!= "NA"), :]

    # Analysing clusters identified at bioregion level
    cluster_analysis_plots(GCM, analysis_layers, rel_cover, dhw_ts, :bioregion, fig_out_dir)
    cluster_analysis_plots(GCM, analysis_layers, rel_cover, dhw_ts, :management_area, fig_out_dir)
    cluster_analysis_plots(GCM, analysis_layers, rel_cover, dhw_ts, :gbr, fig_out_dir)

end

man_area_gcm_cluster_cols = [Symbol("$(GCM)_management_area_clusters") for GCM in GCMs]
bioregion_gcm_cluster_cols = [Symbol("$(GCM)_bioregion_clusters") for GCM in GCMs]
gbr_gcm_cluster_cols = [Symbol("$(GCM)_gbr_clusters") for GCM in GCMs]
analysis_layers_long = stack(
    context_layers[:, [:UNIQUE_ID, :management_area, man_area_gcm_cluster_cols...]], 
    man_area_gcm_cluster_cols
)
analysis_layers_long.GCM = [first(split(name, "_")) for name in analysis_layers_long.variable]

rel_cover_arrays = [GCM_analysis_results(ADRIA.load_results("../outputs/ADRIA_results/HighResCoralStress/GBR_2024_10_15_HighResCoralStress__RCPs_45_$(GCM)")).relative_cover for GCM in GCMs]
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
    (analysis_layers[:, gcm_cluster_cols[1]] .!= 0) .&
    (analysis_layers[:, gcm_cluster_cols[2]] .!= 0) .&
    (analysis_layers[:, gcm_cluster_cols[3]] .!= 0) .&
    (analysis_layers[:, gcm_cluster_cols[4]] .!= 0) .&
    (analysis_layers[:, gcm_cluster_cols[5]] .!= 0)
)
analysis_layers = analysis_layers[reefs_with_clusters, :]
consistent_reefs = (
    (analysis_layers[:, gcm_cluster_cols[1]] .== analysis_layers[:, gcm_cluster_cols[2]]) .&
    (analysis_layers[:, gcm_cluster_cols[1]] .== analysis_layers[:, gcm_cluster_cols[3]]) .&
    (analysis_layers[:, gcm_cluster_cols[1]] .== analysis_layers[:, gcm_cluster_cols[4]]) .&
    (analysis_layers[:, gcm_cluster_cols[1]] .== analysis_layers[:, gcm_cluster_cols[5]])
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
analysis_layers_depth = analysis_layers[analysis_layers.depth_mean .!= 7, :]
groups_too_few_clusters_depth = grouping_counts(
    grouping, 
    analysis_layers_depth, 
    "$(GCM)_$(grouping)_clusters", 
    3,
    5
)
analysis_layers_depth = analysis_layers_depth[
    analysis_layers_depth[:, grouping] .∉ [groups_too_few_clusters_depth], :
]
threshold = 10

bioregion_colors = distinguishable_colors(length(unique(analysis_layers_depth.bioregion)))
bioregion_colors = 

fig = Figure()
ax1 = Axis(
    fig[1,1];
    xlabel="median reef depth",
    ylabel="number of years above threshold $(threshold)"
)
scat1 = scatter!(
    ax1, 
    analysis_layers_depth.depth_med, 
    analysis_layers_depth[:, "$(GCM)_years_above_$(threshold)"];
    color=analysis_layers_depth[:, "$(GCM)_gbr_clusters"],
    alpha=0.5
)
Colorbar(fig[1,1][1,2], scat1, label="GBR-scale cluster")

ax2 = Axis(
    fig[1,2];
    xlabel="Log total connectivity strength",
    ylabel="number of years above threshold $(threshold)"
)
scat2 = scatter!(
    ax2, 
    analysis_layers_depth.log_total_strength,
    analysis_layers_depth[:, "$(GCM)_years_above_$(threshold)"];
    color=log.(analysis_layers_depth.area),
    alpha=0.7
)
Colorbar(fig[1,2][1,2], scat2, label="Log reef area")