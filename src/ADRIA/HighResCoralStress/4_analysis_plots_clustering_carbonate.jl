"""
Create analysis plots for clustering and carbonate budget analyses. Performed for each GCM separately,
and multi-GCM comparison plot created at end of the script.
"""

include("../../common.jl")

CairoMakie.activate!()

# # Select the target GCM from
dhw_scenarios = open_dataset(joinpath(gbr_domain_path, "DHWs/dhwRCP45.nc"))
GCMs = dhw_scenarios.dhw.properties["members"]

context_layers = GDF.read(joinpath(output_path, "analysis_context_layers_carbonate.gpkg"))
context_layers.gbr .= "Great Barrier Reef"
context_layers.gbr_average_latitude .= 0.0
context_layers.log_so_to_si = log10.(context_layers.so_to_si)

gbr_dom = ADRIA.load_domain(gbr_domain_path, "45")
gbr_dom_filtered = gbr_dom.loc_data[gbr_dom.loc_data.UNIQUE_ID.∈[context_layers.UNIQUE_ID], :]
filtered_indices = indexin(gbr_dom_filtered.UNIQUE_ID, gbr_dom.loc_data.UNIQUE_ID)

thresholds = 10:1:20
year_cols = Vector{String}()
conn_cols = Vector{String}()
for GCM in GCMs
    context_layers[:, "$(GCM)_weighted_incoming_conn_log"] = log10.(context_layers[:, "$(GCM)_weighted_incoming_conn"]) 
    context_layers[:, "$(GCM)_weighted_outgoing_conn_log"] = log10.(context_layers[:, "$(GCM)_weighted_outgoing_conn"]) 
    
    push!(conn_cols, "$(GCM)_weighted_incoming_conn_log")
    for threshold in thresholds
        push!(year_cols, "$(GCM)_years_above_$(threshold)")
    end
end

reefs_long = stack(
    context_layers[:, ["UNIQUE_ID", "depth_med", conn_cols..., year_cols...]],
    year_cols
)

for (i_gcm, GCM) in enumerate(GCMs)
    # Select GCM and load relevant results
    @info "Analysing reef clustering for $(GCM)"

    fig_out_dir = joinpath(figs_path, "$(GCM)")

    absolute_cover = readcubedata(open_dataset(joinpath(output_path, "processed_model_outputs/median_cover_$(GCM).nc")).layer)
    absolute_cover = absolute_cover[locations=At(context_layers.UNIQUE_ID)]
    percentage_cover = percentage_cover_timeseries(context_layers.area, absolute_cover)

    dhw_ts = gbr_dom.dhw_scens[:, filtered_indices, i_gcm]
    dhw_ts = rebuild(dhw_ts, dims=percentage_cover.axes, metadata=dhw_timeseries_properties)

    # Analysing clusters identified at bioregion level
    cluster_analysis_plots(GCM, context_layers, percentage_cover, dhw_ts, :bioregion, fig_out_dir)
    cluster_analysis_plots(GCM, context_layers, percentage_cover, dhw_ts, :management_area, fig_out_dir)
    cluster_analysis_plots(GCM, context_layers, percentage_cover, dhw_ts, :gbr, fig_out_dir)

    gcm_reefs_long = reefs_long[contains.(reefs_long.variable, GCM), :]
    gcm_reefs_long.variable = last.(split.(gcm_reefs_long.variable, "_"))

    depth_year_correlation = Dict([
        (thresh, corspearman(
            gcm_reefs_long[gcm_reefs_long.variable.==thresh, :value],
            gcm_reefs_long[gcm_reefs_long.variable.==thresh, :depth_med]
        ))
        for thresh in string.(thresholds)
    ])

    # depth_carbonate_scatter = carbonate_budget_variable_scatter(
    #     gcm_reefs_long,
    #     :depth_med,
    #     :value,
    #     :variable,
    #     depth_year_correlation;
    #     color_label="Median reef depth [m]"
    # )
    # save(
    #     joinpath(fig_out_dir, "depth_carbonate_scatter.png"),
    #     depth_carbonate_scatter,
    #     px_per_unit=dpi
    # )


    log_incoming_conn_year_correlation = Dict([
        (thresh, corspearman(
            gcm_reefs_long[gcm_reefs_long.variable.==thresh, :value],
            gcm_reefs_long[gcm_reefs_long.variable.==thresh, "$(GCM)_weighted_incoming_conn_log"]
        ))
        for thresh in string.(thresholds)
    ])
    # total_connectivity_carbonate_scatter = carbonate_budget_variable_scatter(
    #     gcm_reefs_long,
    #     :log_total_strength,
    #     :value,
    #     :variable,
    #     log_strength_year_correlation;
    #     color_label="Log total connectivity strength"
    # )
    depth_conn_carbonate_comparison = carbonate_budget_variable_scatter(
        gcm_reefs_long,
        :value,
        :variable,
        depth_year_correlation,
        log_incoming_conn_year_correlation;
        conn_var_col="$(GCM)_weighted_incoming_conn_log"
    )
    save(
        joinpath(fig_out_dir, "depth_conn_carbonate_comparison.png"),
        depth_conn_carbonate_comparison,
        px_per_unit=dpi
    )
    # save(
    #     joinpath(fig_out_dir, "log_total_strength_carbonate_scatter.png"),
    #     total_connectivity_carbonate_scatter,
    #     px_per_unit=dpi
    # )
    gcm_threshold_columns = ["$(GCM)_years_above_$(x)" for x in 10:20]
    context_layers[!, "$(GCM)_mean_positive_years"] = vec(mean(Matrix(context_layers[:, gcm_threshold_columns]), dims=2))
    fig, ax, scat = scatter(context_layers.depth_med, context_layers[:, "$(GCM)_weighted_incoming_conn_log"], color=context_layers[:, "$(GCM)_mean_positive_years"], alpha=0.5)
    Colorbar(fig[1,2], scat, label="Mean number of positive carbonate budget years \nacross θ")
    ax.ylabel = "Log10 weighted incoming connectivity"
    ax.xlabel = "Median depth [m]"
    save(
        joinpath(fig_out_dir, "depth_conn_mean_carbonate_budget.png"),
        fig,
        px_per_unit=dpi
    )

    incoming_conn_map = map_gbr_reefs_cont(
        context_layers, 
        "$(GCM)_weighted_incoming_conn_log", 
        "Log10 weighted incoming connectivity"
    )
    save(
        joinpath(fig_out_dir, "incoming_conn_input_map.png"),
        incoming_conn_map,
        px_per_unit=dpi
    )

    outgoing_conn_map = map_gbr_reefs_cont(
        context_layers, 
        "$(GCM)_weighted_outgoing_conn_log", 
        "Log10 weighted outgoing connectivity"
    )
    save(
        joinpath(fig_out_dir, "outgoing_conn_input_map.png"),
        outgoing_conn_map,
        px_per_unit=dpi
    )

    mean_dhw_map = map_gbr_reefs_cont(
        context_layers,
        "$(GCM)_mean_dhw",
        "Mean DHW [\u00B0C - Weeks]"
    )
    save(
        joinpath(fig_out_dir, "mean_dhw_input_map.png"),
        mean_dhw_map,
        px_per_unit=dpi
    )
end

context_layers.abs_k_area = context_layers.area .* context_layers.k ./ 1e6
vars = [:depth_med, :log_so_to_si, :abs_k_area]
var_labels = Dict(
    :depth_med => "Median depth [m]",
    :log_so_to_si => "Log10 weighted incoming connectivity",
    :abs_k_area => "Carrying capacity [km²]"
)
for var in vars 
    var_map = map_gbr_reefs_cont(context_layers, var, var_labels[var])
    save(
        joinpath(figs_path, "$(var)_input_map.png"),
        var_map,
        px_per_unit=dpi
    )
end

man_area_gcm_cluster_cols = [Symbol("$(GCM)_management_area_clusters") for GCM in GCMs]
bioregion_gcm_cluster_cols = [Symbol("$(GCM)_bioregion_clusters") for GCM in GCMs]
gbr_gcm_cluster_cols = [Symbol("$(GCM)_gbr_clusters") for GCM in GCMs]

# Prepare management area grouped clusters and timeseries for plotting
analysis_layers_long = stack(
    context_layers[:, [:UNIQUE_ID, :management_area, :management_area_average_latitude, man_area_gcm_cluster_cols...]],
    man_area_gcm_cluster_cols
)
analysis_layers_long.GCM = [first(split(name, "_")) for name in analysis_layers_long.variable]

# Collate all cover timeseries and plot them grouped by management area
rel_cover_arrays = [
    percentage_cover_timeseries(
        gbr_dom.loc_data.area,
        readcubedata(open_dataset(joinpath(output_path, "processed_model_outputs/median_cover_$(GCM).nc")).layer)
    ) for GCM in GCMs
]
rel_cover_arrays = concatenatecubes(rel_cover_arrays, Dim{:GCM}(GCMs))

dhw_arrays = [
    rebuild(gbr_dom.dhw_scens[:, :, i_gcm], dims=rel_cover_arrays.axes[1:2]) for i_gcm in eachindex(GCMs)
]
dhw_arrays = concatenatecubes(dhw_arrays, Dim{:GCM}(GCMs))
dhw_arrays = rebuild(dhw_arrays, metadata=dhw_timeseries_properties)

rel_cover_arrays = rel_cover_arrays[locations=At(analysis_layers_long.UNIQUE_ID)]
dhw_arrays = dhw_arrays[locations=At(analysis_layers_long.UNIQUE_ID)]

GCM_comparison_cover_plot = grouped_GCM_cluster_timeseries_plots(
    rel_cover_arrays,
    analysis_layers_long,
    :value,
    [:management_area, :GCM],
    1:50
)
save(
    joinpath(figs_path, "GCM_cover_timeseries_plot.png"),
    GCM_comparison_cover_plot,
    px_per_unit=dpi
)

GCM_comparison_dhw_plot = grouped_GCM_cluster_timeseries_plots(
    dhw_arrays,
    analysis_layers_long,
    :value,
    [:management_area, :GCM],
    1:50
)
save(
    joinpath(figs_path, "GCM_dhw_timeseries_plot.png"),
    GCM_comparison_dhw_plot,
    px_per_unit=dpi
)

# analysis_layers = context_layers[context_layers.UNIQUE_ID.∈[unique(analysis_layers_long.UNIQUE_ID)], :]
# reefs_with_clusters = (
#     (analysis_layers[:, gbr_gcm_cluster_cols[1]] .!= 0) .&
#     (analysis_layers[:, gbr_gcm_cluster_cols[2]] .!= 0) .&
#     (analysis_layers[:, gbr_gcm_cluster_cols[3]] .!= 0) .&
#     (analysis_layers[:, gbr_gcm_cluster_cols[4]] .!= 0) .&
#     (analysis_layers[:, gbr_gcm_cluster_cols[5]] .!= 0)
# )
# analysis_layers = analysis_layers[reefs_with_clusters, :]
# consistent_reefs = (
#     (analysis_layers[:, gbr_gcm_cluster_cols[1]] .== analysis_layers[:, gbr_gcm_cluster_cols[2]]) .&
#     (analysis_layers[:, gbr_gcm_cluster_cols[1]] .== analysis_layers[:, gbr_gcm_cluster_cols[3]]) .&
#     (analysis_layers[:, gbr_gcm_cluster_cols[1]] .== analysis_layers[:, gbr_gcm_cluster_cols[4]]) .&
#     (analysis_layers[:, gbr_gcm_cluster_cols[1]] .== analysis_layers[:, gbr_gcm_cluster_cols[5]])
# )
# sum(consistent_reefs) / length(consistent_reefs)

# # Identify how many 'high' cluster reefs move to lower clusters under different GCMs (based on bioregion analyses).
# analysis_layers = context_layers
# reefs_with_clusters = (
#     (analysis_layers[:, bioregion_gcm_cluster_cols[1]] .!= 0) .&
#     (analysis_layers[:, bioregion_gcm_cluster_cols[2]] .!= 0) .&
#     (analysis_layers[:, bioregion_gcm_cluster_cols[3]] .!= 0) .&
#     (analysis_layers[:, bioregion_gcm_cluster_cols[4]] .!= 0) .&
#     (analysis_layers[:, bioregion_gcm_cluster_cols[5]] .!= 0)
# )
# analysis_layers = analysis_layers[reefs_with_clusters, :]
# high_cluster_reefs_any_gcm = analysis_layers[(
#     (analysis_layers[:, bioregion_gcm_cluster_cols[1]] .== 3) .|
#     (analysis_layers[:, bioregion_gcm_cluster_cols[2]] .== 3) .|
#     (analysis_layers[:, bioregion_gcm_cluster_cols[3]] .== 3) .|
#     (analysis_layers[:, bioregion_gcm_cluster_cols[4]] .== 3) .|
#     (analysis_layers[:, bioregion_gcm_cluster_cols[5]] .== 3)
# ), :]
# high_to_low_reefs = high_cluster_reefs_any_gcm[(
#     (high_cluster_reefs_any_gcm[:, bioregion_gcm_cluster_cols[1]] .== 1) .|
#     (high_cluster_reefs_any_gcm[:, bioregion_gcm_cluster_cols[2]] .== 1) .|
#     (high_cluster_reefs_any_gcm[:, bioregion_gcm_cluster_cols[3]] .== 1) .|
#     (high_cluster_reefs_any_gcm[:, bioregion_gcm_cluster_cols[4]] .== 1) .|
#     (high_cluster_reefs_any_gcm[:, bioregion_gcm_cluster_cols[5]] .== 1)
# ), :]

# nrow(switching_high_reefs) / nrow(high_cluster_reefs_any_gcm)

# Identify reefs that have a median depth of 1 to 10m and the proportion of these reefs that are in low-medium cover clusters (based on bioregion analyses).
analysis_layers_long_bioregion = stack(
    context_layers[:, [:UNIQUE_ID, :management_area, bioregion_gcm_cluster_cols...]],
    bioregion_gcm_cluster_cols
)

reefs_1_to_10 = context_layers[
    (context_layers.depth_med.>=1).&(context_layers.depth_med.<=10), :]
reefs_1_to_10.low_medium_clusters .= 0
reefs_1_to_10_clusters = analysis_layers_long_bioregion[analysis_layers_long_bioregion.UNIQUE_ID.∈[reefs_1_to_10.UNIQUE_ID], :]
reefs_1_to_10_clusters = reefs_1_to_10_clusters[reefs_1_to_10_clusters.variable.!=0.0, :]

for reef in eachrow(reefs_1_to_10)
    reef_id = reef.UNIQUE_ID
    reef.low_medium_clusters = all(reefs_1_to_10_clusters[reefs_1_to_10_clusters.UNIQUE_ID.==reef_id, :value] .<= 2.0)
end

sum(reefs_1_to_10.low_medium_clusters) / nrow(reefs_1_to_10)

# Quanitifying and visualising uncertainty in connectivity data
conn_files = readdir(joinpath(gbr_domain_path, "connectivity"))
conn_names = [split(fn, "_")[2] * "_" * split(split(fn, "_")[3], "t")[1] for fn in conn_files]
conn_files = [joinpath(gbr_domain_path, "connectivity", file) for file in conn_files]
mean_conn = gbr_dom.conn

function read_conn_file(fn)
    return Matrix(
        CSV.read(
            fn,
            DataFrame;
            comment="#",
            missingstring="NA",
            transpose=false,
            types=Float64,
            drop=[1]
        )
    )
end

all_conn_data = [
    YAXArray(
        mean_conn.axes,
        read_conn_file(c_file)
    ) for c_file in conn_files
]
all_conn_data = concatenatecubes(all_conn_data, Dim{:conn_name}(conn_names))
all_conn_data = all_conn_data[Source=At(context_layers.UNIQUE_ID), Sink=At(context_layers.UNIQUE_ID)]
mean_conn = mean_conn[Source=At(context_layers.UNIQUE_ID), Sink=At(context_layers.UNIQUE_ID)]

std_conn = mapslices(std, all_conn_data, dims=[:conn_name])
rstd_conn = std_conn ./ mean_conn .* 100

rstd_matrix = rstd_conn
fig = Figure(size=(14.82 * centimetre, 14.82 * centimetre), fontsize=fontsize)
ax = Axis(
    fig[1, 1],
    xlabel="",
    ylabel=""
)
hm = heatmap!(ax, rstd_matrix)
Colorbar(fig[1, 2], hm, label="Relative standard deviation connection strength [%]")
save(
    joinpath(figs_path, "connectivity_rsd_percent.png"),
    fig,
    px_per_unit=dpi
)

# Plot cluster assignment heatmap for bioregions
bioregion_gcm_clusters = gcm_cluster_assignment_heatmap(
    context_layers,
    :bioregion,
    bioregion_gcm_cluster_cols
)
save(
    joinpath(figs_path, "GCM_bioregion_cluster_assignments.png"),
    bioregion_gcm_clusters,
    px_per_unit=dpi
)

# Reefs removed - Indicate cutoff based on depth
per_pct_count = [count((i - 0.01) .<= gbr_dom.loc_data.depth_rast_prop .<= i) for i in 0.01:0.01:1.0] ./ length(gbr_dom.loc_data.depth_rast_prop) * 100.0
f, ax, sp = barplot(per_pct_count)
vlines!(ax, 5; color=(:red, 0.5))
text!((10.0, 2.0); text="Cutoff", align=(:center, :bottom))
ax.ylabel = "Proportion of Reefs [%]"
ax.xlabel = "Bathymetry data coverage [%]"

save(joinpath(figs_path, "reef_depth_included_cutoff.png"), f, px_per_unit=dpi)

# Plot a map of reef cluster assignments across GCMs at management area scale:
consistent_reefs_colors = [:gray, :green, :orange, :blue];
order = Dict("variable" => 1, "low" => 2, "medium" => 3, "high" => 4)
sorted_vals = sort(context_layers, :man_area_consistent_reefs, by=x -> order[x])
gbr_assignment_map = map_gbr_reefs_cat(sorted_vals, :man_area_consistent_reefs, consistent_reefs_colors, "Cluster assignment across GCMs")

save(joinpath(figs_path, "gcm_cluster_assignment_map.png"), gbr_assignment_map, px_per_unit=dpi)

n_consistent_reefs = combine(groupby(context_layers[context_layers.man_area_consistent_reefs.!="variable", :], :management_area), nrow => :nrow)
n_consistent_reefs.percentage = n_consistent_reefs.nrow ./ sum(n_consistent_reefs.nrow) .* 100

n_consistent_reefs = combine(groupby(context_layers[context_layers.man_area_consistent_reefs.!="variable", :], :man_area_consistent_reefs), nrow => :nrow)
n_consistent_reefs.percentage = n_consistent_reefs.nrow ./ sum(n_consistent_reefs.nrow) .* 100
