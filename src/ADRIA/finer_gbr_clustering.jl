"""
Try finer scale temporal timeseries clustering across the whole gbr spatial scale.
"""

using YAXArrays, DimensionalData

using DataFrames, Statistics

using ADRIA, CoralBlox

include("common.jl")

context_layers = find_latest_file("../canonical-reefs/output/")
context_layers = GDF.read(context_layers)

rs = ADRIA.load_results("outputs/ADRIA-out/ReefMod Engine__RCPs_45__2024-06-19_10_58_55_647")

fig_opts = Dict(:size => (1600, 800))



# s_tac = ADRIA.metrics.scenario_total_cover(rs)
tac = ADRIA.metrics.total_absolute_cover(rs)

# reduce all the scenarios down to one series for each reef
tac_sites = mapslices_toFloat64(median, tac, :scenarios)

# Have to remove the first year as there seems to be an issue with that year's data
tac_sites_reduced = tac_sites[timesteps=2:79]

# calculate the relative site cover from the initial cover across timesteps for each reef
rel_cover = mapslices_toFloat64(relative_site_cover, tac_sites_reduced, :timesteps)

n_clusters = 15
gbr_clusters = ADRIA.analysis.cluster_series(rel_cover, n_clusters)
gbr_lagged_clusters = lagged_cluster_analysis(rel_cover, gbr_clusters, 1:10)

reefs_p9_l7 = gbr_lagged_clusters[(gbr_lagged_clusters.lag7 .>= 0.9), :UNIQUE_ID]
f, ga = plot_map(context_layers, :management_area)
plot_map!(ga, context_layers[(context_layers.UNIQUE_ID .∈ [reefs_p9_l7]),:], color=:black)
plot_map(context_layers[(context_layers.UNIQUE_ID .∈ [reefs_p9_l7]), :], :closest_port)
