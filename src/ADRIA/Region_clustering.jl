"""
Initial GBR wide and management-area scale timeseries clustering approach.
"""
using CSV

using GLMakie, GeoMakie, GraphMakie

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
n_clusters = 5

# Whole - GBR clustering
clusters = ADRIA.analysis.cluster_series(rel_cover, n_clusters)

tsc_fig = ADRIA.viz.clustered_scenarios(rel_cover, clusters)


# Far Northern Clustering
FN_reefs = context_layers[(context_layers.management_area .== "Far Northern Management Area"), :UNIQUE_ID]
FN_rel_cover = rel_cover[:,rel_cover.sites .∈ [FN_reefs]]

FN_clusters = ADRIA.analysis.cluster_series(FN_rel_cover, n_clusters)
f = ADRIA.viz.clustered_scenarios(
    FN_rel_cover, FN_clusters; fig_opts=fig_opts, axis_opts=Dict(:ylabel => "Far Northern Clustered Proportion of Initial Cover")
)
save("figs/initial_FN_clusters.png", f)

series!(FN_rel_cover.data', solid_color=:black)


# Cairns/Cooktown Clustering
CC_reefs = context_layers[(context_layers.management_area .== "Cairns/Cooktown Management Area"), :UNIQUE_ID]
CC_rel_cover = rel_cover[:,rel_cover.sites .∈ [CC_reefs]]

CC_clusters = ADRIA.analysis.cluster_series(CC_rel_cover, n_clusters)
f = ADRIA.viz.clustered_scenarios(
    CC_rel_cover, CC_clusters; fig_opts=fig_opts, axis_opts=Dict(:ylabel => "Cairns/Cooktown Clustered Proportion of Initial Cover")
)
save("figs/initial_CC_clusters.png", f)

series!(CC_rel_cover.data', solid_color=:black)


# Townsville/Whitsunday Clustering
TSV_reefs = context_layers[(context_layers.management_area .== "Townsville/Whitsunday Management Area"), :UNIQUE_ID]
TSV_rel_cover = rel_cover[:,rel_cover.sites .∈ [TSV_reefs]]

TSV_clusters = ADRIA.analysis.cluster_series(TSV_rel_cover, n_clusters)
f= ADRIA.viz.clustered_scenarios(
    TSV_rel_cover, TSV_clusters; fig_opts=fig_opts, axis_opts=Dict(:ylabel => "Townsville/Whitsunday Clustered Proportion of Initial Cover")
)
save("figs/initial_TSV_clusters.png", f)

series!(TSV_rel_cover.data', solid_color=:black)


# Mackay/Capricorn Clustering
MCap_reefs = context_layers[(context_layers.management_area .== "Mackay/Capricorn Management Area"), :UNIQUE_ID]
MCap_rel_cover = rel_cover[:,rel_cover.sites .∈ [MCap_reefs]]

MCap_clusters = ADRIA.analysis.cluster_series(MCap_rel_cover, n_clusters)
f= ADRIA.viz.clustered_scenarios(
    MCap_rel_cover, MCap_clusters; fig_opts=fig_opts, axis_opts=Dict(:ylabel => "Mackay/Capricorn Clustered Proportion of Initial Cover")
)
save("figs/initial_MCap_clusters.png", f)

series!(MCap_rel_cover.data', solid_color=:black)

#ADRIA.viz.explore(rs)

FN_lagged_clusters = lagged_cluster_analysis(FN_rel_cover, FN_clusters, [1, 2, 5, 10])
CC_lagged_clusters = lagged_cluster_analysis(CC_rel_cover, CC_clusters, [1, 2, 5, 10])
TSV_lagged_clusters = lagged_cluster_analysis(TSV_rel_cover, TSV_clusters, [1, 2, 5, 10])
MCap_lagged_clusters = lagged_cluster_analysis(MCap_rel_cover, MCap_clusters, 5:15)

FN_lagged_region = lagged_region_analysis(FN_rel_cover, "Far Northern", [1, 2, 5, 10])
CC_lagged_region = lagged_region_analysis(CC_rel_cover, "Cairns-Cooktown", [1, 2, 5, 10])
TSV_lagged_region = lagged_region_analysis(TSV_rel_cover, "Townsville-Whitsunday", [1, 2, 5, 10])
MCap_lagged_region = lagged_region_analysis(MCap_rel_cover, "Mackay-Capricorn", 5:15)

plot_lines()
