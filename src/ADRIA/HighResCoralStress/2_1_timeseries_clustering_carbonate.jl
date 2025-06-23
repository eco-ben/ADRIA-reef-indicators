using Revise, Infiltrator

using GLMakie, GeoMakie, GraphMakie
using Statistics
using YAXArrays
import GeoDataFrames as GDF

using ADRIA

# include("../../common.jl")
# include("../../plotting_functions.jl")

# GCM = "ACCESS-ESM1-5"
# # Select the target GCM from
# dhw_scenarios = open_dataset("../../ADRIA Domains/GBR_2024_10_15_HighResCoralStress/DHWs/dhwRCP45.nc")

# GBR_domain_path = "../../ADRIA Domains/GBR_2024_10_15_HighResCoralStress/"
# gbr_dom = ADRIA.load_domain(GBR_domain_path, "45")
# context_layers = gbr_dom.loc_data

# # Attach bioregion level data to context_layers
# bioregions = GDF.read("../data/GBRMPA_reefal_bioregions.gpkg")
# context_bioregion = find_intersections(context_layers, bioregions, :GBRMPA_ID, :DESCRIP, :SHAPE)
# context_layers = leftjoin(context_layers, context_bioregion; on=:GBRMPA_ID, matchmissing=:notequal, order=:left)
# rename!(context_layers, :area_ID => :bioregion)
# context_layers.bioregion .= ifelse.(ismissing.(context_layers.bioregion), "NA", context_layers.bioregion)
# context_layers.bioregion = convert.(String, context_layers.bioregion)

# n_clusters = 3

# thresholds = collect(2:1:20) ./ 100

# rs = ADRIA.load_results("../outputs/ADRIA_results/HighResCoralStress/GBR_2024_10_15_HighResCoralStress__RCPs_45_$(GCM)")
# GCM_results = GCM_analysis_results(rs)
# # Convert absolute cover to cover relative to carbonate budjet area
# areas = gbr_dom.loc_data.k .* gbr_dom.loc_data.area
# cover = GCM_results.absolute_median_cover

# n_timesteps = size(cover, 1)
# n_locations = size(cover, 2)
# n_thresholds = length(thresholds)
# threshold_cover = YAXArray(
#     (
#         cover.timesteps,
#         cover.locations,
#         Dim{:thresholds}(thresholds)
#     ),
#     zeros(n_timesteps, n_locations, n_thresholds)
# )
# threshold_clusters = YAXArray(
#     (
#         cover.locations,
#         Dim{:thresholds}(thresholds)
#     ),
#     zeros(n_locations, n_thresholds)
# )

# function carbonate_threshold_cover(area, cover; threshold=0.2)
#     carbonate_cover = area * threshold

#     return (cover .- carbonate_cover) ./ carbonate_cover .* 100
# end

# function threshold_cover_timeseries(areas, cover_timeseries, threshold)
#     threshold_cover = YAXArray(
#         (cover_timeseries.timesteps, cover_timeseries.locations),
#         zeros(size(cover_timeseries))
#     )

#     for loc in eachindex(areas)
#         threshold_cover[:, loc] = carbonate_threshold_cover(areas[loc], cover[:, loc]; threshold=threshold)
#     end

#     return threshold_cover
# end

# for (t, threshold_t) in enumerate(thresholds)
#     for loc in eachindex(areas)
#         threshold_cover[:, loc, t] = carbonate_threshold_cover(areas[loc], cover[:, loc]; threshold=threshold_t)
#     end

#     threshold_clusters[:, t] = first(grouped_timeseries_clustering(threshold_cover[:, :, t], ones(n_locations); n_clusters=3, length_t=1:50))
#     #threshold_clusters[:, t] = ADRIA.analysis.cluster_series(threshold_cover[1:50, :, t], 3)
# end

# length_t = 1:50
# fig = Figure()
# layout = [fldmod1(x, 4) for x in eachindex(thresholds)]

# legend_entries = []
# for (i, col) in enumerate([:green, :orange, :blue])
#     LE = LineElement(; color=col, marker=:circle)
#     push!(legend_entries, [LE])
# end
# Legend(
#     fig[5,4],
#     legend_entries,
#     ["Low", "Medium", "High"],
#     nbanks=3
# )

# for t in eachindex(thresholds)
    
#     ADRIA.viz.clustered_scenarios!(
#             fig[layout[t]...],
#             threshold_cover[length_t, :, t],
#             Int64.(threshold_clusters[:, t].data);
#             opts = Dict{Symbol, Any}(:legend => false),
#             axis_opts = Dict(
#                 :title => string.(thresholds[t]),
#                 :xticks => (
#                     first(length_t):10:last(length_t), 
#                     string.(collect((first(threshold_cover.timesteps):10:2065)))
#                 )
#             )
#     )
# end

# loc_cluster_identical = zeros(n_locations)
# for loc in 1:n_locations
#     loc_clusters = threshold_clusters[loc, :]
#     loc_cluster_mode = mode(loc_clusters)
#     loc_cluster_identical[loc] = sum(loc_clusters .== loc_cluster_mode)
# end
# hist(loc_cluster_identical)

# fig = Figure()
# ax = Axis(fig[1,1], ylabel="log so_to_si", xlabel="number of times each reef has the same cluster ID out of 19")
# scatter!(loc_cluster_identical, context_layers.bioregion_average_latitude)

# #threshold_clusters = threshold_clusters .+ rand(-0.1:0.05:0.1, 3806, 19)
# fig = Figure()
# ax = Axis(fig[1,1], ylabel="Cluster ID", xlabel="cover-threshold used", xticks = (1:19, string.(collect(threshold_clusters.thresholds))))
# series!(threshold_clusters.data; solid_color=(:gray, 0.1))


include("../../common.jl")
includet("../../plotting_functions.jl")

# Select the target GCM from
dhw_scenarios = open_dataset("../../ADRIA Domains/GBR_2024_10_15_HighResCoralStress/DHWs/dhwRCP45.nc")

GBR_domain_path = "../../ADRIA Domains/GBR_2024_10_15_HighResCoralStress/"
gbr_dom = ADRIA.load_domain(GBR_domain_path, "45")
context_layers = gbr_dom.loc_data

# Attach bioregion level data to context_layers
bioregions = GDF.read("../data/GBRMPA_reefal_bioregions.gpkg")
context_bioregion = find_intersections(context_layers, bioregions, :GBRMPA_ID, :DESCRIP, :SHAPE)
context_layers = leftjoin(context_layers, context_bioregion; on=:GBRMPA_ID, matchmissing=:notequal, order=:left)
rename!(context_layers, :area_ID => :bioregion)
context_layers.bioregion .= ifelse.(ismissing.(context_layers.bioregion), "NA", context_layers.bioregion)
context_layers.bioregion = convert.(String, context_layers.bioregion)

n_clusters = 3

areas = gbr_dom.loc_data.area

for GCM in dhw_scenarios.dhw.properties["members"]
    
    # Select GCM and load relevant results
    @info "Performing timeseries clustering for $(GCM)"

    rs = ADRIA.load_results("../outputs/ADRIA_results/HighResCoralStress/GBR_2024_10_15_HighResCoralStress__RCPs_45_$(GCM)")
    GCM_results = GCM_analysis_results(rs)
    absolute_cover = GCM_results.absolute_median_cover

    # threshold_cover = threshold_cover_timeseries(areas, absolute_cover, 0.17)
    threshold_cover = percentage_cover_timeseries(areas, absolute_cover)

    (bioregion_clusters, bioregion_cluster_cats) = grouped_timeseries_clustering(threshold_cover.data, context_layers.bioregion; n_clusters=n_clusters, length_t=1:50)
    (man_region_clusters, man_region_cluster_cats) = grouped_timeseries_clustering(threshold_cover, context_layers.management_area; n_clusters=n_clusters, length_t=1:50)
    (gbr_clusters, gbr_cluster_cats) = grouped_timeseries_clustering(threshold_cover, ones(size(context_layers, 1)); n_clusters=n_clusters, length_t=1:50)

    context_layers[:, "$(GCM)_bioregion_clusters"] = bioregion_clusters
    context_layers[:, "$(GCM)_bioregion_cluster_cats"] = bioregion_cluster_cats
    context_layers[:, "$(GCM)_management_area_clusters"] = man_region_clusters
    context_layers[:, "$(GCM)_management_area_cluster_cats"] = man_region_cluster_cats
    context_layers[:, "$(GCM)_gbr_clusters"] = gbr_clusters
    context_layers[:, "$(GCM)_gbr_cluster_cats"] = gbr_cluster_cats

end

GDF.write("../outputs/ADRIA_results/HighResCoralStress/bellwether_reefs_carbonate.gpkg", context_layers)

# threshold_cover = absolute_cover
# for loc in eachindex(absolute_cover.locations)
#     threshold_cover[:, loc] = normalise(threshold_cover[:, loc], (0,1))
# end