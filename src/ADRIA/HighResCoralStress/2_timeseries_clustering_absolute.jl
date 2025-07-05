"""
Cluster reef cover timeseries at three different scales across the GBR. For clustering,
reef cover as a percentage of total reef area is used.
"""

using Revise, Infiltrator

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

    absolute_cover = open_dataset("../outputs/ADRIA_results/HighResCoralStress/processed_model_outputs/median_cover_$(GCM).nc").layer

    # threshold_cover = threshold_cover_timeseries(areas, absolute_cover, 0.17)
    percentage_cover = percentage_cover_timeseries(areas, absolute_cover)

    (bioregion_clusters, bioregion_cluster_cats) = grouped_timeseries_clustering(percentage_cover.data, context_layers.bioregion; n_clusters=n_clusters, length_t=1:50)
    (man_region_clusters, man_region_cluster_cats) = grouped_timeseries_clustering(percentage_cover, context_layers.management_area; n_clusters=n_clusters, length_t=1:50)
    (gbr_clusters, gbr_cluster_cats) = grouped_timeseries_clustering(percentage_cover, ones(size(context_layers, 1)); n_clusters=n_clusters, length_t=1:50)

    context_layers[:, "$(GCM)_bioregion_clusters"] = bioregion_clusters
    context_layers[:, "$(GCM)_bioregion_cluster_cats"] = bioregion_cluster_cats
    context_layers[:, "$(GCM)_management_area_clusters"] = man_region_clusters
    context_layers[:, "$(GCM)_management_area_cluster_cats"] = man_region_cluster_cats
    context_layers[:, "$(GCM)_gbr_clusters"] = gbr_clusters
    context_layers[:, "$(GCM)_gbr_cluster_cats"] = gbr_cluster_cats

end

GDF.write("../outputs/ADRIA_results/HighResCoralStress/bellwether_reefs_carbonate.gpkg", context_layers)