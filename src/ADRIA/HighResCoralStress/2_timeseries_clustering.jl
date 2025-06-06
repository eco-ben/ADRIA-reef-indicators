using Revise, Infiltrator

using GLMakie, GeoMakie, GraphMakie
using Statistics
using YAXArrays
import GeoDataFrames as GDF

using ADRIA

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

for GCM in dhw_scenarios.dhw.properties["members"]
    
    # Select GCM and load relevant results
    @info "Performing timeseries clustering for $(GCM)"

    rs = ADRIA.load_results("../outputs/ADRIA_results/HighResCoralStress/GBR_2024_10_15_HighResCoralStress__RCPs_45_$(GCM)")
    GCM_results = GCM_analysis_results(rs)
    rel_cover = GCM_results.relative_cover

    rel_cover_less_than_5 = [all(rel_cover[:, i].data .< 5) for i in 1:size(rel_cover, 2)]
    rel_cover_less_than_5_ind = findall(rel_cover_less_than_5)
    # rel_cover_less_than_5 = ones(size(context_layers, 1))

    (bioregion_clusters, bioregion_cluster_cats) = grouped_timeseries_clustering(rel_cover[:, rel_cover_less_than_5], context_layers.bioregion[rel_cover_less_than_5]; n_clusters=n_clusters, length_t=1:50)
    (man_region_clusters, man_region_cluster_cats) = grouped_timeseries_clustering(rel_cover[:, rel_cover_less_than_5], context_layers.management_area[rel_cover_less_than_5]; n_clusters=n_clusters, length_t=1:50)
    (gbr_clusters, gbr_cluster_cats) = grouped_timeseries_clustering(rel_cover[:, rel_cover_less_than_5], ones(sum(rel_cover_less_than_5)); n_clusters=n_clusters, length_t=1:50)

    context_layers[:, "$(GCM)_bioregion_clusters"] = zeros(size(context_layers, 1))
    context_layers[:, "$(GCM)_bioregion_cluster_cats"] = Vector{Union{String, Missing}}(missing, size(context_layers, 1))
    context_layers[:, "$(GCM)_management_area_clusters"] = zeros(size(context_layers, 1))
    context_layers[:, "$(GCM)_management_area_cluster_cats"] = Vector{Union{String, Missing}}(missing, size(context_layers, 1))
    context_layers[:, "$(GCM)_gbr_clusters"] = zeros(size(context_layers, 1))
    context_layers[:, "$(GCM)_gbr_cluster_cats"] = Vector{Union{String, Missing}}(missing, size(context_layers, 1))

    context_layers[rel_cover_less_than_5_ind, "$(GCM)_bioregion_clusters"] = bioregion_clusters
    context_layers[rel_cover_less_than_5_ind, "$(GCM)_bioregion_cluster_cats"] = bioregion_cluster_cats
    context_layers[rel_cover_less_than_5_ind, "$(GCM)_management_area_clusters"] = man_region_clusters
    context_layers[rel_cover_less_than_5_ind, "$(GCM)_management_area_cluster_cats"] = man_region_cluster_cats
    context_layers[rel_cover_less_than_5_ind, "$(GCM)_gbr_clusters"] = gbr_clusters
    context_layers[rel_cover_less_than_5_ind, "$(GCM)_gbr_cluster_cats"] = gbr_cluster_cats
end

GDF.write("../outputs/ADRIA_results/HighResCoralStress/bellwether_reefs.gpkg", context_layers)