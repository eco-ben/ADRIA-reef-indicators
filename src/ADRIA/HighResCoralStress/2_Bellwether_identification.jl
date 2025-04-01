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

for GCM in dhw_scenarios.dhw.properties["members"]
    # Select GCM and load relevant results
    @info "identifying bellwether reefs for $(GCM)"

    rs = ADRIA.load_results("../outputs/ADRIA_results/HighResCoralStress/GBR_2024_10_15_HighResCoralStress__RCPs_45_$(GCM)")
    GCM_results = GCM_analysis_results(rs)
    rel_cover = GCM_results.relative_cover

    # Calculate the split timeseries and temporal variability metrics for the relative coral cover timeseries
    temporal_scores = [temporal_variability(naive_split_metric(@view(rel_cover[1:30, i]), 3, mean)) for i in 1:length(rel_cover.locations)]
    context_layers.temporal_scores = temporal_scores

    rm_cover_temporal_outliers = temporal_scores .< 5

    context_layers[:, "$(GCM)_bellwether_reefs_bior"] .= 0
    context_layers[:, "$(GCM)_bellwethers_cat_bior"] .= ""
    context_layers[:, "$(GCM)_bellwether_reefs_man_area"] .= 0
    context_layers[:, "$(GCM)_bellwethers_cat_man_area"] .= ""

    context_layers[!, "$(GCM)_bellwether_reefs_bior"] = cluster_temporal_variability(context_layers.bioregion, context_layers.temporal_scores; quantile_threshold=0.1)
    context_layers[!, "$(GCM)_bellwethers_cat_bior"] = ifelse.(context_layers[:, "$(GCM)_bellwether_reefs_bior"], "bellwether", "non-bellwether")

    context_layers[!, "$(GCM)_bellwether_reefs_man_area"] = cluster_temporal_variability(context_layers.management_area, context_layers.temporal_scores, quantile_threshold=0.1)
    context_layers[!, "$(GCM)_bellwethers_cat_man_area"] = ifelse.(context_layers[:, "$(GCM)_bellwether_reefs_man_area"], "bellwether", "non-bellwether")

end
GDF.write("../outputs/ADRIA_results/HighResCoralStress/bellwether_reefs.gpkg", context_layers)
