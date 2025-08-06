"""
Script to create additional plots required for publication. E.g. maps and GCM plots.
"""

import GeometryOps as GO

include("../../common.jl")

CairoMakie.activate!()

context_layers = GDF.read(joinpath(output_path, "analysis_context_layers_carbonate.gpkg"))

# ECS plot - methods
ecs_values = Dict(
    "EC-Earth3-Veg" => 4.31,
    "ACCESS-ESM1-5" => 3.87,
    "ACCESS-CM2" => 4.72,
    "GFDL-CM4" => 2.9,
    "NorESM2-MM" => 2.5
)
ecs = ecs_plot(collect(values(ecs_values)), [2.5, 5.1], [2.1, 7.7], collect(keys(ecs_values)))
save(joinpath(output_path, "figs/ecs_plot.png"), ecs, pt_per_unit=1)

# GBR map plot - methods
investigation_reefs = context_layers[
    (context_layers.management_area.!="NA").&(context_layers.bioregion.!="NA"),
    :]
regions = GDF.read("../data/GBRMPA_Management_Areas.gpkg")
regions.region_name = ["Mackay/Capricorn", "Cairns/Cooktown", "Far Northern", "Townsville Whitsunday"]
qld = GDF.read("../data/GBRMPA_Reef_Features.gpkg")
qld = qld[qld.FEAT_NAME.=="Mainland", :SHAPE]

map_width = fig_sizes["map_width"]
map_height = fig_sizes["map_height"]
region_col = tuple.(Makie.wong_colors()[1:4], repeat([0.3], 4)) # Manually set alpha value to 0.3

ordered_reefs = sort(investigation_reefs, :bioregion_average_latitude; rev=true)
bioregion_colors_labels = DataFrame(
    bioregion=unique(investigation_reefs.bioregion),
    ref=1:length(unique(investigation_reefs.bioregion)),
    color=distinguishable_colors(length(unique(investigation_reefs.bioregion))),
    label=label_lines.((unique(investigation_reefs.bioregion)); l_length=17)
)

ordered_reefs = leftjoin(ordered_reefs, bioregion_colors_labels; on=:bioregion)
centroids = GO.centroid.(ordered_reefs.geometry)

bgcol = :gray90
fig = Figure(size=(map_width, map_height), fontsize=fontsize, backgroundcolor=bgcol)
ax = Axis(
    fig[1, 1],
    width=map_width - 350,
    height=map_height - 100,
    xgridvisible=false,
    ygridvisible=false,
    limits=((142.5, 154.1), (-25, -10)),
    backgroundcolor=bgcol
)
poly!(qld; color=:darkgray)
poly!(regions.SHAPE, color=region_col)
scatter!(centroids; color=Int64.(ordered_reefs.ref), colormap=unique(ordered_reefs.color), markersize=4)
lines!([(146.25, -20.5), (147.0559, -19.2697)]; color=:black)
text!((146.25, -20.5); text="AIMS - Cape Cleveland", align=(:center, :top))

Legend(
    fig[2, 1],
    [PolyElement(color=col, alpha=0.3) for col in region_col],
    regions.region_name,
    nbanks=2,
    patchsize=(10, 10),
    colgap=8,
    backgroundcolor=bgcol
)
Legend(
    fig[1, 2],
    [MarkerElement(color=bioregion_color, marker=:circle) for bioregion_color in unique(ordered_reefs.color)],
    unique(ordered_reefs.label),
    nbanks=2,
    colgap=6,
    rowgap=1,
    backgroundcolor=bgcol
)
save(joinpath(output_path, "figs/region_map.png"), fig, px_per_unit=dpi)