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
save(joinpath(output_path, "figs/ecs_plot.png"), ecs, px_per_unit=dpi)

# GBR map plot - methods
investigation_reefs = context_layers
regions = GDF.read("../data/GBRMPA_Management_Areas.gpkg")
regions.region_name = ["Mackay/Capricorn", "Cairns/Cooktown", "Far Northern", "Townsville Whitsunday"]
qld = GDF.read("../data/GBRMPA_Reef_Features.gpkg")
qld = qld[qld.FEAT_NAME.=="Mainland", :SHAPE]

map_width = fig_sizes["map_width"]
map_height = fig_sizes["map_height"]
region_col = tuple.(Makie.wong_colors()[1:4], fill(0.3, 4)) # Manually set alpha value to 0.3

ordered_reefs = sort(investigation_reefs, :bioregion_average_latitude; rev=true)

reef_colors = distinguishable_colors(length(unique(investigation_reefs.bioregion)))
reef_colors = tuple.(reef_colors, fill(0.4, length(reef_colors)))

bioregion_colors_labels = DataFrame(
    bioregion=unique(investigation_reefs.bioregion),
    ref=1:length(unique(investigation_reefs.bioregion)),
    color=reef_colors,
    label=label_lines.((unique(investigation_reefs.bioregion)); l_length=17)
)

ordered_reefs = leftjoin(ordered_reefs, bioregion_colors_labels; on=:bioregion)
centroids = GO.centroid.(ordered_reefs.geometry)

# GeoMakie currently has a bug where x/y labels never display.
# We adjust map size to roughly align with the correct projection.
bgcol = :gray90
fig = Figure(size=(map_width + 75, map_height), fontsize=fontsize, backgroundcolor=bgcol)
ax = Axis(
    fig[1, 1],
    width=map_width - 275,
    height=map_height - 120,
    xgridvisible=false,
    ygridvisible=false,
    limits=((142.5, 154.1), (-25, -10)),
    backgroundcolor=bgcol
)
poly!(qld; color=:darkgray)
poly!(regions.SHAPE, color=region_col)

scatter!(centroids; color=ordered_reefs.color, markersize=4)
ax.ylabel = "Latitude"
ax.xlabel = "Longitude"

# Note key cities/towns for reference
scatter!((145.754120, -16.925491); color=:black)
text!((145.1, -16.925491); text="Cairns", align=(:center, :top))

scatter!((146.8057, -19.2664); color=:black)
text!((146.0, -19.2664); text="Townsville", align=(:center, :top))

scatter!((149.182147, -21.142496); color=:black)
text!((148.432147, -21.142496); text="Mackay", align=(:center, :top))

scatter!((150.733333, -23.133333); color=:black)
text!((150.0, -23.133333); text="Yeppoon", align=(:center, :top))

# lines!([(146.25, -20.5), (147.0559, -19.2697)]; color=:black)
# text!((146.25, -20.5); text="AIMS - Cape Cleveland", align=(:center, :top))

Legend(
    fig[2, 1],
    [PolyElement(color=col, alpha=0.3) for col in region_col],
    regions.region_name,
    "Management regions",
    nbanks=2,
    patchsize=(10, 10),
    colgap=8,
    backgroundcolor=bgcol
)
Legend(
    fig[1, 2],
    [MarkerElement(color=bioregion_color, marker=:circle) for bioregion_color in unique(ordered_reefs.color)],
    unique(ordered_reefs.label),
    "Bioregions",
    nbanks=2,
    colgap=6,
    rowgap=1,
    backgroundcolor=bgcol
)
save(joinpath(output_path, "figs/region_map.png"), fig, px_per_unit=dpi)

# GBR - wide DHW figure
dhw_scenarios = open_dataset(joinpath(gbr_domain_path, "DHWs/dhwRCP45.nc"))
GCMs = dhw_scenarios.dhw.properties["members"]

context_layers = GDF.read(joinpath(output_path, "analysis_context_layers_carbonate.gpkg"))
gbr_dom = ADRIA.load_domain(gbr_domain_path, "45")
gbr_dom_filtered = gbr_dom.loc_data[gbr_dom.loc_data.UNIQUE_ID.∈[context_layers.UNIQUE_ID], :]
filtered_indices = indexin(gbr_dom_filtered.UNIQUE_ID, gbr_dom.loc_data.UNIQUE_ID)

dhw_arrays = [
    rebuild(gbr_dom.dhw_scens[:, filtered_indices, i_gcm], dims=(gbr_dom.dhw_scens.timesteps, Dim{:sites}(context_layers.UNIQUE_ID))) for i_gcm in eachindex(GCMs)
]
dhw_arrays = concatenatecubes(dhw_arrays, Dim{:GCM}(GCMs))
dhw_arrays = rebuild(dhw_arrays, metadata=dhw_timeseries_properties)

context_layers = sort(context_layers, :management_area)

dhw_arrays = dhw_arrays[sites=At(context_layers.UNIQUE_ID)]

fig = Figure(
    size=(fig_sizes["gcm_timeseries_width"], fig_sizes["gcm_timeseries_height"]),
    fontsize=fontsize
)
plot_layout = [(x, 1) for x in eachindex(GCMs)]
man_areas_categorical = CategoricalArray(context_layers.management_area)
for (g, GCM) in enumerate(GCMs)
    group_timeseries = dhw_arrays[1:50, :, g]
    plot_layout_xi = plot_layout[g]
    # ax = _setup_grouped_axes(
    #         fig,
    #         plot_layout_xi,
    #         ([1, 10, 20, 30, 40, 50], string.([2025, 2035, 2045, 2055, 2065, 2075]));
    #         ylabel="DHW [\u00B0 C]",
    #         xlabel="",
    #         title="",
    #         xsize=nothing,
    #         ysize=nothing,
    #         background_color=:white,
    #         xticklabelrotation=45.0,
    #         spinewidth=0.5
    # )

    ADRIA.viz.clustered_scenarios!(
        fig[plot_layout_xi...],
        group_timeseries,
        Vector{Int64}(man_areas_categorical.refs);
        opts=Dict{Symbol,Any}(:legend => false),
        axis_opts=Dict(
            :title => GCM,
            :xticks => ([1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50], string.([2025, 2030, 2035, 2040, 2045, 2050, 2055, 2060, 2065, 2070, 2075])),
            :spinewidth => 0.5,
            :ylabelpadding => 2,
            :xlabelpadding => 2,
            :xticksize => 2,
            :yticksize => 2,
            # :ylabel => ""
        )
    )
end

axes_above_bottom = filter(x -> x isa Axis, fig.content)[1:end-1]
map(x -> hidexdecorations!(x; grid=false, ticks=false), axes_above_bottom)
map(x -> hideydecorations!(x, ticks=false, ticklabels=false, grid=false), filter(x -> x isa Axis, fig.content))
linkaxes!(filter(x -> x isa Axis, fig.content)...)

legend_entries = []
for (i, col) in enumerate([:green, :orange, :blue, :red])
    LE = LineElement(; color=col, marker=:circle)
    push!(legend_entries, [LE])
end

Legend(
    fig[length(GCMs)+1, 1],
    legend_entries,
    levels(man_areas_categorical),
    nbanks=2,
    tellwidth=false,
    padding=(2.0, 2.0, 2.0, 2.0),             # shrink padding inside legend box
    labelsize=fontsize,    # smaller font
    framevisible=false,        # optional: remove box
    # markerlabelgap = 3,          # reduce space between marker and label
    rowgap=0,                  # reduce vertical spacing between items
    colgap=4,                  # reduce horizontal spacing
    patchsize=(5, 5)
)
rowsize!(fig.layout, maximum(first.(plot_layout)) + 1, Relative(0.03))
feat = dhw_arrays.properties[:metric_feature]
unit_text = dhw_arrays.properties[:metric_unit]
Label(
    fig[:, 0], "$feat [$unit_text]",
    rotation=π / 2,  # Rotate 90 degrees
    tellheight=false,
    tellwidth=true
)
save(joinpath(output_path, "figs/methods_dhw_timeseries.png"), fig, px_per_unit=dpi)
