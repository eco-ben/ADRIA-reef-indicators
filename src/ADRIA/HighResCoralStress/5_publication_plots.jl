"""
Script to create additional plots required for publication. E.g. maps and GCM plots.
"""

include("../../common.jl")

CairoMakie.activate!()

context_layers = GDF.read(joinpath(output_path, "analysis_context_layers_carbonate.gpkg"))

# ECS plot - methods
ecs_values = Dict(
    "EC-Earth3-Veg" => 4.30,  # Wyser et al. (2020). On the increased climate sensitivity in the EC-Earth model from CMIP5 to CMIP6
    "ACCESS-ESM1-5" => 3.87,  # Ziehn et al. (2020). The Australian Earth System Model: ACCESS-ESM1.5
    "ACCESS-CM2" => 4.66,  # Grose et al. (2020). A CMIP6-based multi-model downscaling ensemble to underpin climate change services in Australia
    "GFDL-CM4" => 4.1,  # Sentman et al. (2025). Quantifying Equilibrium Climate Sensitivity to Atmospheric Chemistry and Composition Representations in GFDL-CM4.0 and GFDL-ESM4.1
    "NorESM2-MM" => 2.5  # Seland et al. (2020). Overview of the Norwegian Earth System Model (NorESM2) and key climate response of CMIP6 DECK, historical, and scenario simulations
)
ecs = ecs_plot(collect(values(ecs_values)), [2.5, 5.1], [2.1, 7.7], collect(keys(ecs_values)))
save(joinpath(figs_path, "ecs_plot.png"), ecs, px_per_unit=dpi)

# GBR map plot - methods
bioregion_colors = distinguishable_colors(length(unique(context_layers.bioregion)));
gbr_methods_map = map_gbr_reefs(context_layers, :bioregion, bioregion_colors, "Bioregions")

save(joinpath(figs_path, "region_map.png"), gbr_methods_map, px_per_unit=dpi)

# GBR - wide DHW figure
dhw_scenarios = open_dataset(joinpath(gbr_domain_path, "DHWs/dhwRCP45.nc"))
GCMs = dhw_scenarios.dhw.properties["members"]

context_layers = GDF.read(joinpath(output_path, "analysis_context_layers_carbonate.gpkg"))
gbr_dom = ADRIA.load_domain(gbr_domain_path, "45")
gbr_dom_filtered = gbr_dom.loc_data[gbr_dom.loc_data.UNIQUE_ID.∈[context_layers.UNIQUE_ID], :]
filtered_indices = indexin(gbr_dom_filtered.UNIQUE_ID, gbr_dom.loc_data.UNIQUE_ID)

dhw_arrays = [
    rebuild(gbr_dom.dhw_scens[:, filtered_indices, i_gcm],
        dims=(gbr_dom.dhw_scens.timesteps, Dim{:sites}(context_layers.UNIQUE_ID))
    ) for i_gcm in eachindex(GCMs)
]
dhw_arrays = concatenatecubes(dhw_arrays, Dim{:GCM}(GCMs))
dhw_arrays = rebuild(dhw_arrays, metadata=dhw_timeseries_properties)

context_layers = sort(context_layers, :management_area)

dhw_arrays = dhw_arrays[sites=At(context_layers.UNIQUE_ID)]

fig = Figure(
    size=(fig_sizes["gcm_timeseries_width"], fig_sizes["gcm_timeseries_height"] * 1.3),
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
            :xticks => (
                [1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50],
                string.([2025, 2030, 2035, 2040, 2045, 2050, 2055, 2060, 2065, 2070, 2075])
            ),
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
save(joinpath(figs_path, "methods_dhw_timeseries.png"), fig, px_per_unit=dpi)
