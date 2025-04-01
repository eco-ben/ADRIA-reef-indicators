using Revise, Infiltrator

using GLMakie, GeoMakie, GraphMakie
using Statistics
using YAXArrays
using CairoMakie

using ADRIA

include("../../common.jl")
includet("../../plotting_functions.jl")

CairoMakie.activate!()

# Select the target GCM from
dhw_scenarios = open_dataset("../../ADRIA Domains/GBR_2024_10_15_HighResCoralStress/DHWs/dhwRCP45.nc")
dhw_scenarios.dhw.properties["members"]

context_layers = GDF.read("../outputs/ADRIA_results/HighResCoralStress/analysis_context_layers.gpkg")
context_layers.gbr .= "Great Barrier Reef"

gbr_dom = ADRIA.load_domain("../../ADRIA Domains/GBR_2024_10_15_HighResCoralStress/", "45")

for (i_gcm, GCM) in enumerate(dhw_scenarios.dhw.properties["members"])
    # Select GCM and load relevant results
    @info "analysing reef clustering for $(GCM)"

    fig_out_dir = "../outputs/ADRIA_results/HighResCoralStress/figs/$(GCM)"

    rs = ADRIA.load_results("../outputs/ADRIA_results/HighResCoralStress/GBR_2024_10_15_HighResCoralStress__RCPs_45_$(GCM)")
    GCM_results = GCM_analysis_results(rs)
    rel_cover = GCM_results.relative_cover

    rel_cover_less_than_5 = [all(rel_cover[:, i].data .< 5) for i in 1:size(rel_cover, 2)]
    rel_cover_less_than_5_ind = findall(rel_cover_less_than_5)

    dhw_ts = gbr_dom.dhw_scens[:,:,i_gcm]
    ax = (
        Dim{:timesteps}(2025:2099),
        Dim{:locations}(collect(GCM_results.relative_cover.locations))
    )
    dhw_ts = rebuild(dhw_ts, dims=ax)

    # Subset layers to remove reefs that start with very low coral cover
    analysis_layers = context_layers[rel_cover_less_than_5, :]
    rel_cover = rel_cover[:, rel_cover_less_than_5]

    # Remove reefs that are not assigned to a management region or bioregion
    analysis_layers = analysis_layers[(analysis_layers.management_area .!= "NA") .& (analysis_layers.bioregion .!= "NA"), :]
    rel_cover = rel_cover[:, rel_cover.locations .∈ [analysis_layers.UNIQUE_ID]]
    dhw_ts = dhw_ts[:, dhw_ts.locations .∈ [rel_cover.locations]]

    # Analysing clusters identified at bioregion level
    cluster_analysis_plots(analysis_layers, rel_cover, dhw_ts, :bioregion, fig_out_dir)
    cluster_analysis_plots(analysis_layers, rel_cover, dhw_ts, :management_area, fig_out_dir)
    cluster_analysis_plots(analysis_layers, rel_cover, dhw_ts, :gbr, fig_out_dir)
end

    

# Add code to calculate coral evenness and produce further figures for evenness analyses:
## context_layers = context_layers[temporal_scores .< 5, :]
## context_layers.temporal_scores = temporal_scores[temporal_scores .< 5]

# location 1
rel_cover_y1 = Float64.(mapslices(relative_site_cover, absolute_cover[:, :], dims=[:timesteps]))
rel_cover_y2 = Float64.(mapslices(relative_site_cover, absolute_cover[2:end, :], dims=[:timesteps]))

fig = Figure()
ax1 = Axis(fig[1,1];
    ylabel = "Absolute coral cover"
)
ax2 = Axis(fig[2,1];
    ylabel = "Relative coral cover from Y1"
)
ax3 = Axis(fig[3,1];
    ylabel = "Relative coral cover from Y3"
)
ax4 = Axis(fig[4,1];
    ylabel = "DHW"
)


series!(ax1, loc_cover[1:40,500:510].data'; solid_color = (:gray, 0.5))
series!(ax2, rel_cover_y1[1:40,500:510].data'; solid_color = (:gray, 0.5))
series!(ax3, rel_cover_y2[1:40,500:510].data'; solid_color = (:gray, 0.5))
series!(ax4, gbr_dom.dhw_scens[1:40,500:510,1].data'; solid_color = (:gray, 0.5))