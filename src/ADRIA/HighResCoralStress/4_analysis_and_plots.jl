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
context_layers = context_layers[context_layers.management_area .!= "NA", :]

gbr_dom = ADRIA.load_domain("../../ADRIA Domains/GBR_2024_10_15_HighResCoralStress/", "45")

for (i_gcm, GCM) in enumerate(dhw_scenarios.dhw.properties["members"])
    # Select GCM and load relevant results
    @info "analysing bellwether reefs for $(GCM)"

    fig_out_dir = "../outputs/ADRIA_results/HighResCoralStress/figs/$(GCM)"

    rs = ADRIA.load_results("../outputs/ADRIA_results/HighResCoralStress/GBR_2024_10_15_HighResCoralStress__RCPs_45_$(GCM)")
    GCM_results = GCM_analysis_results(rs)

    dhw_ts = gbr_dom.dhw_scens[:,:,i_gcm]
    ax = (
        Dim{:timesteps}(2025:2099),
        Dim{:locations}(collect(GCM_results.relative_cover.locations))
    )
    dhw_ts = rebuild(dhw_ts, dims=ax)

    # Analysing bellwether reefs identified at bioregion level
    if sum(context_layers[:, "$(GCM)_bellwether_reefs_bior"]) > 1
        # Filter to only include reefs that are within the same bioregion/closest_port subregion as target reefs
        filtered_bior = context_layers[context_layers.bioregion .∈ [unique(context_layers[context_layers[:, "$(GCM)_bellwether_reefs_bior"], :bioregion])], :]
        filtered_bior = filtered_bior[filtered_bior.so_to_si .< quantile(filtered_bior.so_to_si, 0.95), :]
        filtered_bior = filtered_bior[filtered_bior.bioregion .∈ [unique(filtered_bior[filtered_bior[:, "$(GCM)_bellwethers_cat_bior"] .== "bellwether", :bioregion])], :]

        removed_bioregions = grouping_counts(
            :bioregion,
            filtered_bior,
            "$(GCM)_bellwether_reefs_bior",
            "../outputs/ADRIA_results/HighResCoralStress/$(GCM)_bioregion_reef_counts.csv",
            GCM
        )

        filtered_bior = filtered_bior[
            (filtered_bior.bioregion .∉ [removed_bioregions]) .&
            (filtered_bior.bioregion .!= "NA"), :
        ]

        @info "$(removed_bioregions) bioregions have been removed from $(GCM) cover analysis as they are lacking bellwether reef numbers"
        
        filtered_bior_bioregions = length(unique(filtered_bior.bioregion))
        # Plotting bioregion violin plots for each variable
        mean_dhw_violin = grouped_violin_plots(
            filtered_bior,
            Symbol("$(GCM)_bellwethers_cat_bior"),
            :bioregion, Symbol("$(GCM)_mean_dhw");
            ylabel="mean DHW", xlabel="Bellwether Reefs"
        );
        save(
            joinpath(fig_out_dir, "mean_dhw_bioregion_violin.png"), 
            mean_dhw_violin
        )

        so_to_si_violin = grouped_violin_plots(
            filtered_bior,
            Symbol("$(GCM)_bellwethers_cat_bior"),
            :bioregion, :so_to_si;
            ylabel="Source to sink ratio", xlabel="Bellwether Reefs"
        );
        save(
            joinpath(fig_out_dir, "so_to_si_bioregion_violin.png"), 
            so_to_si_violin
        )

        total_strength_violin = grouped_violin_plots(
            filtered_bior,
            Symbol("$(GCM)_bellwethers_cat_bior"),
            :bioregion, :total_strength;
            ylabel="Total connectivity strength", xlabel="Bellwether Reefs"
        );
        save(
            joinpath(fig_out_dir, "total_strength_bioregion_violin.png"), 
            total_strength_violin
        )

        initial_cover_violin = grouped_violin_plots(
            filtered_bior,
            Symbol("$(GCM)_bellwethers_cat_bior"),
            :bioregion, Symbol("initial_coral_cover");
            ylabel="Initial coral cover", xlabel="Bellwether Reefs"
        );
        save(
            joinpath(fig_out_dir, "initial_cover_bioregion_violin.png"), 
            initial_cover_violin
        )

        initial_proportion_violin = grouped_violin_plots(
            filtered_bior,
            Symbol("$(GCM)_bellwethers_cat_bior"),
            :bioregion, Symbol("initial_proportion");
            ylabel="Initial proportion of habitable area occupied", xlabel="Bellwether Reefs"
        );
        save(
            joinpath(fig_out_dir, "initial_proportion_bioregion_violin.png"), 
            initial_proportion_violin
        )

        dhw_cover_cor_violin = grouped_violin_plots(
            filtered_bior,
            Symbol("$(GCM)_bellwethers_cat_bior"),
            :bioregion, Symbol("$(GCM)_dhw_cover_cor");
            ylabel="Total coral cover - DHW correlation", xlabel="Bellwether Reefs"
        );
        save(
            joinpath(fig_out_dir, "dhw_cover_cor_bioregion_violin.png"), 
            dhw_cover_cor_violin
        )

        filtered_bior_depth = filtered_bior[filtered_bior.depth_mean .!= 7, :]
        filtered_bior_depth = filtered_bior_depth[filtered_bior_depth.bioregion .∈ [unique(filtered_bior_depth[filtered_bior_depth[:, "$(GCM)_bellwethers_cat_bior"] .== "bellwether", :bioregion])], :]
        depth_median_violin = grouped_violin_plots(
            filtered_bior_depth,
            Symbol("$(GCM)_bellwethers_cat_bior"),
            :bioregion, Symbol("depth_med");
            ylabel="Median Depth (m)", xlabel="Bellwether Reefs"
        );
        save(
            joinpath(fig_out_dir, "depth_cover_bioregion_violin.png"), 
            depth_median_violin
        )

        rel_cover = GCM_results.relative_cover
        rel_cover_reefs = extract_timeseries(rel_cover, filtered_bior, ["bioregion", "$(GCM)_bellwethers_cat_bior", "cor", "management_area"])

        bioregion_grouped_timeseries = grouped_timeseries_plots(
            rel_cover_reefs,
            "$(GCM)_bellwethers_cat_bior",
            :bioregion,
            (1,30), 0;
            ylabel="Total coral cover",
            xlabel="Year"
        )
        save(
            joinpath(fig_out_dir, "bioregion_cover_timeseries.png"), 
            bioregion_grouped_timeseries
        )

        bioregion_grouped_timeseries_lagged = grouped_timeseries_plots(
            rel_cover_reefs,
            "$(GCM)_bellwethers_cat_bior",
            :bioregion,
            (1,30), 4;
            ylabel="Total coral cover",
            xlabel="Year"
        )
        save(
            joinpath(fig_out_dir, "bioregion_cover_lagged_timeseries.png"), 
            bioregion_grouped_timeseries_lagged
        )

        bior_reefs_dhw = extract_timeseries(dhw_ts, filtered_bior, ["bioregion", "$(GCM)_bellwethers_cat_bior", "cor", "management_area"])
        bior_reefs_dhw_plot = grouped_timeseries_plots(
            bior_reefs_dhw,
            "$(GCM)_bellwethers_cat_bior",
            :bioregion,
            (1,30), 0;
            ylabel="Degree Heating Weeks (°C)",
            xlabel="Year"
        )
        save(
            joinpath(fig_out_dir, "bioregion_dhw_timeseries.png"), 
            bior_reefs_dhw_plot
        )
    end

    # Analysing bellwether reefs identified at management_area level
    if sum(context_layers[:, "$(GCM)_bellwether_reefs_man_area"]) > 1
        # Filter to only include reefs that are within the same man_areaegion/closest_port subregion as target reefs
        filtered_man_area = context_layers[context_layers.management_area .∈ [unique(context_layers[context_layers[:, "$(GCM)_bellwether_reefs_man_area"], :management_area])], :]
        filtered_man_area = filtered_man_area[filtered_man_area.so_to_si .< quantile(filtered_man_area.so_to_si, 0.95), :]
        filtered_man_area = filtered_man_area[filtered_man_area.management_area .∈ [unique(filtered_man_area[filtered_man_area[:, "$(GCM)_bellwethers_cat_man_area"] .== "bellwether", :management_area])], :]

        removed_management_areas = grouping_counts(
            :management_area,
            filtered_man_area,
            "$(GCM)_bellwether_reefs_man_area",
            "../outputs/ADRIA_results/HighResCoralStress/$(GCM)_management_area_reef_counts.csv",
            GCM
        )
        filtered_man_area = filtered_man_area[filtered_man_area.management_area .∉ [removed_management_areas], :]

        @info "$(removed_management_areas) management_areas have been removed from $(GCM) cover analysis as they are lacking bellwether reef numbers"
        
        filtered_man_area_management_areas = length(unique(filtered_man_area.management_area))
        # Plotting management_area violin plots for each variable
        mean_dhw_violin = grouped_violin_plots(
            filtered_man_area,
            Symbol("$(GCM)_bellwethers_cat_man_area"),
            :management_area, Symbol("$(GCM)_mean_dhw");
            ylabel="mean DHW", xlabel="Bellwether Reefs"
        );
        save(
            joinpath(fig_out_dir, "mean_dhw_management_area_violin.png"), 
            mean_dhw_violin
        )

        so_to_si_violin = grouped_violin_plots(
            filtered_man_area,
            Symbol("$(GCM)_bellwethers_cat_man_area"),
            :management_area, :so_to_si;
            ylabel="Source to sink ratio", xlabel="Bellwether Reefs"
        );
        save(
            joinpath(fig_out_dir, "so_to_si_management_area_violin.png"), 
            so_to_si_violin
        )

        total_strength_violin = grouped_violin_plots(
            filtered_man_area,
            Symbol("$(GCM)_bellwethers_cat_man_area"),
            :management_area, :total_strength;
            ylabel="Total connectivity strength", xlabel="Bellwether Reefs"
        );
        save(
            joinpath(fig_out_dir, "total_strength_management_area_violin.png"), 
            total_strength_violin
        )

        initial_cover_violin = grouped_violin_plots(
            filtered_man_area,
            Symbol("$(GCM)_bellwethers_cat_man_area"),
            :management_area, Symbol("initial_coral_cover");
            ylabel="Initial coral cover", xlabel="Bellwether Reefs"
        );
        save(
            joinpath(fig_out_dir, "initial_cover_management_area_violin.png"), 
            initial_cover_violin
        )

        initial_proportion_violin = grouped_violin_plots(
            filtered_man_area,
            Symbol("$(GCM)_bellwethers_cat_man_area"),
            :management_area, Symbol("initial_proportion");
            ylabel="Initial proportion of habitable area occupied", xlabel="Bellwether Reefs"
        );
        save(
            joinpath(fig_out_dir, "initial_proportion_management_area_violin.png"), 
            initial_proportion_violin
        )

        dhw_cover_cor_violin = grouped_violin_plots(
            filtered_man_area,
            Symbol("$(GCM)_bellwethers_cat_man_area"),
            :management_area, Symbol("$(GCM)_dhw_cover_cor");
            ylabel="Total coral cover - DHW correlation", xlabel="Bellwether Reefs"
        );
        save(
            joinpath(fig_out_dir, "dhw_cover_cor_management_area_violin.png"), 
            dhw_cover_cor_violin
        )

        filtered_man_area_depth = filtered_man_area[filtered_man_area.depth_mean .!= 7, :]
        filtered_man_area_depth = filtered_man_area_depth[filtered_man_area_depth.management_area .∈ [unique(filtered_man_area_depth[filtered_man_area_depth[:, "$(GCM)_bellwethers_cat_bior"] .== "bellwether", :management_area])], :]
        depth_median_violin = grouped_violin_plots(
            filtered_man_area_depth,
            Symbol("$(GCM)_bellwethers_cat_man_area"),
            :management_area, Symbol("depth_med");
            ylabel="Median Depth (m)", xlabel="Bellwether Reefs"
        );
        save(
            joinpath(fig_out_dir, "depth_cover_man_area_violin.png"), 
            depth_median_violin
        )

        rel_cover = GCM_results.relative_cover
        rel_cover_reefs = extract_timeseries(rel_cover, filtered_man_area, ["management_area", "$(GCM)_bellwethers_cat_man_area", "cor"])

        management_area_grouped_timeseries = grouped_timeseries_plots(
            rel_cover_reefs,
            "$(GCM)_bellwethers_cat_man_area",
            :management_area,
            (1,30), 0;
            ylabel="Total coral cover",
            xlabel="Year"
        )
        save(
            joinpath(fig_out_dir, "management_area_cover_timeseries.png"), 
            management_area_grouped_timeseries
        )

        management_area_grouped_timeseries_lagged = grouped_timeseries_plots(
            rel_cover_reefs,
            "$(GCM)_bellwethers_cat_man_area",
            :management_area,
            (1,30), 4;
            ylabel="Total coral cover",
            xlabel="Year"
        )
        save(
            joinpath(fig_out_dir, "management_area_cover_lagged_timeseries.png"), 
            management_area_grouped_timeseries_lagged
        )

        man_area_reefs_dhw = extract_timeseries(dhw_ts, filtered_man_area, ["management_area", "$(GCM)_bellwethers_cat_man_area", "cor"])
        man_area_reefs_dhw_plot = grouped_timeseries_plots(
            man_area_reefs_dhw,
            "$(GCM)_bellwethers_cat_man_area",
            :management_area,
            (1,30), 0;
            ylabel="Degree Heating Weeks (°C)",
            xlabel="Year"
        )
        save(
            joinpath(fig_out_dir, "management_area_dhw_timeseries.png"), 
            man_area_reefs_dhw_plot
        )
    end
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