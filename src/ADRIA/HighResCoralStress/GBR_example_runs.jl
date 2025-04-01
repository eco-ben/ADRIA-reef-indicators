using Revise, Infiltrator

using GLMakie, GeoMakie, GraphMakie
using Statistics
using YAXArrays

using ADRIA

GBR_domain_path = "../../ADRIA Domains/GBR_2024_10_15_HighResCoralStress/"
dhw_scenarios = open_dataset("../../ADRIA Domains/GBR_2024_10_15_HighResCoralStress/DHWs/dhwRCP45.nc")
gcms = dhw_scenarios.dhw.properties["members"]

# GBR wide domain
gbr_dom = ADRIA.load_domain(GBR_domain_path, "45")
context_layers = gbr_dom.loc_data

# generate 4096 sample scenarios from counterfactual scenarios
scens = ADRIA.sample_cf(gbr_dom, 4096)

# change dhw scenarios to GFDL-CM4
scens[!, :dhw_scenario] .= 5

# repeat(scens, 16)

# scens_repeat = ADRIA.param_table(gbr_dom)
# scens_repeat = repeat(scens_repeat, 5)
# scens_repeat[!, :dhw_scenario] .= [1, 2, 3, 4, 5]
# rs = ADRIA.run_scenarios(gbr_dom, scens_repeat, "45")

# Run sampled scenarios for a given RCP
#rs = ADRIA.run_scenarios(gbr_dom, scens, "45")
rs = ADRIA.load_results("Outputs/GBR_2024_10_15_HighResCoralStress__RCPs_45_NorESM2-MM")

# Extract metric from scenarios
tac = ADRIA.metrics.total_absolute_cover(rs)

# Get a timeseries summarizing the scenarios for each site
absolute_cover = ADRIA.metrics.loc_trajectory(median, tac)

include("../common.jl")
includet("../plotting_functions.jl")
rel_cover = Float64.(mapslices(relative_site_cover, absolute_cover[:, :], dims=[:timesteps]))

# Filter out reefs with less than 5% coral cover in initial timesteps
# icc = absolute_cover[1,:].data
# initial_cover_indices = icc .> 5

# Calculate the split timeseries and temporal variability metrics for the relative coral cover timeseries
temporal_scores = [temporal_variability(naive_split_metric(@view(rel_cover[1:50, i]), 3, mean)) for i in 1:length(rel_cover.locations)]
context_layers.temporal_scores = temporal_scores
hist(temporal_scores)

# Plot lines and the distribution of temporal variability scores
rel_cover_vec = [vec(rel_cover.data[1:50, i]) for i in 1:length(rel_cover.locations)]
sub_2_5 = (temporal_scores .< 2.5)
sub_2_5_rel_cover_sub_5 = (rel_cover.data[10, :] .< 2) .& sub_2_5

line_scores_zip = zip(rel_cover_vec[sub_2_5_rel_cover_sub_5], temporal_scores[sub_2_5_rel_cover_sub_5])
extremes = extrema(temporal_scores[sub_2_5_rel_cover_sub_5])
lines_color((a, c)) = lines!(ax1, a; color=c, colorrange=extremes, alpha = 0.1)

fig = Figure()
ax1=Axis(fig[1,1], ylabel="Relative Coral Cover")
ax2=Axis(fig[2,1], ylabel="DHW")

map(lines_color, line_scores_zip)
Colorbar(fig[1,2], limits = extremes)
series!(ax2, gbr_dom.dhw_scens[:,:,1].data'; solid_color=(:gray, 0.2))

dhw_lines = [vec(gbr_dom.dhw_scens.data[:, i,1]) for i in 1:length(rel_cover.locations)]
line_dhw_zip = zip(dhw_lines[sub_2_5_rel_cover_sub_5], (1:3806)[sub_2_5_rel_cover_sub_5])
extremes = extrema((1:3806)[sub_2_5_rel_cover_sub_5])


"""

"""
function cluster_temporal_variability(reef_clusters, temporal_variability_results, threshold)
    bellwether_reefs = zeros(Bool, length(reef_clusters))

    for (ind, cluster) in enumerate(reef_clusters)
        indices_not_target = findall(reef_clusters .== cluster)
        indices_not_target = indices_not_target[indices_not_target .!= ind]

        target_temporal_variability = temporal_variability_results[ind]
        non_target_temp_variability = median(temporal_variability_results[indices_not_target])

        bellwether_reefs[ind] = (non_target_temp_variability - target_temporal_variability) > threshold
    end

    return bellwether_reefs
end

context_layers = context_layers[temporal_scores .< 5, :]
context_layers.temporal_scores = temporal_scores[temporal_scores .< 5]

context_layers.lag_0_correlation = cluster_correlation(context_layers.bioregion, rel_cover, 0, pearsons_cor)
context_layers.lag_4_correlation =  cluster_correlation(context_layers.bioregion, rel_cover, 4, pearsons_cor)
diff_correlation = context_layers.lag_4_correlation .- context_layers.lag_0_correlation
context_layers.bellwether_reefs = diff_correlation .> 0
context_layers.bellwethers_cat = ifelse.(context_layers.bellwether_reefs, "bellwether", "non-bellwether")

#context_layers.moore .= "moore"
context_layers.bellwether_reefs = cluster_temporal_variability(context_layers.bioregion, context_layers.temporal_scores, 0.075)
context_layers.bellwethers_cat = ifelse.(context_layers.bellwether_reefs, "bellwether", "non-bellwether")

# consistent_bellwether_reefs = bellwethers_gbr .& bellwethers_management .& bellwethers_bioregion
# context_layers.consistent_bellwethers = consistent_bellwether_reefs
# context_layers.consistent_bellwethers_cat = ifelse.(context_layers.consistent_bellwethers, "bellwether", "non-bellwether")

dhw_ts = Float64.(mapslices(median, gbr_dom.dhw_scens; dims=[:scenarios]))
dhw_ts = dhw_ts[sites = temporal_scores .< 5]
consistent_reefs_dhw = extract_timeseries(dhw_ts, context_layers, ["bioregion", "temporal_scores", "bellwether_reefs", "bellwethers_cat","management_area"])
consistent_reefs_rel_cover = extract_timeseries(rel_cover, context_layers, ["bioregion", "temporal_scores","bellwether_reefs", "bellwethers_cat", "management_area"])
#consistent_reefs_absolute_cover = extract_timeseries(cover_ts, context_layers, ["bioregion", "diff_correlation", "bellwether_reefs", "management_area"])

subset_reef_sites = consistent_reefs_rel_cover[consistent_reefs_rel_cover.bioregion .∈ [context_layers[context_layers.bellwether_reefs, :bioregion]], :]
high_cover_reefs = subset_reef_sites[:, "5"] .< 1.5
subset_reef_sites = subset_reef_sites[high_cover_reefs, :]
subset_reef_sites = subset_reef_sites[subset_reef_sites.bioregion .∈ [subset_reef_sites[subset_reef_sites.bellwether_reefs, :bioregion]], :]


consistent_reefs_dhw_plot = grouped_timeseries_plots(
    subset_reef_sites,
    "bellwethers_cat",
    "temporal_scores",
    :bioregion,
    (1,50), 0;
    ylabel="Relative_site_cover",
    xlabel="Year",
    x_fig_size = 800,
    y_fig_size = 600
)


fig = Figure()

ax1 = Axis(
    fig[1,1],
    ylabel="absolute_reef_cover",
    title="Reef 1 $(context_layers.GBRMPA_ID[1])"
)
lines!(ax1, 2025:2099, absolute_cover[:, 1].data)

ax2 = Axis(
    fig[2,1],
    ylabel="Relative_reef_cover"
)
lines!(ax2, 2025:2099, rel_cover[:, 1].data)


ax3 = Axis(
    fig[3,1],
    ylabel="DHW"
)
lines!(ax3, 2025:2099, dhw_ts[:, 1].data)
save("Outputs/cover_and_dhw_reef_6.png", fig)