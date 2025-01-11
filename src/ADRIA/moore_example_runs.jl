using Revise, Infiltrator

using GLMakie, GeoMakie, GraphMakie

using ADRIA

moore_domain_path = "c:/Users/grier/Documents/Projects/ADRIA Domains/Moore_2024-08-09_v070_rc1/"

# GBR wide domain
moore_dom = ADRIA.load_domain(moore_domain_path, "45")
context_layers = moore_dom.loc_data

# generate 4096 sample scenarios from counterfactual scenarios
scens = ADRIA.sample_cf(moore_dom, 512)

# repeat(scens, 16)

# scens_repeat = ADRIA.param_table(gbr_dom)
# scens_repeat = repeat(scens_repeat, 5)
# scens_repeat[!, :dhw_scenario] .= [1, 2, 3, 4, 5]
# rs = ADRIA.run_scenarios(gbr_dom, scens_repeat, "45")

# Run sampled scenarios for a given RCP
rs = ADRIA.run_scenarios(moore_dom, scens, "45")

# Extract metric from scenarios
tac = ADRIA.metrics.total_absolute_cover(rs)

# Get a timeseries summarizing the scenarios for each site
absolute_cover = ADRIA.metrics.loc_trajectory(median, tac)

include("../common.jl")
includet("../plotting_functions.jl")
rel_cover = Float64.(mapslices(relative_site_cover, absolute_cover, dims=[:timesteps]))

# Filter out reefs with less than 5% coral cover in initial timesteps
# icc = absolute_cover[1,:].data
# initial_cover_indices = icc .> 5

# Calculate the split timeseries and temporal variability metrics for the relative coral cover timeseries
temporal_scores = [temporal_variability(naive_split_metric(@view(rel_cover[:, i]), 3, mean)) for i in 1:length(rel_cover.locations)]
hist(temporal_scores)

# Plot lines and the distribution of temporal variability scores
rel_cover_vec = [vec(rel_cover.data[:, i]) for i in 1:length(rel_cover.locations)]
sub_2_5 = (temporal_scores .< 2.5)

line_scores_zip = zip(rel_cover_vec[sub_2_5], temporal_scores[sub_2_5])
extremes = extrema(temporal_scores[sub_2_5])
lines_color((a, c)) = lines!(a; color=c, colorrange=extremes, alpha = 0.3)

fig = Figure()
Axis(fig[1,1])
map(lines_color, line_scores_zip)
Colorbar(fig[1,2], limits = extremes)

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

context_layers.lag_0_correlation = cluster_correlation(context_layers.Reef, rel_cover, 0, pearsons_cor)
context_layers.lag_4_correlation =  cluster_correlation(context_layers.Reef, rel_cover, 4, pearsons_cor)
diff_correlation = context_layers.lag_4_correlation .- context_layers.lag_0_correlation
context_layers.bellwether_reefs = diff_correlation .> 0
context_layers.bellwethers_cat = ifelse.(context_layers.bellwether_reefs, "bellwether", "non-bellwether")

context_layers.moore .= "moore"
bellwethers_moore = cluster_temporal_variability(context_layers.moore, context_layers.temporal_scores, 0.1)
bellwethers_reef = cluster_temporal_variability(context_layers.reef, context_layers.temporal_scores, 0.1)

consistent_bellwether_reefs = bellwethers_gbr .& bellwethers_management .& bellwethers_bioregion
context_layers.consistent_bellwethers = consistent_bellwether_reefs
context_layers.consistent_bellwethers_cat = ifelse.(context_layers.consistent_bellwethers, "bellwether", "non-bellwether")


consistent_reefs_dhw = extract_timeseries(dhw_ts, context_layers, ["bioregion", "$(GCM)_lag4_bior", "consistent_bellwethers_cat", "management_area", "gbr"])
consistent_reefs_rel_cover = extract_timeseries(rel_cover, context_layers, ["Reef", "lag_4_correlation", "bellwethers_cat", "management_area"])
consistent_reefs_absolute_cover = extract_timeseries(cover_ts, context_layers, ["bioregion", "$(GCM)_lag4_bior", "consistent_bellwethers_cat", "management_area", "gbr"])

subset_reef_sites = consistent_reefs_rel_cover[consistent_reefs_rel_cover.Reef .∈ [context_layers[context_layers.bellwether_reefs, :Reef]], :]

consistent_reefs_dhw_plot = grouped_timeseries_plots(
    subset_reef_sites,
    "bellwethers_cat",
    "lag_4_correlation",
    :Reef,
    (1,70), 0;
    ylabel="Relative_site_cover",
    xlabel="Year",
    x_fig_size = 1500,
    y_fig_size = 1000
)

