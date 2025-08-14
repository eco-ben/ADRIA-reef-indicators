
using CSV
using CategoricalArrays
using StatsModels
using Turing, Distributions, StatsFuns, LinearAlgebra
using MLDataUtils: shuffleobs, stratifiedobs, rescale!
using JLSO
using ReverseDiff

include("../../common.jl")

context_layers = GDF.read(joinpath(output_path, "analysis_context_layers_carbonate.gpkg"))
context_layers.log_so_to_si = log10.(context_layers.so_to_si)
context_layers.log_total_strength = log10.(context_layers.total_strength)
context_layers.log_area = log10.(context_layers.area)
context_layers.log_depth = log10.(context_layers.depth_med)

dhw_scenarios = open_dataset(joinpath(gbr_domain_path, "DHWs/dhwRCP45.nc"))
GCMs = dhw_scenarios.dhw.properties["members"]

# Prepare long-form of data to support multinomial analysis.
# Each row should be unique attribute combination of reef properties:
# reef, depth, connectivity, dhw, cluster and GCM

gcm_dhw_cols = [Symbol("$(GCM)_mean_dhw") for GCM in GCMs]
gcm_bioregion_cluster_cols = [Symbol("$(GCM)_bioregion_cluster_cats") for GCM in GCMs]
gcm_dhw_tol_cols = [Symbol("$(GCM)_mean_dhw_tol") for GCM in GCMs]

collated_reef_properties = Vector{DataFrame}(undef, length(GCMs))
target_cols = [:UNIQUE_ID, :GBRMPA_ID, :depth_med, :log_total_strength, :log_so_to_si, :bioregion, :log_area, :log_depth, :initial_proportion]
reef_properties = context_layers[:, target_cols]
for (i, GCM) in enumerate(GCMs)
    dhw_col = Symbol("$(GCM)_mean_dhw")
    bio_cluster_col = Symbol("$(GCM)_bioregion_cluster_cats")
    dhw_tol_col = Symbol("$(GCM)_mean_DHW_tol")

    bio_cluster_details = context_layers[:, [dhw_col, bio_cluster_col, dhw_tol_col]]
    rename!(bio_cluster_details, dhw_col => :mean_dhw, bio_cluster_col => :cluster, dhw_tol_col => :mean_DHW_tol)
    bio_cluster_details.GCM .= GCM
    bio_cluster_details.cluster .= categorical(bio_cluster_details.cluster)

    collated_reef_properties[i] = hcat(
        copy(reef_properties),
        bio_cluster_details
    )
end

reef_properties = vcat(collated_reef_properties...)
reef_properties.GCM .= categorical(reef_properties.GCM)
reef_properties.bioregion .= categorical(reef_properties.bioregion)

# plot check of predictor covariance:
@df reef_properties corrplot(
    cols([:depth_med, :log_total_strength, :log_so_to_si, :mean_dhw, :log_area, :mean_DHW_tol]), 
    grid=false, 
    tickfontsize=6,
    titlefontsize=8,
    labelfontsize=6,
    compact=true,
    size=(500,500)
)

# # lhs = Term(:cluster)
# # rhs_terms = [
# #     Term(:depth_med),
# #     Term(:log_total_strength),
# #     Term(:log_so_to_si),
# #     Term(:mean_dhw),
# #     (1|Term(:GCM)),
# #     (1|Term(:bioregion))
# # ]

# f = FormulaTerm(Term(:cluster), Term(:depth_med) + Term(:log_total_strength) + Term(:log_so_to_si) + Term(:mean_dhw) + (ConstantTerm(1) | Term(:bioregion)) + (ConstantTerm(1) | Term(:GCM)))
# f = FormulaTerm(Term(:cluster), Term(:depth_med) + Term(:log_total_strength) + Term(:log_so_to_si) + Term(:mean_dhw) + (Term(:bioregion))|Term(:GCM))

# fm2 = lmm(f, reef_properties)

# model = fit(
#     EconometricModel,
#     f,
#     reef_properties,
#     contrasts=Dict(:cluster => DummyCoding(base="low"))
# )

# coeftab = DataFrame(coeftable(model))

# regression_save_path = joinpath(
#     output_path,
#     "cluster_regression_outputs/bioregion_cluster_regression_coefficients.csv"
# )
# CSV.write(regression_save_path, coeftab)

# """
#     ordinal_random_intercepts(y, X, grp1, grp2, n_grp1, n_grp2)

# # Arguments
# - `y` : Integer vector of reef cluster assignment levels.
# - `X` : n-samples * n-predictors Matrix of independent variables
# - `grp1` : Integer vector of the first categorical random effects term
# - `grp2` : Integer vector of the second categorical random effects term
# - `n_grp1` : Number of unique levels in grp1
# - `n_grp2` : Number of unique levels in grp2
# """
# @model function ordinal_random_intercepts(y, X::Matrix{Float64}, grp1::Vector{UInt32}, grp2::Vector{UInt32}, n_grp1::Int64, n_grp2::Int64)
#     N, p = size(X)
#     # @assert p == 4

#     # priors for fixed effects
#     # weakly-informative prior on beta
#     # β ~ MvNormal(zeros(p), 5.0^2 * I)
#     β ~ filldist(Normal(0, 5), p)

#     # --- RANDOM EFFECTS: non-centered parameterization ---
#     # grp1
#     σ1 ~ truncated(Cauchy(0, 2), 0, Inf)
#     # z1 ~ filldist(Normal(0,1), n_grp1)
#     u1 ~ filldist(Normal(0, σ1), n_grp1)

#     # grp2
#     σ2 ~ truncated(Cauchy(0, 2), 0, Inf)
#     # z2 ~ filldist(Normal(0,1), n_grp2)
#     u2 ~ filldist(Normal(0, σ2), n_grp2)

#     # # transform to actual random intercepts
#     # u1 = σ1 .* z1
#     # u2 = σ2 .* z2

#     # estimates for cutpoints - required for ordinal bayesian model to separate the 3 levels of response variable
#     c_raw ~ filldist(Normal(0, 5), 2)
#     c = sort(c_raw)  # c[1] < c[2]


#     # Likelihood
#     # for i in 1:N
#     #     η = dot(X[i, :], β) + u1[grp1[i]] + u2[grp2[i]]
#     #     p1 = cdf(Normal(), c[1] - η)
#     #     p2 = cdf(Normal(), c[2] - η)
#     #     ps = [p1, p2 - p1, 1 - p2]  # probabilities for categories 1,2,3
#     #     ps = clamp.(ps, 1e-10, 1.0) # Prevent errors if probabilities appear as 0 or -ve
#     #     ps ./= sum(ps)

#     #     y[i] ~ Distributions.Categorical(ps)
#     # end

#     # Likelihood
#     η = X * β .+ u1[grp1] .+ u2[grp2]
#     p1 = cdf.(Normal(), c[1] .- η)
#     p2 = cdf.(Normal(), c[2] .- η)

#     ps = hcat(p1, p2 .- p1, 1 .- p2)
#     ps = clamp.(ps, 1e-10, 1.0)
#     ps ./= sum(ps, dims=2)

#     for i in 1:N
#         y[i] ~ Distributions.Categorical(vec(ps[i, :]))
#     end
#     # y ~ arraydist([Categorical(ps[i, :]) for i in 1:N])
# end

@model function ordinal_mixed_informed_priors(y, X::Matrix{Float64}, grp1::Vector{UInt32}, grp2::Vector{UInt32}, n_grp1::Int64, n_grp2::Int64, pos_idx, neg_idx)
    N, fp = size(X)
    # @assert p == 4

    other_fp = setdiff(1:fp, (pos_idx, neg_idx))
    β = Vector{Float64}(undef, fp)
    for fp_i in other_fp
        β[fp_i] ~ Normal(0, 5)
    end
    β[pos_idx] ~ truncated(Normal(0, 5), 0, Inf)
    β[neg_idx] ~ truncated(Normal(0, 5), -Inf, 0)

    # --- RANDOM EFFECTS: non-centered parameterization ---
    # grp1
    σ1 ~ truncated(Cauchy(0, 2), 0, Inf)  
    # z1 ~ filldist(Normal(0,1), n_grp1)
    u1 ~ filldist(Normal(0, σ1), n_grp1)

    # grp2
    σ2 ~ truncated(Cauchy(0, 2), 0, Inf)
    # z2 ~ filldist(Normal(0,1), n_grp2)
    u2 ~ filldist(Normal(0, σ2), n_grp2)

    # # transform to actual random intercepts
    # u1 = σ1 .* z1
    # u2 = σ2 .* z2

    # estimates for cutpoints - required for ordinal bayesian model to separate the 3 levels of response variable
    c_raw ~ filldist(Normal(0, 5), 2)
    c = sort(c_raw)  # c[1] < c[2]

    # Likelihood
    η = X * β .+ u1[grp1] .+ u2[grp2]
    p1 = cdf.(Normal(), c[1] .- η)
    p2 = cdf.(Normal(), c[2] .- η)

    ps = hcat(p1, p2 .- p1, 1 .- p2)
    ps = clamp.(ps, 1e-10, 1.0)
    ps ./= sum(ps, dims=2)

    for i in 1:N
        y[i] ~ Distributions.Categorical(vec(ps[i, :]))
    end
    # y ~ arraydist([Categorical(ps[i, :]) for i in 1:N])
end

# @model function ordinal_mixed_informed_priors_1grp(y, X::Matrix{Float64}, grp2::Vector{UInt32}, n_grp2::Int64, pos_idx, neg_idx)
#     N, fp = size(X)
#     # @assert p == 4

#     other_fp = setdiff(1:fp, (pos_idx, neg_idx))
#     β = Vector{Float64}(undef, fp)
#     for fp_i in other_fp
#         β[fp_i] ~ Normal(0, 5)
#     end
#     β[pos_idx] ~ truncated(Normal(0, 5), 0, -Inf)
#     β[neg_idx] ~ truncated(Normal(0, 5), -Inf, 0)

#     # --- RANDOM EFFECTS: non-centered parameterization ---
#     # grp1
#     # σ1 ~ truncated(Cauchy(0, 2), 0, Inf)
#     # # z1 ~ filldist(Normal(0,1), n_grp1)
#     # u1 ~ filldist(Normal(0, σ1), n_grp1)

#     # grp2
#     σ2 ~ truncated(Cauchy(0, 2), 0, Inf)
#     # z2 ~ filldist(Normal(0,1), n_grp2)
#     u2 ~ filldist(Normal(0, σ2), n_grp2)

#     # # transform to actual random intercepts
#     # u1 = σ1 .* z1
#     # u2 = σ2 .* z2

#     # estimates for cutpoints - required for ordinal bayesian model to separate the 3 levels of response variable
#     c_raw ~ filldist(Normal(0, 5), 2)
#     c = sort(c_raw)  # c[1] < c[2]

#     # Likelihood
#     η = X * β .+ u2[grp2]
#     p1 = cdf.(Normal(), c[1] .- η)
#     p2 = cdf.(Normal(), c[2] .- η)

#     ps = hcat(p1, p2 .- p1, 1 .- p2)
#     ps = clamp.(ps, 1e-10, 1.0)
#     ps ./= sum(ps, dims=2)

#     for i in 1:N
#         y[i] ~ Distributions.Categorical(vec(ps[i, :]))
#     end
#     # y ~ arraydist([Categorical(ps[i, :]) for i in 1:N])
# end

# @model function multilevel_model(y, X, Z_rand_slope, grp2, n_grp2, pos_idx, neg_idx)
#   # number of predictors
#     N, fp = size(X)  # including intercept (column of 1s) if desired
#     num_predictors_Z = size(Z_rand_slope, 2)

#     # Prior for standard deviation for errors.
#     other_fp = setdiff(1:fp, (pos_idx, neg_idx))
#     β = Vector{Float64}(undef, fp)
#     for fp_i in other_fp
#         β[fp_i] ~ Normal(0, 5)
#     end
#     β[pos_idx] ~ truncated(Normal(0, 5), 0, Inf)
#     β[neg_idx] ~ truncated(Normal(0, 5), -Inf, 0)

#     # Prior for variance of slope group effects. Usually requires thoughtful specification.
#     s2 ~ truncated(InverseGamma(3, 2), 0.05, 10.0)
#     s = sqrt(s2)
#     z1 ~ filldist(Normal(0,1), num_predictors_Z)
#     u1 = z1 .* s

#     # grp2
#     σ2 ~ truncated(Cauchy(0, 2), 0, Inf)
#     u2 ~ filldist(Normal(0,σ2), n_grp2)

#     # estimates for cutpoints - required for ordinal bayesian model to separate the 3 levels of response variable
#     c1 ~ Normal(0, 5)
#     δ  ~ Truncated(Normal(0, 2), 0.1, Inf)  # minimum 0.1 gap
#     c = (c1, c1 + δ)
#     # c_raw ~ filldist(Normal(0, 5), 2)
#     # c = sort(c_raw)  # c[1] < c[2]

#     # Likelihood
#     η = (X * β) .+ (Z_rand_slope * u1) .+ u2[grp2]
#     p1 = cdf.(Normal(), c[1] .- η)
#     p2 = cdf.(Normal(), c[2] .- η)

#     ps = hcat(p1, p2 .- p1, 1 .- p2)
#     ps = clamp.(ps, 1e-10, 1.0)
#     ps ./= sum(ps, dims=2)

#     for i in 1:N
#         y[i] ~ Distributions.Categorical(vec(ps[i, :]))
#     end
# end

n_GCMs = length(unique(reef_properties.GCM))
n_bioregions = length(unique(reef_properties.bioregion))

# Encode random effects groups (GCM and bioregion) as int
GCM_idx = categorical(reef_properties.GCM).refs
bioregion_idx = categorical(reef_properties.bioregion).refs

predictors = [:depth_med, :log_total_strength, :log_so_to_si, :mean_dhw, :log_area, :mean_DHW_tol, :initial_proportion]

# function split_data(df, target; at=0.70)
#     shuffled = shuffleobs(df)
#     return trainset, testset = stratifiedobs(row -> row[target], shuffled; p=at)
# end

# target = :cluster_vals

# trainset, testset = split_data(reef_properties, target; at=0.50)
for feature in predictors
    μ, σ = rescale!(reef_properties[!, feature]; obsdim=1)
    # rescale!(testset[!, feature], μ, σ; obsdim=1)
end

# random_intercept_design = zeros(nrow(reef_properties), n_bioregions)
# for i in 1:nrow(reef_properties)
#     random_intercept_design[i, bioregion_idx[i]] = 1
# end

# random_slope_depth = random_intercept_design .* reef_properties.depth_med

# Z_random_int_slope_design = hcat(random_intercept_design, random_slope_depth)




# Encode response levels as increasing integers
reef_properties.cluster_vals = categorical(reef_properties.cluster; ordered=true, levels=["low", "medium", "high"]).refs
#reef_properties.row_n = 1:nrow(reef_properties)



# train_idx = indexin(trainset.row_n, reef_properties.row_n)
# test_idx = indexin(testset.row_n, reef_properties.row_n)

# # Turing requires data in matrix form, not dataframe
# train = Matrix(trainset[:, predictors])
# test = Matrix(testset[:, predictors])
# train_cluster = trainset[:, target]
# test_cluster = testset[:, target];
# train_GCM = GCM_idx[train_idx]
# train_bioregion = bioregion_idx[train_idx]
# test_GCM = GCM_idx[test_idx]
# test_bioregion = bioregion_idx[test_idx]

X = Matrix(reef_properties[:, predictors])
clusters = reef_properties.cluster_vals

# model_rand_slopes = multilevel_model(
#     clusters,
#     X,
#     Z_random_int_slope_design,
#     GCM_idx,
#     n_GCMs,
#     findfirst(predictors .== :log_depth),
#     findfirst(predictors .== :mean_dhw)
# )

# rand_slopes_chain = sample(model_rand_slopes, NUTS(0.8; adtype=AutoReverseDiff(; compile=true)), MCMCThreads(), 2000, 4)

model_1_rand_var = ordinal_mixed_informed_priors(
    clusters, 
    X,
    bioregion_idx,
    GCM_idx,
    n_bioregions,
    n_GCMs,
    findfirst(predictors .== :depth_med),
    findfirst(predictors .== :mean_dhw)
)
# model_constrained_priors = ordinal_mixed_informed_priors(
#     clusters,
#     X,
#     bioregion_idx,
#     GCM_idx,
#     n_bioregions,
#     n_GCMs,
#     findfirst(predictors .== :depth_med), # We know that depth is likely to have a positive impact
#     findfirst(predictors .== :mean_dhw) # We know that it is unlikely for DHW to have a positive impact on reefs
# )

# chain = sample(model, NUTS(), MCMCThreads(), 1000, 4) # 4 chains using Threads, 1000 samples each
# JLSO.save(joinpath(output_path, "ordinal_mixed_regression_chain_$(now).jlso"), chain)

# StatsPlots.plot(chain)

small_chain_1_rvar = sample(model_1_rand_var, NUTS(), 200; warmup=100, chains=1)
informed_prior_small_chain = sample(model_constrained_priors, NUTS(), 150; warmup=100, chains=1)

description = describe(small_chain)

summary = leftjoin(DataFrame(summarize(small_chain)), DataFrame(quantile(small_chain)), on=:parameters)

beta_params = summary[1:length(predictors), :]
fig = Figure()
ax = Axis(fig[1,1], ylabel = "Parameter", xlabel = "Estimates mean (2.5-97.5% quantiles)",
    yticks = (1:length(predictors), String.(predictors[[2,3,5,6,1,4]])))
Makie.scatter!(ax, beta_params.mean, 1:length(predictors); color=:blue)
lines_x = [[beta_params[x, "2.5%"], beta_params[x, "97.5%"]] for x in 1:length(predictors)]
lines_y = [[x, x] for x in 1:length(predictors)]

Makie.lines!.(lines_x, lines_y)

function predict_probs(chain, X, grp1, grp2, n_grp1, n_grp2)
    N, p = size(X)

    # Extract posterior draws
    βs = Array(chain[["β[$(i)]" for i in 1:p]])
    u2s = Array(chain[["u2[$(i)]" for i in 1:n_grp2]])

    cs   = Array(chain[["c_raw[1]", "c_raw[2]"]]) # Will need sorting per draw

    n_draws = size(βs, 1)
    probs = zeros(n_draws, N, 3)  # 3 categories

    for d in 1:n_draws
        β   = βs[d, :]
        c   = sort(cs[d, :])
        u2 = u2s[d, :]

        η   = X * β .+ u2[grp2]
        p1  = cdf.(Normal(), c[1] .- η)
        p2  = cdf.(Normal(), c[2] .- η)

        probs[d, :, 1] .= p1
        probs[d, :, 2] .= p2 .- p1
        probs[d, :, 3] .= 1 .- p2
    end
    return probs
end

function predict_probs_slopes(chain, X, Z_rand_slope,  grp2, n_grp2)
    N, p = size(X)

    # Extract posterior draws
    βs = Array(chain[["β[$(i)]" for i in 1:p]])
    u2s = Array(chain[["u2[$(i)]" for i in 1:n_grp2]])

    cs   = Array(chain[["c_raw[1]", "c_raw[2]"]]) # Will need sorting per draw

    n_draws = size(βs, 1)
    probs = zeros(n_draws, N, 3)  # 3 categories

    for d in 1:n_draws
        β   = βs[d, :]
        c   = sort(cs[d, :])
        u2 = u2s[d, :]

        η   = X * β .+ u2[grp2]
        p1  = cdf.(Normal(), c[1] .- η)
        p2  = cdf.(Normal(), c[2] .- η)

        probs[d, :, 1] .= p1
        probs[d, :, 2] .= p2 .- p1
        probs[d, :, 3] .= 1 .- p2
    end
    return probs
end

probs_weak_prior = predict_probs(weak_prior_small_chain, X, bioregion_idx, GCM_idx, n_bioregions, n_GCMs)
pred_cats_weak_prior = mapslices(p -> argmax(p), probs_weak_prior; dims=3)
mean_probs_weak_prior = dropdims(mean(probs_weak_prior, dims=1); dims=1)  # Average over draws
mean_pred_cats_weak_prior = [argmax(mean_probs_weak_prior[i, :]) for i in eachindex(reef_properties.cluster)]

n_draws, N = size(probs_weak_prior)
# Accuracy per draw
acc_per_draw = [mean(pred_cats_weak_prior[d, :] .== reef_properties.cluster_vals) for d in 1:n_draws]
# Summarize accuracy distribution
mean_acc = mean(acc_per_draw)
acc_ci = quantile(acc_per_draw, [0.025, 0.975])


probs_constrained_prior = predict_probs(informed_prior_small_chain, X, bioregion_idx, GCM_idx, n_bioregions, n_GCMs)
pred_cats_constrained_prior = mapslices(p -> argmax(p), probs_constrained_prior; dims=3)
mean_probs_constrained_prior = dropdims(mean(probs_constrained_prior, dims=1); dims=1)  # Average over draws
mean_pred_cats_constrained_prior = [argmax(mean_probs_constrained_prior[i, :]) for i in eachindex(reef_properties.cluster)]

n_draws, N = size(probs_constrained_prior)
# Accuracy per draw
acc_per_draw = [mean(pred_cats_constrained_prior[d, :] .== reef_properties.cluster_vals) for d in 1:n_draws]
# Summarize accuracy distribution
mean_acc = mean(acc_per_draw)
acc_ci = quantile(acc_per_draw, [0.025, 0.975])

df = DataFrame(
    bioregion = reef_properties.bioregion,
    y = reef_properties.cluster_vals,
    yhat = mean_pred_cats_weak_prior,
    reef_size = reef_properties.log_area,
    reef_depth = reef_properties.depth_med,
    dhw = reef_properties.mean_dhw,
    conn_str = reef_properties.log_total_strength
)

perf_by_bio = combine(groupby(df, :bioregion)) do sdf
    (; n = nrow(sdf),
       acc = mean(sdf.y .== sdf.yhat),
       size = mean(sdf.reef_size),
       depth = mean(sdf.reef_depth),
       dhw = mean(sdf.dhw),
       conn_str = mean(sdf.conn_str)
    )
end

first(sort!(perf_by_bio, :n, rev=true), 10)  # inspect small groups
show(perf_by_bio)

perf_by_bio[!, :max_int] .= vec(maximum(u1s; dims=1))
perf_by_bio[!, :min_int] .= vec(minimum(u1s; dims=1))
perf_by_bio[!, :int_range] .= perf_by_bio.max_int .- perf_by_bio.min_int

bioregion_colors = distinguishable_colors(nrow(perf_by_bio))

fig = Figure(size = (900, 900), fontsize=10)
ax = Axis(
    fig[1,1],
    ylabel = "Mean bioregion accuracy",
    xlabel = "Mean bioregion depth"
)
Makie.scatter!(ax, perf_by_bio.depth, perf_by_bio.acc, color=bioregion_colors)
ax = Axis(
    fig[1,2],
    ylabel = "Mean bioregion accuracy",
    xlabel = "Mean bioregion reef size"
)
Makie.scatter!(ax, perf_by_bio.size, perf_by_bio.acc, color=bioregion_colors)
ax = Axis(
    fig[2,1],
    ylabel = "Mean bioregion accuracy",
    xlabel = "Mean bioregion DHW"
)
Makie.scatter!(ax, perf_by_bio.dhw, perf_by_bio.acc, color=bioregion_colors)
ax = Axis(
    fig[2,2],
    ylabel = "Mean bioregion accuracy",
    xlabel = "Mean bioregion connectivity strength"
)
Makie.scatter!(ax, perf_by_bio.n, perf_by_bio.acc, color=bioregion_colors)
# ax = Axis(
#     fig[3,1],
#     ylabel = "Mean bioregion accuracy",
#     xlabel = "Bioregion intercept value range"
# )
# Makie.scatter!(ax, perf_by_bio.int_range, perf_by_bio.acc, color=bioregion_colors)
ax = Axis(
    fig[3,2],
    ylabel = "Mean bioregion accuracy",
    xlabel = "Bioregion number of reefs"
)
Makie.scatter!(ax, perf_by_bio.n, perf_by_bio.acc, color=bioregion_colors)

y = reef_properties.cluster_vals
log_scores = zeros(N)
for i in 1:N
    log_scores[i] = mean(-log.(probs_weak_prior[:, i, y[i]]))
end
mean_log_score = mean(log_scores)


function loo_predictor_contribution(X, predictors, y, grp1, grp2, n_grp1, n_grp2, pos_idx, neg_idx)
    # N, p = size(X)
    # n_categories = size(preds_full, 2)
    # contribution = Dict{Int, Float64}()

    for (j, predictor) in enumerate(predictors)
        X_loo = copy(X)
        X_loo[:, j] .= 0.0  # leave-one-out: zero out predictor j

        # recompute linear predictor for all posterior samples
        # here you may need to plug in your model's function to get predicted probabilities
        # Example: preds_loo = model_predict(X_loo, ...)

        # for illustration, let's assume you have a function returning p_hat given X:
        chain_loo = sample(
            ordinal_mixed_informed_priors(
                y, 
                X_loo, 
                grp1, 
                grp2, 
                n_grp1, 
                n_grp2,
                pos_idx,
                neg_idx
            ),
            NUTS(), 
            500;
            warmup=250, 
            chains=1
        )

        Serialization.serialize(joinpath(output_path, "cluster_regression_outputs/ordinal_chain_loo_less_$(predictor).jls"), chain_loo)

        # preds_loo = predict_probs(chain_loo, X_loo, bioregion_idx, GCM_idx, n_bioregions, n_GCMs)  # your existing prediction function

        # # compute log score change
        # logscore_full = -log.(preds_full[CartesianIndex(i, y[i])] for i in 1:N)
        # logscore_loo  = -log.(preds_loo[CartesianIndex(i, y[i])] for i in 1:N)

        # contribution[predictor] = mean(logscore_loo) - mean(logscore_full)
    end

    return nothing
end

loo_predictor_contribution(
    X, 
    predictors,
    clusters, 
    bioregion_idx, 
    GCM_idx, 
    n_bioregions, 
    n_GCMs,
    findfirst(predictors .== :depth_med),
    findfirst(predictors .== :mean_dhw)
)