
using CSV
using CategoricalArrays
using StatsModels
using Turing, Distributions, StatsFuns, LinearAlgebra
using MLDataUtils: shuffleobs, stratifiedobs, rescale!

include("../../common.jl")

context_layers = GDF.read(joinpath(output_path, "analysis_context_layers_carbonate.gpkg"))
context_layers.log_so_to_si = log10.(context_layers.so_to_si)
context_layers.log_total_strength = log10.(context_layers.total_strength)
context_layers.log_area = log10.(context_layers.area)

dhw_scenarios = open_dataset(joinpath(gbr_domain_path, "DHWs/dhwRCP45.nc"))
GCMs = dhw_scenarios.dhw.properties["members"]
GCMs = ["EC-Earth3-Veg", "ACCESS-ESM1-5", "ACCESS-CM2", "NorESM2-MM", "GFDL-CM4"]

# Prepare long-form of data to support multinomial analysis.
# Each row should be unique attribute combination of reef properties:
# reef, depth, connectivity, dhw, cluster and GCM

gcm_dhw_cols = [Symbol("$(GCM)_mean_dhw") for GCM in GCMs]
gcm_bioregion_cluster_cols = [Symbol("$(GCM)_bioregion_cluster_cats") for GCM in GCMs]

collated_reef_properties = Vector{DataFrame}(undef, length(GCMs))
target_cols = [:UNIQUE_ID, :GBRMPA_ID, :depth_med, :log_total_strength, :log_so_to_si, :bioregion, :log_area]
reef_properties = context_layers[:, target_cols]
for (i, GCM) in enumerate(GCMs)
    dhw_col = Symbol("$(GCM)_mean_dhw")
    bio_cluster_col = Symbol("$(GCM)_bioregion_cluster_cats")

    bio_cluster_details = context_layers[:, [dhw_col, bio_cluster_col]]
    rename!(bio_cluster_details, dhw_col => :mean_dhw, bio_cluster_col => :cluster)
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

"""
    ordinal_random_intercepts(y, X, grp1, grp2, n_grp1, n_grp2)

# Arguments
- `y` : Integer vector of reef cluster assignment levels.
- `X` : n-samples * n-predictors Matrix of independent variables
- `grp1` : Integer vector of the first categorical random effects term
- `grp2` : Integer vector of the second categorical random effects term
- `n_grp1` : Number of unique levels in grp1
- `n_grp2` : Number of unique levels in grp2
"""
@model function ordinal_random_intercepts(y::Vector{UInt32}, X::Matrix{Float64}, grp1::Vector{UInt32}, grp2::Vector{UInt32}, n_grp1::Int64, n_grp2::Int64)
    N, p = size(X)
    # @assert p == 4

    # priors for fixed effects
    # weakly-informative prior on beta
    # β ~ MvNormal(zeros(p), 5.0^2 * I)
    β ~ filldist(Normal(0, 5), p)

    # --- RANDOM EFFECTS: non-centered parameterization ---
    # grp1
    σ1 ~ truncated(Cauchy(0, 2), 0, Inf)
    z1 ~ filldist(Normal(0,1), n_grp1)
    # grp2
    σ2 ~ truncated(Cauchy(0, 2), 0, Inf)
    z2 ~ filldist(Normal(0,1), n_grp2)

    # transform to actual random intercepts
    u1 = σ1 .* z1
    u2 = σ2 .* z2

    # estimates for cutpoints - required for ordinal bayesian model to separate the 3 levels of response variable
    c_raw ~ filldist(Normal(0, 5), 2)
    c = sort(c_raw)  # c[1] < c[2]


    # Likelihood
    # for i in 1:N
    #     η = dot(X[i, :], β) + u1[grp1[i]] + u2[grp2[i]]
    #     p1 = cdf(Normal(), c[1] - η)
    #     p2 = cdf(Normal(), c[2] - η)
    #     ps = [p1, p2 - p1, 1 - p2]  # probabilities for categories 1,2,3
    #     ps = clamp.(ps, 1e-10, 1.0) # Prevent errors if probabilities appear as 0 or -ve
    #     ps ./= sum(ps)

    #     y[i] ~ Distributions.Categorical(ps)
    # end

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

reef_properties.GCM = categorical(reef_properties.GCM)
reef_properties.bioregion = categorical(reef_properties.bioregion)

# Encode random effects groups (GCM and bioregion) as int
GCM_idx = reef_properties.GCM.refs
bioregion_idx = reef_properties.bioregion.refs

n_GCMs = length(levels(reef_properties.GCM))
n_bioregions = length(levels(reef_properties.bioregion))

# Encode response levels as increasing integers
reef_properties.cluster_vals = categorical(reef_properties.cluster; ordered=true, levels=["low", "medium", "high"]).refs
K = length(unique(reef_properties.cluster))
reef_properties.row_n = 1:nrow(reef_properties)

# Design matrix (without intercept — model will handle that)
predictors = [:depth_med, :log_total_strength, :log_so_to_si, :mean_dhw, :log_area]

function split_data(df, target; at=0.70)
    shuffled = shuffleobs(df)
    return trainset, testset = stratifiedobs(row -> row[target], shuffled; p=at)
end

target = :cluster_vals

trainset, testset = split_data(reef_properties, target; at=0.50)
for feature in predictors
    μ, σ = rescale!(trainset[!, feature]; obsdim=1)
    rescale!(testset[!, feature], μ, σ; obsdim=1)
end

train_idx = indexin(trainset.row_n, reef_properties.row_n)
test_idx = indexin(testset.row_n, reef_properties.row_n)

# Turing requires data in matrix form, not dataframe
train = Matrix(trainset[:, predictors])
test = Matrix(testset[:, predictors])
train_cluster = trainset[:, target]
test_cluster = testset[:, target];
train_GCM = GCM_idx[train_idx]
train_bioregion = bioregion_idx[train_idx]
test_GCM = GCM_idx[test_idx]
test_bioregion = bioregion_idx[test_idx]

model = ordinal_random_intercepts(
    train_cluster, 
    train, 
    train_bioregion, 
    train_GCM,
    n_bioregions, 
    n_GCMs
)
# chain = sample(model, NUTS(), MCMCThreads(), 1500, 4)  # 4 chains using Threads, 2000 samples each

# StatsPlots.plot(chain)

small_chain = sample(model, NUTS(), 150; warmup=100, chains=1)

description = describe(small_chain)

summary = leftjoin(DataFrame(summarize(small_chain)), DataFrame(quantile(small_chain)), on=:parameters)

beta_params = summary[1:5, :]
fig = Figure()
ax = Axis(fig[1,1], ylabel = "Parameter", xlabel = "Estimates mean (2.5-97.5% quantiles)",
    yticks = (1:5, String.(predictors)))
Makie.scatter!(ax, beta_params.mean, 1:5; color=:blue)
lines_x = [[beta_params[x, "2.5%"], beta_params[x, "97.5%"]] for x in 1:5]
lines_y = [[x, x] for x in 1:5]

Makie.lines!.(lines_x, lines_y)