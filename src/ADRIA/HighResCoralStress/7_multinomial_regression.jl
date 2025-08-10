
using CSV
using CategoricalArrays
using Econometrics, StatsModels
using MixedModels

include("../../common.jl")

context_layers = GDF.read(joinpath(output_path, "analysis_context_layers_carbonate.gpkg"))
context_layers.log_so_to_si = log10.(context_layers.so_to_si)
context_layers.log_total_strength = log10.(context_layers.total_strength)

dhw_scenarios = open_dataset(joinpath(gbr_domain_path, "DHWs/dhwRCP45.nc"))
GCMs = dhw_scenarios.dhw.properties["members"]

# Prepare long-form of data to support multinomial analysis.
# Each row should be unique attribute combination of reef properties:
# reef, depth, connectivity, dhw, cluster and GCM

gcm_dhw_cols = [Symbol("$(GCM)_mean_dhw") for GCM in GCMs]
gcm_bioregion_cluster_cols = [Symbol("$(GCM)_bioregion_clusters") for GCM in GCMs]

collated_reef_properties = Vector{DataFrame}(undef, length(GCMs))
target_cols = [:UNIQUE_ID, :GBRMPA_ID, :depth_med, :log_total_strength, :log_so_to_si, :bioregion]
reef_properties = context_layers[:, target_cols]
for (i, GCM) in enumerate(GCMs)
    dhw_col = Symbol("$(GCM)_mean_dhw")
    bio_cluster_col = Symbol("$(GCM)_bioregion_clusters")

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

using Turing, Distributions, StatsFuns, LinearAlgebra

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
@model function ordinal_random_intercepts(y, X, grp1, grp2, n_grp1, n_grp2)
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
    for i in 1:N
        η = dot(X[i, :], β) + u1[grp1[i]] + u2[grp2[i]]
        p1 = cdf(Normal(), c[1] - η)
        p2 = cdf(Normal(), c[2] - η)
        ps = [p1, p2 - p1, 1 - p2]  # probabilities for categories 1,2,3
        ps = clamp.(ps, 1e-10, 1.0) # Prevent errors if probabilities appear as 0 or -ve
        ps ./= sum(ps)

        y[i] ~ Distributions.Categorical(ps)
    end
end

reef_properties.GCM = categorical(reef_properties.GCM)
reef_properties.bioregion = categorical(reef_properties.bioregion)

# Encode random effects groups (GCM and bioregion) as int
GCM_idx = reef_properties.GCM.refs
bioregion_idx = reef_properties.bioregion.refs

n_GCMs = length(levels(reef_properties.GCM))
n_bioregions = length(levels(reef_properties.bioregion))

# Encode response levels as increasing integers
y = reef_properties.cluster
K = length(unique(reef_properties.cluster))

# Design matrix (without intercept — model will handle that)
predictors = [:depth_med, :log_total_strength, :log_so_to_si, :mean_dhw]
predictors = [:depth_med, :log_total_strength]
X = Matrix(DataFrames.select(reef_properties, predictors))

model = ordinal_random_intercepts(
    y, 
    X, 
    bioregion_idx, 
    GCM_idx,
    n_bioregions, 
    n_GCMs
)
chain = sample(model, NUTS(), MCMCThreads(), 2000, 4)  # 4 chains using Threads, 2000 samples each

small_chain = sample(model, NUTS(), 150; warmup=100, chains=1)