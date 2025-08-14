
using CSV
using CategoricalArrays
using Econometrics, StatsModels

include("../../common.jl")

context_layers = GDF.read(joinpath(output_path, "analysis_context_layers_carbonate.gpkg"))
context_layers.log_so_to_si = log10.(context_layers.so_to_si)
context_layers.log_total_strength = log10.(context_layers.total_strength)
context_layers.log_area = log10.(context_layers.area)

dhw_scenarios = open_dataset(joinpath(gbr_domain_path, "DHWs/dhwRCP45.nc"))
GCMs = dhw_scenarios.dhw.properties["members"]

# Prepare long-form of data to support multinomial analysis.
# Each row should be unique attribute combination of reef properties:
# reef, depth, connectivity, dhw, cluster and GCM

gcm_dhw_cols = [Symbol("$(GCM)_mean_dhw") for GCM in GCMs]
gcm_bioregion_cluster_cols = [Symbol("$(GCM)_bioregion_cluster_cats") for GCM in GCMs]

collated_reef_properties = Vector{DataFrame}(undef, length(GCMs))
target_cols = [:UNIQUE_ID, :GBRMPA_ID, :depth_med, :log_total_strength, :log_so_to_si, :log_area, :bioregion]
reef_properties = context_layers[:, target_cols]
for (i, GCM) in enumerate(GCMs)
    dhw_col = Symbol("$(GCM)_mean_dhw")
    bio_cluster_col = Symbol("$(GCM)_bioregion_cluster_cats")

    bio_cluster_details = context_layers[:, [dhw_col, bio_cluster_col]]
    rename!(bio_cluster_details, dhw_col => :mean_dhw, bio_cluster_col => :cluster)
    bio_cluster_details.GCM .= GCM
    bio_cluster_details.cluster .= categorical(bio_cluster_details.cluster, levels=["low", "medium", "high"])

    collated_reef_properties[i] = hcat(
        copy(reef_properties),
        bio_cluster_details
    )
end

reef_properties = vcat(collated_reef_properties...)
reef_properties.GCM .= categorical(reef_properties.GCM)

lhs = Term(:cluster)
rhs_terms = [
    Term(:depth_med),
    Term(:log_total_strength),
    Term(:log_so_to_si),
    Term(:mean_dhw),
    Term(:log_area),
    Term(:GCM)
]

f = FormulaTerm(lhs, sum(rhs_terms))

model = fit(
    EconometricModel,
    f,
    reef_properties,
    contrasts=Dict(:cluster => DummyCoding(base="low"))
)

coeftab = DataFrame(coeftable(model))

cluster_predictions = [argmax(predict(model)[i, :]) for i in 1:nrow(reef_properties)]
mean(cluster_predictions .== reef_properties.cluster.refs)

regression_save_path = joinpath(
    output_path,
    "cluster_regression_outputs/bioregion_cluster_regression_coefficients.csv"
)
CSV.write(regression_save_path, coeftab)

df = DataFrame(
    bioregion = reef_properties.bioregion,
    y = reef_properties.cluster.refs,
    pred_y = cluster_predictions,
    reef_size = reef_properties.log_area,
    reef_depth = reef_properties.depth_med,
    dhw = reef_properties.mean_dhw,
    conn_str = reef_properties.log_total_strength
)

perf_by_bio = combine(groupby(df, :bioregion)) do sdf
    (; n = nrow(sdf),
       acc = mean(sdf.pred_y .== sdf.y),
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

function loo_predictor_contribution(base_model_predictors, reef_properties, base_model_predictions)
    y = reef_properties.cluster.refs
    # N, p = size(X)
    # n_categories = size(preds_full, 2)
    # contribution = Dict{Int, Float64}()
    N = nrow(reef_properties)
    logscore_full = -log.(base_model_predictions[CartesianIndex(i, y[i])] for i in 1:N)

    contributions = Dict(zip(base_model_predictors, zeros(length(base_model_predictors))))

    for (j, predictor) in enumerate(base_model_predictors)            
        lhs = Term(:cluster)
        rhs_terms = [
            Term(:depth_med),
            Term(:log_total_strength),
            Term(:log_so_to_si),
            Term(:mean_dhw),
            Term(:log_area),
            Term(:GCM)
        ]
        rhs_terms = rhs_terms[rhs_terms .!= Term(predictor)]

        f = FormulaTerm(lhs, sum(rhs_terms))

        model = fit(
            EconometricModel,
            f,
            reef_properties,
            contrasts=Dict(:cluster => DummyCoding(base="low"))
        )

        coeftab = DataFrame(coeftable(model))

        cluster_predictions = predict(model)
        # mean(cluster_predictions .== reef_properties.cluster.refs)
            
            # # compute log score change

        logscore_loo  = -log.(cluster_predictions[CartesianIndex(i, y[i])] for i in 1:N)

        contributions[predictor] = mean(logscore_loo) - mean(logscore_full)
    end

    return nothing
end