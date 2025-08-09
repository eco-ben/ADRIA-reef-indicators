using Econometrics, CategoricalArrays, CSV

include("../../common.jl")

context_layers = GDF.read(joinpath(output_path, "analysis_context_layers_carbonate.gpkg"))
context_layers.log_so_to_si = log10.(context_layers.so_to_si)
context_layers.log_total_strength = log10.(context_layers.total_strength)


dhw_scenarios = open_dataset(joinpath(gbr_domain_path, "DHWs/dhwRCP45.nc"))
GCMs = dhw_scenarios.dhw.properties["members"]

for GCM in GCMs
    context_layers[!, "$(GCM)_bioregion_cluster_cats"] = categorical(context_layers[:, "$(GCM)_bioregion_cluster_cats"])

    lhs = Term(Symbol("$(GCM)_bioregion_cluster_cats"))
    rhs_terms = [
        Term(:depth_med),
        Term(:log_total_strength),
        Term(:log_so_to_si),
        Term(Symbol("$(GCM)_mean_dhw"))
    ]

    f = FormulaTerm(lhs, sum(rhs_terms))

    model = fit(
        EconometricModel,
        f,
        context_layers,
        contrasts = Dict(Symbol("$(GCM)_bioregion_cluster_cats") => DummyCoding(base = "low"))
    )

    coeftab = DataFrame(coeftable(model))
    CSV.write(joinpath(output_path, "cluster_regression_outputs/$(GCM)_bioregion_cluster_regression_coefficients.csv"), coeftab)
end
