
using CSV
using CategoricalArrays
using Econometrics, StatsModels

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
gcm_bioregion_cluster_cols = [Symbol("$(GCM)_bioregion_cluster_cats") for GCM in GCMs]

collated_reef_properties = Vector{DataFrame}(undef, length(GCMs))
target_cols = [:UNIQUE_ID, :GBRMPA_ID, :depth_med, :log_total_strength, :log_so_to_si]
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

lhs = Term(:cluster)
rhs_terms = [
    Term(:depth_med),
    Term(:log_total_strength),
    Term(:log_so_to_si),
    Term(:mean_dhw),
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

regression_save_path = joinpath(
    output_path,
    "cluster_regression_outputs/bioregion_cluster_regression_coefficients.csv"
)
CSV.write(regression_save_path, coeftab)
