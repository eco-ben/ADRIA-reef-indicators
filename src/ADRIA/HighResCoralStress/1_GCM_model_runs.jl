using YAXArrays, NetCDF
using ArchGDAL
using CSV

include("../../common.jl")

# Ensure ADRIA debug mode is set to false to allow parallel processing.
change_ADRIA_debug(false)

using ADRIA

dhw_scenarios = open_dataset(joinpath(gbr_domain_path, "DHWs/dhwRCP45.nc"))
gcms = dhw_scenarios.dhw.properties["members"]

# GBR wide domain
gbr_dom = ADRIA.load_domain(gbr_domain_path, "45")
context_layers = gbr_dom.loc_data

# Need to recast gbr_dom geometry col for ADRIA.run_scenarios().
# Possibly issue with GeoDataFrames version.
gbr_dom.loc_data.geometry = Vector{ArchGDAL.IGeometry}(gbr_dom.loc_data.geometry)

for (g, GCM) in enumerate(gcms)

    # generate 4096 sample scenarios from counterfactual scenarios
    scens = ADRIA.sample_cf(gbr_dom, 4096)

    # change dhw scenarios to selected GCM
    scens[!, :dhw_scenario] .= g

    CSV.write(joinpath(output_path, "HighResCoralStress_ADRIA_scens_$(GCM).csv"), scens)

    # Run sampled scenarios for a given RCP
    rs = ADRIA.run_scenarios(gbr_dom, scens, "45")

    # Rename result store directory according to GCM that is chosen
end

# Save output absolute median cover to allow for data sharing.
for GCM in gcms
    processed_outputs = joinpath(output_path, "processed_model_outputs/median_cover_$(GCM).nc")
    if !isfile(processed_outputs)

        rs = ADRIA.load_results(joinpath(output_path, "GBR_2024_10_15_HighResCoralStress__RCPs_45_$(GCM)"))
        scenario_cover = ADRIA.metrics.total_absolute_cover(rs)
        median_cover = ADRIA.metrics.loc_trajectory(median, scenario_cover)

        dims = (
            Dim{:timesteps}(Int64.(scenario_cover.timesteps)),
            Dim{:locations}(string.(scenario_cover.locations))
        )
        properties = median_cover.properties
        pop!(properties, :is_relative)
        properties = Dict{String,Any}([(string.(k), v) for (k, v) in properties])

        median_cover = rebuild(median_cover, dims=dims, metadata=properties)
        savecube(median_cover, processed_outputs, driver=:netcdf, overwrite=true)
    end
end
