using Revise, Infiltrator
using YAXArrays

change_ADRIA_debug(false) # Ensure ADRIA debug mode is set to false to allow parallel processing.

using ADRIA

GBR_domain_path = "../../ADRIA Domains/GBR_2024_10_15_HighResCoralStress/"
dhw_scenarios = open_dataset("../../ADRIA Domains/GBR_2024_10_15_HighResCoralStress/DHWs/dhwRCP45.nc")
gcms = dhw_scenarios.dhw.properties["members"]

# GBR wide domain
gbr_dom = ADRIA.load_domain(GBR_domain_path, "45")
context_layers = gbr_dom.loc_data

for (g, GCM) in enumerate(gcms)

    # generate 4096 sample scenarios from counterfactual scenarios
    scens = ADRIA.sample_cf(gbr_dom, 4096)

    # change dhw scenarios to GFDL-CM4
    scens[!, :dhw_scenario] .= g

    # Run sampled scenarios for a given RCP
    rs = ADRIA.run_scenarios(gbr_dom, scens, "45")

    # Rename result store directory according to GCM that is chosen
end
