"""
Load initial result set from GBR-wide ReefMod domain.
"""

using GLMakie, GeoMakie, GraphMakie

using Revise, Infiltrator

using ADRIA

#moore_domain_path = "c:/Users/bgrier/Documents/Projects/ADRIA_Domains/Moore_2024-02-14_v060_rc1/"
# Local scale domain
#moore_dom = ADRIA.load_domain(moore_domain_path)

include("common.jl")

# GBR wide domain
gbr_dom = ADRIA.load_domain(RMEDomain, gbr_domain_path, "45")

# gbr_dom.conn
# for reef_id in context_layers[context_layers.target_reefs, :RME_GBRMPA_ID]
#     if any(collect(gbr_dom.conn.Source) .== reef_id)
#         gbr_dom.conn[Source = At(reef_id)] .= 0.0
#     end
# end

#ADRIA.fix_factor!(gbr_dom; dhw_scenario=0)

# generate 4096 sample scenarios from counterfactual scenarios
scens = ADRIA.sample_cf(gbr_dom, 4096)

# repeat(scens, 16)

# scens_repeat = ADRIA.param_table(gbr_dom)
# scens_repeat = repeat(scens_repeat, 5)
# scens_repeat[!, :dhw_scenario] .= [1, 2, 3, 4, 5]
# rs = ADRIA.run_scenarios(gbr_dom, scens_repeat, "45")

# Run sampled scenarios for a given RCP
rs = ADRIA.run_scenarios(gbr_dom, scens, "45")

# ... or repeat scenario runs across multiple RCPs
#rs = ADRIA.run_scenarios(dom, scens, ["45", "60", "85"])

# ADRIA.viz.dhw_scenario(gbr_dom, 1)
# test_dhw = gbr_dom.dhw_scens
# test_dhw = Float64.(mapslices(median, test_dhw, dims=[:scenarios]))
# test_dhw = median(test_dhw; dims=:scenarios)
