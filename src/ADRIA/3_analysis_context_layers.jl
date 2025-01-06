"""
Access the domain connectivity and calculate the strength and number of outgoing/incoming
connections for each reef. Run after running 3_.jl to include the bioregions and identified
reef layers.
"""

using DataFrames, Statistics, YAXArrays

using GLMakie, GeoMakie, GraphMakie

using ADRIA

import GeoDataFrames as GDF

import GeoFormatTypes as GFT

include("common.jl")

# Load context layers with bioregions and target reefs
context_layers = GDF.read("data/context_layers_targetted.gpkg")
context_layers.RME_GBRMPA_ID = ifelse.((context_layers.RME_GBRMPA_ID .== "20198"), "20-198", context_layers.RME_GBRMPA_ID)

# GBR wide domain
gbr_domain_path = "c:/Users/bgrier/Documents/Projects/ADRIA_Domains/rme_ml_2024_01_08/"
gbr_dom = ADRIA.load_domain(RMEDomain, gbr_domain_path, "45")

# 1. Attach connectivity data to context_layers
connectivity_matrix = gbr_dom.conn

# Calculate connectivity metrics for context analysis
reefs = collect(getAxis("Source", connectivity_matrix).val)

col_types = [
    String,
    Float64,
    Int64,
    Float64,
    Float64,
    Int64,
    Float64,
    Float64,
    Int64,
    Float64,
    Float64
]

source_to_sink = DataFrame(
    [T[] for T in col_types],
    [
    :RME_GBRMPA_ID,
    :income_strength,
    :income_count,
    :income_comb,
    :out_strength,
    :out_count,
    :out_comb,
    :total_strength,
    :total_count,
    :total_comb,
    :so_to_si
    ]
)

for reef in eachindex(reefs)
    reef_id = reefs[reef]
    outgoing = connectivity_matrix[connectivity_matrix.Source .== reef_id, :]
    incoming = connectivity_matrix[:, connectivity_matrix.Sink .== reef_id]

    income_strength = sum(incoming)
    income_count = count(incoming .> 0)
    income_comb = income_strength / income_count

    out_strength = sum(outgoing)
    out_count = count(outgoing .> 0)
    out_comb = out_strength / out_count

    total_strength = income_strength + out_strength
    total_count = income_count + out_count
    total_comb = total_strength / total_count

    so_to_si = out_comb / income_comb

    push!(
        source_to_sink,
        [
        reef_id,
        income_strength,
        income_count,
        income_comb,
        out_strength,
        out_count,
        out_comb,
        total_strength,
        total_count,
        total_comb,
        so_to_si
        ]
    )
end
context_layers = leftjoin(context_layers, source_to_sink; on=:RME_GBRMPA_ID, order=:left)

# 2. Attach initial coral cover data to context_layers
initial_cc = gbr_dom.init_coral_cover
total_icc = Float64.(mapslices(sum, initial_cc, dims=[:species]))
context_layers.initial_coral_cover = total_icc.data

# species =
# [
# "tabular Acropora_2_1", "tabular Acropora_2_2",
# "tabular Acropora_2_3", "tabular Acropora_2_4", "tabular Acropora_2_5", "tabular Acropora_2_6", "tabular Acropora_2_7",
# "corymbose Acropora_3_1", "corymbose Acropora_3_2", "corymbose Acropora_3_3", "corymbose Acropora_3_4",
# "corymbose Acropora_3_5", "corymbose Acropora_3_6", "corymbose Acropora_3_7", "corymbose non-Acropora_4_1", "corymbose non-Acropora_4_2",
# "corymbose non-Acropora_4_3", "corymbose non-Acropora_4_4", "corymbose non-Acropora_4_5", "corymbose non-Acropora_4_6", "corymbose non-Acropora_4_7",
# "Small massives_5_1", "Small massives_5_2", "Small massives_5_3", "Small massives_5_4", "Small massives_5_5",
# "Small massives_5_6", "Small massives_5_7", "Large massives_6_1", "Large massives_6_2", "Large massives_6_3",
# "Large massives_6_4", "Large massives_6_5", "Large massives_6_6", "Large massives_6_7"
# ]
# spec_groups = getindex.(species, [1:15])

# icc = DataFrame(Matrix(initial_cc.data), collect(getAxis("locs", initial_cc).val))
# icc.spec_groups = spec_groups
# gdf = DataFrames.groupby(icc, :spec_groups)
# mapcols(gdf, sum)


# 3. Attach DHW data to context_layers
# mean DHW data for sites
dhw_time_scens = gbr_dom.dhw_scens
dhw_time = Float64.(mapslices(median, dhw_time_scens, dims=[:scenarios]))
dhw_locs = Float64.(mapslices(mean, dhw_time, dims=[:timesteps]))
dhw_locs = DataFrame(mean_dhw = dhw_locs.data, RME_GBRMPA_ID = collect(getAxis("locs", dhw_locs).val))

context_layers = leftjoin(context_layers, dhw_locs, on=:RME_GBRMPA_ID, order=:left)

# DHW correlation for sites
rs = ADRIA.load_results("outputs/ADRIA-out/ReefMod Engine__RCPs_45__2024-06-19_10_58_55_647")
#rs = ADRIA.load_results("outputs/ADRIA-out/ReefMod Engine__RCPs_45__2024-07-26_12_05_07_844")
#rs = ADRIA.load_results("outputs/ADRIA-out/ReefMod Engine__RCPs_85__2024-07-12_09_49_08_988")
tac = ADRIA.metrics.total_absolute_cover(rs)
tac_sites = Float64.(mapslices(median, tac, dims=[:scenarios]))
tac_sites_reduced = tac_sites[timesteps=2:79]
rel_cover = Float64.(mapslices(relative_site_cover, tac_sites_reduced, dims=[:timesteps]))

#rel_dhw = Float64.(mapslices(relative_site_cover, dhw_time, dims=[:timesteps]))
dhw_cor = DataFrame(UNIQUE_ID = [], dhw_cor = [])
for (ind, reef) in enumerate(eachcol(rel_cover))

    reef_cover = rel_cover[:,ind]
    id = collect(getAxis("sites", rel_cover).val)[ind]
    reef_dhw = dhw_time[:,ind]
    reef_dhw.data .= -1 * reef_dhw.data

    correlation = cross_correlation(reef_cover, reef_dhw, [0], true)
    push!(dhw_cor, [id; correlation])

end
context_layers = leftjoin(context_layers, dhw_cor, on=:UNIQUE_ID, order=:left)

# 4. Attach site data/area to context_layers
gbr_site_data = gbr_dom.site_data[:,[:UNIQUE_ID, :area]]
context_layers = leftjoin(context_layers, gbr_site_data; on=:UNIQUE_ID, order=:left)

# 5. Calculate the number of reefs in each bioregion
context_layers.bioregion_count .= 1
for reef in eachrow(context_layers)
    reef.bioregion_count = count(context_layers.bioregion .== [reef.bioregion])
end

# 6. Calculate eigenvector_centrality scores for each reef - GBR-wide connections weighted by subregions
mean_conn, std_conn = load_connectivity(conn_path, context_layers.RME_UNIQUE_ID)

mean_conn[diagind(mean_conn)] .= 0.0

conn_whole_gbr = connectivity_scoring(mean_conn)
conn_mgmt_area = connectivity_scoring(mean_conn; gdf=context_layers, context_layer=:management_area, conn_col_name=:conn_score_mgmt, by_layer=true)
conn_subr_area = connectivity_scoring(mean_conn; gdf=context_layers, context_layer=:closest_port, conn_col_name=:conn_score_subr, by_layer=true)
conn_bior_area = connectivity_scoring(mean_conn; gdf=context_layers, context_layer=:bioregion, conn_col_name=:conn_score_bior, by_layer=true)

context_layers = leftjoin(context_layers, conn_whole_gbr, on=:RME_UNIQUE_ID, order=:left)
context_layers = leftjoin(context_layers, conn_mgmt_area, on=:RME_UNIQUE_ID, order=:left)
context_layers = leftjoin(context_layers, conn_subr_area, on=:RME_UNIQUE_ID, order=:left)
context_layers = leftjoin(context_layers, conn_bior_area, on=:RME_UNIQUE_ID, order=:left)

context_layers = weight_by_context(context_layers, :conn_score_mgmt, :management_area, :conn_ranked_mgmt)
context_layers = weight_by_context(context_layers, :conn_score_subr, :closest_port, :conn_ranked_subr)
context_layers = weight_by_context(context_layers, :conn_score_bior, :bioregion, :conn_ranked_bior)

# Format columns for writing to geopackage
context_layers.income_strength = Float64.(context_layers.income_strength)
context_layers.income_count .= convert.(Int64, context_layers.income_count)
context_layers.income_comb .= convert.(Float64, context_layers.income_comb)
context_layers.out_strength .= convert.(Float64, context_layers.out_strength)
context_layers.out_count .= convert.(Int64, context_layers.out_count)
context_layers.out_comb .= convert.(Float64, context_layers.out_comb)
context_layers.total_strength .= convert.(Float64, context_layers.total_strength)
context_layers.total_count .= convert.(Int64, context_layers.total_count)
context_layers.total_comb .= convert.(Float64, context_layers.total_comb)
context_layers.so_to_si .= convert.(Float64, context_layers.so_to_si)
context_layers.mean_dhw .= convert.(Float64, context_layers.mean_dhw)
context_layers.dhw_cor = ifelse.(isnan.(context_layers.dhw_cor), 0.0, context_layers.dhw_cor)
context_layers.area .= convert.(Float64, context_layers.area)
context_layers.conn_score .= convert.(Float64, context_layers.conn_score)
context_layers.conn_score_mgmt .= convert.(Float64, context_layers.conn_score_mgmt)
context_layers.conn_score_subr .= convert.(Float64, context_layers.conn_score_subr)
context_layers.conn_score_bior .= convert.(Float64, context_layers.conn_score_bior)

GDF.write("data/analysis_context_layers.gpkg", context_layers; crs=GFT.EPSG(7844))
