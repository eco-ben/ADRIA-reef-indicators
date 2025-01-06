"""
Perfom lagged correlation analysis across reefs within bioregions. There are ~30 bioregions
across the reefs in the GBR. This analysis is slightly finer spatial scale than Subregion.
If 3_Subregion.jl has been run immediately before then the result reefs can be compared.
"""

using GLMakie, GeoMakie, GraphMakie

using YAXArrays, DimensionalData

using DataFrames, Statistics, YAXArrays

using ADRIA, CoralBlox

import GeoDataFrames as GDF

include("common.jl")
include("4_Bioregion_correlation.jl")

results_allgbr = DataFrame(gbr_reef=[], target_reef=[], lag3=[], lag4=[], lag5=[], lag6=[], lag7=[], lag8=[])
for t_reef in target_reefs
    t_reef_cover = rel_cover[:, rel_cover.sites .== [t_reef]]
    t_reef_id = t_reef_cover.sites[1]
    t_reef_cover = t_reef_cover[:, 1]
    gbr_cover = rel_cover[:, rel_cover.sites .!= [t_reef]]

    for (ind, x) in enumerate(eachcol(gbr_cover))
        corr = cross_correlation(t_reef_cover, x, 3:8)
        gbr_name = gbr_cover.sites[ind]

        push!(results_allgbr, [gbr_name; t_reef_id; corr])
    end
end

hist(results.lag8)


ports_s_tsv = ["Hay Point", "Port of Rockhampton", "Port of Gladstone", "Bundaberg", "Maryborough"]
reefs_s_tsv  = context_layers[(context_layers.closest_port .∈ [ports_s_tsv]), :UNIQUE_ID]
results_s_tsv = DataFrame(closest_port=[], target_reef=[], lag3=[], lag4=[], lag5=[], lag6=[], lag7=[], lag8=[], lag9=[], lag10=[])

for t_reef in target_reefs
    t_reef_cover = rel_cover[:, rel_cover.sites .== [t_reef]]
    t_reef_id = t_reef_cover.sites[1]
    t_reef_cover = t_reef_cover[:, 1]
    #reefs_s_cover = rel_cover[:, rel_cover.sites .∈ [reefs_s_tsv]]

    for x in ports_s_tsv
        port_reefs = context_layers[context_layers.closest_port .== [x], :UNIQUE_ID]
        port_cover = rel_cover[:, rel_cover.sites .∈ [port_reefs]]
        port_median = Float64.(mapslices(median, port_cover, :sites))

        corr = cross_correlation(t_reef_cover, port_median, 3:10)

        push!(results_s_tsv, [x; t_reef_id; corr])
    end
end

hist(results_s_tsv.lag8)
