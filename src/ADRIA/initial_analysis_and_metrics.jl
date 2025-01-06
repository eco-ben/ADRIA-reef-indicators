"""
Load the initially created result set and calculate the metrics needed for later analyses.
(Currently these metrics are recalculated in later scripts.)
"""

using GLMakie, GeoMakie, GraphMakie

using YAXArrays, DimensionalData

using DataFrames, Statistics, CSV

using ADRIA, CoralBlox

import GeoDataFrames as GDF

include("common.jl")

#rs = ADRIA.load_results("outputs/ADRIA-out/ReefMod Engine__RCPs_45__2024-06-19_10_58_55_647")
rs = ADRIA.load_results("outputs/ADRIA-out/ReefMod Engine__RCPs_45__2024-07-26_12_05_07_844")

context_layers = find_latest_file("../canonical-reefs/output/")
context_layers = GDF.read(context_layers)

# Analysis
s_tac = ADRIA.metrics.scenario_total_cover(rs)
ADRIA.viz.scenarios(rs, s_tac)

tac = ADRIA.metrics.total_absolute_cover(rs)
#tac_sites = ADRIA.metrics.loc_trajectory(median, tac)

tac_sites = Float64.(mapslices(mean, tac, dims=[:scenarios]))

df = DataFrame(tac_sites.data, collect(getAxis("sites", tac_sites).val))
#CSV.write("outputs/tac_site_series_19_06_2024_depth_new.csv", df)
#df = CSV.read("outputs/tac_site_series_19_06_2024_depth_new.csv", DataFrame)

rel_cover = mapcols(relative_site_cover, df[2:79,:])
rel_cover.year = [string(i+1) for i in 1:size(rel_cover,1)]
select!(rel_cover, :year, Not(:year))

data = permutedims(rel_cover, 1, "UNIQUE_ID")

data = leftjoin(data, context_layers[:, [:UNIQUE_ID, :target_reefs]], on=:UNIQUE_ID)
filtered = data[(data.UNIQUE_ID .∈ [filtered_bior.UNIQUE_ID]),:]

data = leftjoin(data, context_layers[:, [:UNIQUE_ID, :management_area]], on=:UNIQUE_ID)

data.target_reefs = ifelse.(data.target_reefs, "target", "non_target")
f, ax = plot_lines(filtered, :target_reefs, 2, 15, "Year", "Proportion of reef cover Year2", 0.2)
hlines!(1, color=:green)
save("figs/all_reefs_initial_man_areas.png", f)


f, ax = plot_lines(data, :management_area, 2, 17, "Year", "Proportion of reef cover Year2", 0.1)
hlines!(1, color=:green)
hlines!(0.1, color=:blue)
save("figs/initial_man_areas_15_years.png", f)

no_conn_crash = collect(getAxis("sites", rel_cover).val)[rel_cover.data[2,:] .> 0.5]
no_conn_crash = findall(collect(getAxis("sites", rel_cover).val) .∈ [no_conn_crash])
no_conn = rel_cover[sites = no_conn_crash]
no_conn = Float64.(mapslices(median, no_conn, dims=[:sites]))

w_conn_crash = collect(getAxis("sites", rel_cover).val)[rel_cover.data[2,:] .> 0.5]
w_conn_crash = findall(collect(getAxis("sites", rel_cover).val) .∈ [w_conn_crash])
w_conn = rel_cover[sites = w_conn_crash]
w_conn = Float64.(mapslices(median, w_conn, dims=[:sites]))

diff = w_conn - no_conn

series(no_conn.data', solid_color=:blue)
series!(w_conn.data', solid_color=:red)
series!(diff.data', solid_color=:black)


tac_rsa = ADRIA.sensitivity.rsa(rs, rel_cover; S=10)
