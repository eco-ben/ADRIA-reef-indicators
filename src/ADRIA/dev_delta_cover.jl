"""
Investigate the series of instantaneous change in coral cover (delta cover) between timesteps.
"""
using YAXArrays, DimensionalData

using DataFrames, Statistics

using ADRIA, CoralBlox

include("common.jl")

context_layers = find_latest_file("../canonical-reefs/output/")
context_layers = GDF.read(context_layers)

rs = ADRIA.load_results("outputs/ADRIA-out/ReefMod Engine__RCPs_45__2024-06-19_10_58_55_647")

fig_opts = Dict(:size => (1600, 800))



# s_tac = ADRIA.metrics.scenario_total_cover(rs)
tac = ADRIA.metrics.total_absolute_cover(rs)

# reduce all the scenarios down to one series for each reef
tac_sites = mapslices_toFloat64(median, tac, :scenarios)

# Have to remove the first year as there seems to be an issue with that year's data
tac_sites_reduced = tac_sites[timesteps=2:79]

# calculate the relative site cover from the initial cover across timesteps for each reef
rel_cover = mapslices_toFloat64(relative_site_cover, tac_sites_reduced, :timesteps)
# series(rel_cover.data', solid_color=(:black, 0.05))
# hlines!(1, color=:red)


test = rel_cover[:,1]

function delta_cover(x)
    y = Vector()
    for i in eachindex(x)
        if i < lastindex(x)
            push!(y, x[i+1] - x[i])
        end
    end

    return y
end



test_delta = Float64.(mapslices(delta_cover, rel_cover, dims=[:timesteps]))
series(test_delta.data', solid_color=:black)
hlines!(0, color=:red)


series(rel_cover.data', solid_color=:black)


df = DataFrame(rel_cover.data, collect(getAxis("sites", rel_cover).val))

df.year = [string(i+1) for i in 1:size(df,1)]
select!(df, :year, Not(:year))

data = permutedims(df, 1, "UNIQUE_ID")

data = leftjoin(data, context_layers[:, [:UNIQUE_ID, :management_area]], on=:UNIQUE_ID)
plot_lines(data, :management_area, 2, 79, "Year", "Proportion of reef cover Year2", 0.01)
hlines!(1, color=:green)
