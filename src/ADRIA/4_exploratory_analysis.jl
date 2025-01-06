"""
Access the domain connectivity and calculate the strength and number of outgoing/incoming
connections for each reef. Run immediately after 4_.jl to use the updated context_layers with
bioregions.
"""

using DataFrames, Statistics, YAXArrays

using GLMakie, GeoMakie, GraphMakie

using GLM, MixedModels

using ADRIA

import GeoDataFrames as GDF

# If previous scripts have never been run or have been updated:
include("1_initial_results.jl");
include("2_Subregion_correlation.jl");
include("3_analysis_context_layers.jl")

# Load updated context layer data
context_layers = GDF.read("data/analysis_context_layers.gpkg",)

new_target_reefs_context = context_layers[context_layers.target_reefs , :]
old_target_reefs_context = context_layers[context_layers.target_reefs , :]
# Comparing the distribution of gbr wide connectivity metrics to those found in target_reefs

# hist(context_layers.total_comb; bins=[0,0.5,1,1.5,2,2.5,3,3.5,4,5,7.5,10,12.5,15,20,25,30,35,40,50,60,70,80,90,150])
# vlines!(1, color=:red)

# show(target_reefs_context.total_comb)
# hist(target_reefs_context.total_comb, bins=[0,0.5,1,1.5,2,2.5,3,3.5,4,5,7.5,10,12.5,15,20,25,30,35,40,50,60,70,80,90,150])
# vlines!(1, color=:red)
# hist(context_layers.total_connect)
# hist(target_reefs_context.total_connect)

# hist(context_layers.income_strength)
# hist(target_reefs_context.income_strength)

# hist(context_layers.income_count)
# hist(target_reefs_context.income_count)

# hist(context_layers_context.initial_coral_cover)
# show(target_reefs_context.initial_coral_cover)
# hist(target_reefs_context.initial_coral_cover)

# hist(context_layers_context[(context_layers.closest_port .== "Mourilyan Harbour"),:income_strength])
# hist(target_reefs_context[(target_reefs_context.closest_port .== "Mourilyan Harbour"),:income_strength])

# hist(context_layers[(context_layers.closest_port .== "Quintell Beach"),:total_comb], bins = 40)
# hist(target_reefs_context[(target_reefs_context.closest_port .== "Quintell Beach"),:total_comb])

# hist(context_layers[(context_layers.bioregion .== "Far Northern Protected Mid Shelf Reefs and Shoals"),:total_comb])
# hist(target_reefs_context[(target_reefs_context.bioregion .== "Far Northern Protected Mid Shelf Reefs and Shoals"),:total_comb])

# hist(context_layers[(context_layers.bioregion .== "Far Northern Outer Mid Shelf Reefs"),:total_comb])
# hist(target_reefs_context[(target_reefs_context.bioregion .== "Far Northern Outer Mid Shelf Reefs"),:total_comb])

# hist(context_layers[(context_layers.bioregion .== "Outer Shelf Reefs"),:total_comb])
# hist(target_reefs_context[(target_reefs_context.bioregion .== "Outer Shelf Reefs"),:total_comb])



# Load rel_cover data and identify reefs that crash in the first timeseries year
#rs = ADRIA.load_results("outputs/ADRIA-out/ReefMod Engine__RCPs_45__2024-06-19_10_58_55_647")
rs = ADRIA.load_results("outputs/ADRIA-out/ReefMod Engine__RCPs_45__2024-07-26_12_05_07_844")
tac = ADRIA.metrics.total_absolute_cover(rs)
tac_sites = Float64.(mapslices(median, tac, dims=[:scenarios]))
tac_sites_reduced = tac_sites[timesteps=2:79]
rel_cover = Float64.(mapslices(relative_site_cover, tac_sites_reduced, dims=[:timesteps]))
crash_reefs = collect(getAxis("sites", rel_cover).val)[rel_cover.data[2,:] .< 0.5]

# Remove crashing reefs
no_crash_gbr = context_layers[context_layers.UNIQUE_ID .∉ [crash_reefs],:]

# Filter to only include reefs that are within the same bioregion/closest_port subregion as target reefs
filtered_bior = no_crash_gbr[no_crash_gbr.bioregion .∈ [unique(no_crash_gbr[no_crash_gbr.target_reefs_bior, :bioregion])], :]
filtered_subr = no_crash_gbr[no_crash_gbr.closest_port .∈ [unique(no_crash_gbr[no_crash_gbr.target_reefs_subr, :closest_port])], :]

# Basic exploratory models of factors
glm_allreefs = glm(@formula(target_reefs ~ out_comb + dhw_cor ), no_crash_gbr, Binomial())
glm_bioregions = glm(@formula(target_reefs_bior ~ initial_coral_cover + out_comb + 0), filtered_bior, Binomial())
glm_subregions = glm(@formula(target_reefs_subr ~ total_comb + so_to_si + initial_coral_cover + bioregion_count), filtered_subr, Binomial())

aic(glm_allreefs), aic(glm_bioregions), aic(glm_subregions)

glmm_form = @formula(target_reefs ~ out_comb + income_comb + total_comb + dhw_cor + mean_dhw + (1|bioregion))
glmm_fit = fit(MixedModel, glmm_form, no_crash_gbr, Binomial(), ProbitLink())
aic(glmm_fit) # seems to be an improvement when (1|bioregion) is used

# Initial plot of target_reefs and their surrounding reefs from filtered_bior
df = DataFrame(rel_cover.data, collect(getAxis("sites", rel_cover).val))
df.year = [string(i+1) for i in 1:size(df,1)]
select!(df, :year, Not(:year))

data = permutedims(df, 1, "UNIQUE_ID")
data = data[data.UNIQUE_ID .∈ [filtered_bior.UNIQUE_ID],:]
data = leftjoin(data, filtered_bior[:, [:UNIQUE_ID, :target_reefs, :bioregion, :management_area]], on=:UNIQUE_ID)
data.target_reefs = ifelse.(data.target_reefs, "target", "non-target")
target_reefs = Matrix(DataFrames.select(data[data.target_reefs .== "target", :], DataFrames.Between(Symbol(2), Symbol(15))))'
target_reefs = mean(target_reefs, dims=2)
non_target_reefs = Matrix(DataFrames.select(data[data.target_reefs .== "non-target", :], DataFrames.Between(Symbol(2), Symbol(15))))'
non_target_reefs = mean(non_target_reefs, dims=2)

f, ax = plot_lines(data, :management_area, 2, 15, "Year", "Proportion of reef cover Year2", 0.2)
series!(target_reefs', solid_color=:red)
series!(non_target_reefs', solid_color=:blue)
f, ax = plot_lines(data[data.UNIQUE_ID .∈ [target_reefs_context.UNIQUE_ID],:], :bioregion, 2, 15, "Year", "Proportion of reef cover Year2", 0.2)


# Looking within a single bioregion and visualising the lagged target timeseries
unique(filtered_bior.bioregion)
bioregion_only = filtered_bior[(filtered_bior.bioregion .== ["Far Northern Protected Mid Shelf Reefs and Shoals"]),:]
df = DataFrame(rel_cover.data, collect(getAxis("sites", rel_cover).val))
df.year = [string(i+1) for i in 1:size(df,1)]
select!(df, :year, Not(:year))

data = permutedims(df, 1, "UNIQUE_ID")
data = data[data.UNIQUE_ID .∈ [filtered_bior.UNIQUE_ID],:]
data = leftjoin(data, filtered_bior[:, [:UNIQUE_ID, :target_reefs]], on=:UNIQUE_ID)
data.target_cat = ifelse.(data.target_reefs, "target", "non-target")

lagged_data = data[data.target_cat .== ["non-target"],:]
select!(lagged_data, :UNIQUE_ID, Not(:2, :3, :4, :5, :6))
#rename!(lagged_data, [:UNIQUE_ID; [Symbol(i) for i in 7:84]; :target_cat])

#lagged_data_mat = zeros(Union{Float64, Missing}, 78, 116)
lagged_data_mat = Matrix(DataFrames.select(lagged_data, DataFrames.Between(Symbol(7), Symbol(79))))'
#lagged_data_mat[1:5, :] .= missing
lagged_data_mat = lagged_data_mat[1:15,:]

test_mean_lagged = mean(lagged_data_mat, dims=2)
test_mean_target = mean(Matrix(DataFrames.select(data[data.target_cat .== ["target"],:], DataFrames.Between(Symbol(2), Symbol(15)))), dims=1)

f, ax = plot_lines(data[data.target_cat .== ["target"],:], :target_cat, 2, 15, "Year", "Proportion of reef cover Year2", 0.2)
series!(ax, lagged_data_mat', solid_color=(:orange, 0.2))
series!(ax, test_mean_lagged', solid_color=:red)
series!(ax, test_mean_target, solid_color=:blue)



# Try plotting with long dataframe format - not currently working as intended
FNOpen_LagoonR = filtered_bior[(filtered_bior.bioregion .== ["Far Northern Open Lagoon Reefs"]),:]
FNOL_targets = FNOpen_LagoonR[FNOpen_LagoonR.target_reefs, :UNIQUE_ID]
mean(lagged_analysis_bior[lagged_analysis_bior.RME_UNIQUE_ID .∈ [FNOL_targets],:lag6])

df.year = parse.(Int,df.year)

test = stack(df, names(df[:,2:3807]))
rename!(test, :variable => :UNIQUE_ID)
lines(test.year, test.value)

test_2 = test[test.year .∉ [[2, 3,4,5,6,7,8,9]],:]
lines!(test_2.year, test_2.value)
