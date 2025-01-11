using YAXArrays
using NetCDF

test_gbr_sst = open_dataset("../../HighResCoralStress_sst_EC-Earth3-Veg_ssp245_Australia_01011985_31122100_v2.0.nc")

lat = readcubedata(test_gbr_sst.lat)
lon = readcubedata(test_gbr_sst.lon)

# Extract rough GBR region
target_lat = (lat.data .< -9) .& (lat.data .> -25)
target_lon = (lon.data .< 155) .& (lon.data .> 140)
target_pixels = target_lat .& target_lon

target_lat_vals = lat[target_pixels]
target_lon_vals = lon[target_pixels]

# Extract 2025-2099 daily data
target_time = (test_gbr_sst.sst.Ti .> Date("2025", "yyyy-mm-dd")) .& (test_gbr_sst.sst.Ti .< Date("2100", "yyyy-mm-dd"))

sst_loc1 = test_gbr_sst.sst[pixels = 1]
historical_Ti = (test_gbr_sst.sst.Ti .>= Date("1985", "yyyy-mm-dd")) .& (test_gbr_sst.sst.Ti .< Date("2013", "yyyy-mm-dd"))
historical_sst = @view(sst_loc1[Ti = historical_Ti].data)
historical_Ti = sst_loc1[Ti = historical_Ti].Ti

MMM = zeros(Float64, (2013-1985) * 12)

for (i, year_i) in enumerate(1985:1:2013-1)
    year_sst = historical_sst[year.(historical_Ti) .== year_i]
    year_ti = historical_Ti[year.(historical_Ti) .== year_i]
    # Need to input the linear regression (Heron et al., 2014 Improvements to and continuity of operational global thermal stress monitoring for coral bleaching.).
    for j in 1:12
        MMM[((i-1)*12)+j] = mean(year_sst[month.(year_ti) .== j])
    end
end



# Subset and write out to .nc file
subset_gbr = test_gbr_sst.sst[pixels = target_pixels, Ti = target_time]
arrays = Dict(:lon => target_lon_vals, :lat => target_lat_vals, :sst => subset_gbr)
ds = Dataset(; arrays...)
savedataset(ds, path="../../subset_gbr_sst.nc", driver=:netcdf)