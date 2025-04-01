using YAXArrays
using Statistics
using NetCDF
using Dates
using GLM
using ProgressMeter
import GeoDataFrames as GDF
using Rasters

include("common.jl")

# Rolling sum function 
function rolling_sum(a, n::Int)
    @assert 1<=n<=length(a)
    out = similar(a, length(a)-n+1)
    out[1] = sum(a[1:n])
    for i in eachindex(out)[2:end]
        out[i] = out[i-1]-a[i-1]+a[i+n-1]
    end
    return out
end

"""
    identify_search_pixels(input_raster::Raster, criteria_function)::DataFrame

Identifies all pixels in an input raster that return true for the function `criteria_function`.

# Arguments
- `input_raster` : Raster containing pixels for the target region.
- `criteria_function` : Function that returns a boolean value for each pixel in `input_raster`.
                        Pixels that return true will be targetted in analysis.

# Returns
DataFrame containing indices, lon and lat for each pixel that is intended for further analysis.
"""
function identify_search_pixels(input_raster::Raster, criteria_function)::DataFrame
    pixels = map(criteria_function, input_raster)
    indices::Vector{CartesianIndex{2}} = findall(pixels)
    indices_lon::Vector{Float64} = lookup(pixels, X)[first.(Tuple.(indices))]
    indices_lat::Vector{Float64} = lookup(pixels, Y)[last.(Tuple.(indices))]

    return DataFrame(; indices=indices, lons=indices_lon, lats=indices_lat)
end

"""
    subset_highrescoralstress_nc(nc_dataset::Dataset, reef_outlines::DataFrame)

Identify the target pixels from `nc_dataset` that touch `reef_outlines` polygons.
"""
function subset_highrescoralstress_nc(nc_dataset::Dataset, reef_outlines::DataFrame)
    lat = readcubedata(nc_dataset.lat)
    lon = readcubedata(nc_dataset.lon)
    all_coords = tuple.(lon.data, lat.data)
    
    # # Extract rough GBR region
    target_lat = (lat.data .< -9) .& (lat.data .> -25)
    target_lon = (lon.data .< 155) .& (lon.data .> 140)
    target_pixels = target_lat .& target_lon
    
    gbr_lon = sort(unique(lon[target_pixels].data))
    gbr_lat = sort(unique(lat[target_pixels].data); rev=true)
    
    base_raster = Raster(
        ones(
            Bool, 
            X(minimum(gbr_lon):0.01:maximum(gbr_lon); sampling=Rasters.Intervals(Rasters.Center())), 
            Y(maximum(gbr_lat):-0.01:minimum(gbr_lat); sampling=Rasters.Intervals(Rasters.Center()))
        )
    )
    base_raster = Rasters.mask(
        Rasters.crop(base_raster; to=reef_outlines.geometry); 
        with=reef_outlines.geometry, 
        boundary=:touches,
        missingval=0
    )
    base_raster = replace_missing(base_raster, 0)
    
    outline_pixels = identify_search_pixels(Bool.(base_raster), x -> x)
    outline_pixels = tuple.(outline_pixels.lons, outline_pixels.lats)
    target_pixels = all_coords .∈ [outline_pixels]
    
    return target_pixels
end

"""
    calculate_annual_dhw(nc_dataset::Dataset, target_pixels::BitVector)::YAXArray

Calculate annual DHW timeseries for each pixel in `nc_dataset.sst` that is selected 
by `target_pixels`. 
The method for DHW calculation follows NOAA DHW calculation, this includes:
- The use of historical period 1985-2012 (inclusive)
- Recentering the monthly mean data to 1988.2857 by using regression.
- Calculating DHW as the rolling sum of HotSpot heat stress anomaly values from past 12 weeks.

# Arguments
- `nc_dataset` : Daily SST dataset with sst, lon and lat variables. 
- `target_pixels` : Vector containing the target pixels for calculation from nc_dataset.

# References
- MMM climatology methodology : https://www.ncei.noaa.gov/data/oceans/coris/library/NOAA/CRCP/project/915/TR_NESDIS_145.pdf
- HotSpot and DHW calculation methodology (and MMM) : https://coralreefwatch.noaa.gov/product/5km/methodology.php#ref_her
"""
function calculate_annual_dhw(nc_dataset::Dataset, target_pixels::BitVector)::YAXArray
    gbr_sst = nc_dataset.sst[pixels = target_pixels]
    required_lon = nc_dataset.lon[target_pixels]
    required_lat = nc_dataset.lat[target_pixels]
    gbr_sst = gbr_sst[:, :]

    gbr_annual_DHW = zeros(Float64, (first(size(gbr_sst)), length(1985:2100)))

    unique_pixels = gbr_sst.pixels

    @showprogress for (p, pixel) in enumerate(unique_pixels)

        sst_loc = gbr_sst[pixels = At(pixel)]
        historical_Ti = (nc_dataset.sst.Ti .>= Date("1985", "yyyy-mm-dd")) .& (nc_dataset.sst.Ti .< Date("2013", "yyyy-mm-dd"))
        historical_sst = sst_loc[Ti = historical_Ti].data
        historical_Ti = sst_loc[Ti = historical_Ti].Ti

        hist_years = 1985:2012

        # Calculating the Monthly Means
        MMs = zeros(Float64, 12)

        # j refers to a unique month and i refers to a unique year.
        for j in 1:12

            month_sst = historical_sst[month.(historical_Ti) .== j]
            month_Ti = historical_Ti[month.(historical_Ti) .== j]
            MM_j = zeros(Float64, length(hist_years))

            for (i, year_i) in enumerate(hist_years)
                # Calculate the mean SST for year_i and month j
                year_sst = month_sst[year.(month_Ti) .== year_i]
                MM_j[i] = mean(year_sst)
            end

            coefficients = coef(lm(@formula(MM_j ~ hist_years), (;hist_years, MM_j)))
            intercept_j, slope_j = (first(coefficients), last(coefficients))
            pred_j = (slope_j * 1988.2857) + intercept_j
            MMs[j] = pred_j
        end

        # Maximum monthly mean
        MMM = maximum(MMs)

        # HotSpot SST anomaly
        # Only consider HotSpot heat stress values greater than 1! 
        # Divide by 7 to express DHW values in degree Celsius weeks. 
        Hs = sst_loc.data .- MMM
        Hs = ifelse.(Hs .< 1, 0, Hs) ./ 7 # !This must be set to  Hs < 1 according to NOAA methodology!

        # DHW calculation
        n = 7 * 12
        DHW = rolling_sum(Hs, n) # Calculate DHW as rolling sum of last 12 weeks.
        timesteps = sst_loc.Ti[n:end]

        annual_DHW = zeros(Float64, length(1985:2100))

        for (i, year_i) in enumerate(1985:2100)

            year_daily_dhw = DHW[year.(timesteps) .== year_i]
            year_timesteps = timesteps[year.(timesteps) .== year_i]
            month_mean_DHWs = [mean(year_daily_dhw[month.(year_timesteps) .== j]) for j in 1:12]
            month_mean_DHWs = ifelse.(isnan.(month_mean_DHWs), -Inf, month_mean_DHWs)

            annual_DHW[i] = maximum(month_mean_DHWs)
        end

        gbr_annual_DHW[p, :] = annual_DHW

    end

    axlist = (
        gbr_sst.pixels,
        Dim{:timesteps}(1985:2100)
    )
    yax_annual_DHW = YAXArray(axlist, gbr_annual_DHW)

    return yax_annual_DHW
end

"""
    extract_reef_DHW_timeseries(
        annual_DHW_pixels::YAXArray, 
        target_reef_outlines::DataFrame, 
        target_pixels::BitVector, 
        source_nc::Dataset
    )::YAXArray

Extract mean DHW values for each reef polygon in `target_reef_outlines`.
"""
function extract_reef_DHW_timeseries(
    annual_DHW_pixels::YAXArray, 
    target_reef_outlines::DataFrame, 
    target_pixels::BitVector, 
    source_nc::Dataset
)::YAXArray

    # Converting DHW timeseries to raster objects to extract mean values for a reef outline
    axlist = (
        Dim{:locations}(target_reef_outlines.UNIQUE_ID),
        Dim{:timesteps}(2025:2099)
    )
    loc_DHW_timeseries = YAXArray(axlist, Matrix{Union{Missing, Float64}}(missing, (length(target_reef_outlines.UNIQUE_ID), length(2025:2099))))

    required_lon = source_nc.lon[target_pixels]
    required_lat = source_nc.lat[target_pixels]

    unique_lon = sort(unique(required_lon))
    unique_lat = sort(unique(required_lat); rev=true)

    for (t, time) in enumerate(2025:2099)
        time_raster = Raster(
            zeros(
                Float64, 
                X(minimum(unique_lon):0.01:maximum(unique_lon); sampling=Rasters.Intervals(Rasters.Center())), 
                Y(maximum(unique_lat):-0.01:minimum(unique_lat); sampling=Rasters.Intervals(Rasters.Center()))
            )
        )
        time_raster .= -9999.0
        time_raster = rebuild(time_raster, missingval=-9999.0)

        for p in eachindex(annual_DHW_pixels.pixels)
            lon_p = required_lon.data[p]
            lat_p = required_lat.data[p]
            DHW = annual_DHW_pixels[timesteps = At(time)][p]
            time_raster[X = At(lon_p), Y = At(lat_p)] = DHW
        end

        for (l, loc) in enumerate(loc_DHW_timeseries.locations)
            println("extracting DHW data for time $(time) location $(l)")
            geom = target_reef_outlines[target_reef_outlines.UNIQUE_ID .== loc, :geometry]
            mean_dhw = Float64(first(zonal(mean, time_raster; of=geom, boundary=:touches)));

            loc_DHW_timeseries[locations = At(loc)][t] = mean_dhw
        end
    end

    return loc_DHW_timeseries
end

function extract_reef_DHW_timeseries(annual_DHW_pixels, target_reef_outlines, source_nc; pixel_step=0.05)

    # Converting DHW timeseries to raster objects to extract mean values for a reef outline
    axlist = (
        Dim{:locations}(target_reef_outlines.UNIQUE_ID),
        Dim{:timesteps}(2025:2099)
    )
    loc_DHW_timeseries = YAXArray(axlist, Matrix{Union{Missing, Float64}}(missing, (length(target_reef_outlines.UNIQUE_ID), length(2025:2099))))

    lon_max = maximum(round.(source_nc.lon.data; digits=3))
    lon_min = minimum(round.(source_nc.lon.data; digits=3))
    lat_max = maximum(round.(source_nc.lat.data; digits=3))
    lat_min = minimum(round.(source_nc.lat.data; digits=3))

    unique_lon = X(lon_min:pixel_step:lon_max; sampling=Rasters.Intervals(Rasters.Center()))
    unique_lat = Y(lat_max:-pixel_step:lat_min; sampling=Rasters.Intervals(Rasters.Center()))
    if (length(unique_lon), length(unique_lat)) != size(annual_DHW_pixels[:,:,1].data)
        throw("Error, lon/lat coordinates are likely not evenly distributed to create a raster grid.")
    end

    missing_val = annual_DHW_pixels.properties["_FillValue"]

    for (t, time) in enumerate(2025:2099)
        # time_raster = Raster(
        #     zeros(
        #         Float64, 
        #         X(minimum(unique_lon):0.01:maximum(unique_lon); sampling=Rasters.Intervals(Rasters.Center())), 
        #         Y(maximum(unique_lat):-0.01:minimum(unique_lat); sampling=Rasters.Intervals(Rasters.Center()))
        #     )
        # )
        # time_raster .= -1.0e20
        # time_raster = rebuild(time_raster, missingval=-9999.0)

        # for p in eachindex(annual_DHW_pixels.pixels)
        #     lon_p = required_lon.data[p]
        #     lat_p = required_lat.data[p]
        #     DHW = annual_DHW_pixels[timesteps = At(time)][p]
        #     time_raster[X = At(lon_p), Y = At(lat_p)] = DHW
        # end
        time_index = findfirst(year.(annual_DHW_pixels.Ti) .== time)
        time_raster = Raster(
            reverse(annual_DHW_pixels[Ti = time_index].data; dims=2), 
            (unique_lon, unique_lat); 
            missingval=missing_val
        )
        time_raster = replace_missing(time_raster; missingval=missing)

        for (l, loc) in enumerate(loc_DHW_timeseries.locations)
            println("extracting DHW data for time $(time) location $(l)")
            geom = target_reef_outlines[target_reef_outlines.UNIQUE_ID .== loc, :geometry]
            mean_dhw = Float64(first(zonal(mean, time_raster; of=geom, boundary=:touches)));

            loc_DHW_timeseries[locations = At(loc)][t] = mean_dhw
        end
    end

    return loc_DHW_timeseries
end

function recreate_adriadomain_dhw_dataset(
    template_gbr_dhw_nc::String,
    reef_outlines::DataFrame, 
    reef_DHW_timeseries::YAXArray, 
    result_nc::String
)::Nothing
    GBR_template = open_dataset(template_gbr_dhw_nc)
    template_dhw = GBR_template.dhw

    longitude = GBR_template.longitude
    latitude = GBR_template.latitude
    ax = (Dim{:sites}(1:1:3806),)
    UNIQUE_ID = YAXArray(ax, reef_outlines.UNIQUE_ID)

    axlist = (
        Dim{:timesteps}(2025:1:2099),
        Dim{:sites}(1:1:3806),
        Dim{:member}(1:2)
    )

    loc_DHW_timeseries_data = reshape(reef_DHW_timeseries.data', (75, 3806, 1))
    loc_DHW_timeseries_data = cat(loc_DHW_timeseries_data, loc_DHW_timeseries_data; dims=3)
    loc_DHW_timeseries_data = YAXArray(axlist, loc_DHW_timeseries_data, template_dhw.properties)

    data_arrays = Dict(:longitude => longitude, :latitude => latitude, :dhw => Float64.(loc_DHW_timeseries_data))
    output_dhw_data = Dataset(; properties=GBR_template.properties, data_arrays...)
    savedataset(output_dhw_data, path=result_nc, driver=:netcdf, overwrite=true)

    return nothing
end

function recreate_adriadomain_dhw_yaxarray(
    template_gbr_dhw_nc::String, 
    reef_DHW_timeseries::YAXArray
)::YAXArray
    GBR_template = open_dataset(template_gbr_dhw_nc)
    template_dhw_properties = GBR_template.dhw.properties

    axlist = (
        Dim{:timesteps}(2025:1:2099),
        Dim{:sites}(1:1:3806),
        #Dim{:member}([1])
    )

    loc_DHW_timeseries_data = reef_DHW_timeseries.data #loc_DHW_timeseries_data = reshape(reef_DHW_timeseries.data', (75, 3806, 1))
    loc_DHW_timeseries_data = YAXArray(axlist, loc_DHW_timeseries_data', template_dhw_properties)

    return loc_DHW_timeseries_data
end