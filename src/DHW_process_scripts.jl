
include("common.jl")
include("DHW_processing_functions.jl")

# Required files:
reef_outlines = GDF.read(find_latest_file("../../canonical-reefs/output/"))
ADRIA_template_fn = "../../ADRIA Domains/GBR_2024_10_15/DHWs/dhwRCP45.nc" # ADRIA template dhw.nc file
output_dhw_dataset_fn = "../../ADRIA Domains/GBR_2024_10_15_HighResCoralStress/DHWs/dhwRCP45.nc" # ADRIADomain location to output final DHW file

GCMs = [
    "ACCESS-ESM1-5",
    "ACCESS-CM2",
    "EC-Earth3-Veg",
    "GFDL-CM4",
    "NorESM2-MM"
]
DHW_data = Dict{String, Any}(gcm => missing for gcm in GCMs)

for gcm in GCMs
    # Load original Australia region SST data from HighResCoralStress
    sst_dataset_fn = "../../DHW data/HighResCoralStress/HighResCoralStress_sst_$(gcm)_ssp245_Australia_01011985_31122100_v2.0.nc"
    gcm_output_fn = "../outputs/$(gcm)_dhw_trajectories.nc"

    # Check if GCM DHW already exists
    if isfile(gcm_output_fn)
        @info "$(gcm) already processed, continuing to next GCM"
        DHW_data[gcm] = open_dataset(gcm_output_fn).dhw
        continue
    end

    aus_sst = open_dataset(sst_dataset_fn)

    # Identify pixels that touch reef polygons
    gbr_reef_pixels = subset_highrescoralstress_nc(aus_sst, reef_outlines)

    # Calculate annual DHW trajectories for each relevant pixel
    annual_dhw_data = calculate_annual_dhw(aus_sst, gbr_reef_pixels)

    # Extract timeseries for reef polygons
    reef_timeseries = extract_reef_DHW_timeseries(
        annual_dhw_data, 
        reef_outlines, 
        gbr_reef_pixels, 
        aus_sst
    )

    # Conform to expected format of DHW arrays in ADRIADomain
    DHW_data[gcm] = recreate_adriadomain_dhw_yaxarray(
        ADRIA_template_fn,
        reef_timeseries
    )
    savecube(DHW_data[gcm], gcm_output_fn, driver=:netcdf) # Save intermediate GCM file
end

# Combine GCM YAXArrays and add GCM labels in properties
cubes = [
    DHW_data["EC-Earth3-Veg"],
    DHW_data["ACCESS-ESM1-5"],
    DHW_data["ACCESS-CM2"],
    DHW_data["NorESM2-MM"],
    DHW_data["GFDL-CM4"]
]
member_axis = Dim{:member}(1:5)
combined_dhw = concatenatecubes(cubes, member_axis)

axlist = (
    Dim{:timesteps}(2025:1:2099),
    Dim{:sites}(1:1:3806),
    Dim{:member}(1:1:5)
)
props = combined_dhw.properties
props["members"] = [
    "EC-Earth3-Veg", 
    "ACCESS-ESM1-5", 
    "ACCESS-CM2",
    "NorESM2-MM",
    "GFDL-CM4"
]
combined_dhw = rebuild(combined_dhw, dims=axlist, properties=props)

# Recreate final ADRIADomain DHW timeseries file
GBR_template = open_dataset(ADRIA_template_fn)
longitude = GBR_template.longitude
latitude = GBR_template.latitude

data_arrays = Dict(:longitude => longitude, :latitude => latitude, :dhw => Float64.(combined_dhw))

new_properties = GBR_template.properties
new_properties["source_path"] = "DHW data/HighResCoralStress/"
new_properties["source_location"] = "highrescoralstress.org"
new_properties["source_desc"] = "DHW calculated from HighResCoralStress downscaled daily SST data."

output_dhw_data = Dataset(; properties = new_properties, data_arrays...)
# Save DHW timeseries file into an ADRIADomain location that is set up.
savedataset(output_dhw_data, path = output_dhw_dataset_fn, driver=:netcdf, overwrite=true)

