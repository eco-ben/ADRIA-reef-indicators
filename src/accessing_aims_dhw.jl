using GLMakie

include("common.jl")
include("DHW_processing_functions.jl")

aims_dhw = open_dataset("../../DHW data/Jessica-Benthuysen/CoralSea_GBR_CNRM-ESM2-1_ssp245_r1i1p1f2_dhw_2015-2100.nc")
reef_outlines = GDF.read(find_latest_file("../../canonical-reefs/output/"))


reef_timeseries = extract_reef_DHW_timeseries(aims_dhw.dhw_max[:,:,:], reef_outlines, aims_dhw)
recreate_adriadomain_dhw_dataset(
    "../../ADRIA Domains/GBR_2024_10_15/DHWs/dhwRCP45.nc",
    reef_outlines,
    reef_timeseries,
    "../../ADRIA Domains/GBR_2024_10_15_Jessica/DHWs/dhwRCP45.nc"
)