include("../../common.jl")

context_layers = GDF.read(joinpath(output_path, "analysis_context_layers_carbonate.gpkg"))

dhw_scenarios = open_dataset(joinpath(gbr_domain_path, "DHWs/dhwRCP45.nc"))
GCMs = dhw_scenarios.dhw.properties["members"]

gbr_dom = ADRIA.load_domain(gbr_domain_path, "45")
gbr_dom_filtered = gbr_dom.loc_data[gbr_dom.loc_data.UNIQUE_ID .∈ [context_layers.UNIQUE_ID], :]
filtered_indices = indexin(gbr_dom_filtered.UNIQUE_ID, gbr_dom.loc_data.UNIQUE_ID)
areas = gbr_dom.loc_data.area

all_levels = vec(collect(Iterators.product(GCMs, unique(context_layers.management_area))))
distances = DataFrame(GCM = first.(all_levels), management_area = last.(all_levels), CID_cover = zeros(length(all_levels)), CID_dhw = zeros(length(all_levels)))

for (i_gcm, GCM) in enumerate(GCMs)

    absolute_cover = readcubedata(open_dataset(joinpath(output_path, "processed_model_outputs/median_cover_$(GCM).nc")).layer)
    relative_cover = percentage_cover_timeseries(areas, absolute_cover)[1:50, :]
    relative_cover = relative_cover[locations=At(context_layers.UNIQUE_ID)]

    dhw_ts = gbr_dom.dhw_scens[1:50, filtered_indices, i_gcm]
    dhw_ts = rebuild(dhw_ts, dims=relative_cover.axes, metadata=dhw_timeseries_properties)

    for man_area in unique(context_layers.management_area)
        man_area_reefs = context_layers[context_layers.management_area .== man_area, :]

        low_cluster_reefs = man_area_reefs[man_area_reefs[:, "$(GCM)_management_area_cluster_cats"] .== "low", :]
        low_cluster_cover = relative_cover[locations=At(low_cluster_reefs.UNIQUE_ID)]
        low_cluster_dhw = dhw_ts[locations=At(low_cluster_reefs.UNIQUE_ID)]

        high_cluster_reefs = man_area_reefs[man_area_reefs[:, "$(GCM)_management_area_cluster_cats"] .== "high", :]
        high_cluster_cover = relative_cover[locations=At(high_cluster_reefs.UNIQUE_ID)]
        high_cluster_dhw = dhw_ts[locations=At(high_cluster_reefs.UNIQUE_ID)]

        man_area_cid_cover = CID(
            dropdims(median(low_cluster_cover, dims=2), dims=2),
            dropdims(median(high_cluster_cover, dims=2), dims=2)
        )
        man_area_cid_dhw = CID(
            dropdims(median(low_cluster_dhw, dims=2), dims=2),
            dropdims(median(high_cluster_dhw, dims=2), dims=2)
        )

        distances[
            (distances.GCM .== GCM) .& 
            (distances.management_area .== man_area),
        :CID_cover] .= round(man_area_cid_cover, digits=2) 

        distances[
            (distances.GCM .== GCM) .& 
            (distances.management_area .== man_area),
        :CID_dhw] .= round(man_area_cid_dhw, digits=2)
    end
end

sort!(distances, [:management_area, :GCM])
distances.management_area = replace.(distances.management_area, ["Management Area" => ""])
