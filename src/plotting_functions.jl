"""
File includes the plotting functions for ADRIA-reef-indicators analysis.
"""

using Printf
using Colors
using Random
using TOML
import GeometryOps as GO


config = TOML.parsefile("config.toml")

# Define figure formatting constants
fontsize = 6
dpi = 300

# Define figure sizing constants for different figure types
# Maximum figure size in word document with 2.54cm margins is 15.9cm (if images are larger, word will shrink them)
fig_sizes = Dict{String,Union{Float64,Int64}}(
    "violin_width" => 14.82,
    "carb_width" => 14.82,
    "map_height" => 14.82,
    "map_width" => 14.82,
    "ecs_width" => 10,
    "ecs_height" => 13,
    "gcm_timeseries_width" => 14.82,
    "violin_height" => 12,
    "carb_height" => 11,
    "gcm_timeseries_height" => 12,
    "cluster_hm_width" => 14.82,
    "cluster_hm_height" => 18,
    "grouped_timeseries_width" => 14.82,
    "grouped_timeseries_height" => 14.82
)

# Convert figure sizes from cm to pixel measurement

# Size of 1cm in pixels relative to 1 CSS px, see:
# - https://docs.makie.org/dev/how-to/match-figure-size-font-sizes-and-dpi
# - https://docs.makie.org/dev/explanations/figure#Figure-size-and-resolution
centimetre = 37.7952755906
map!(x -> x * centimetre, values(fig_sizes))

# Convert fontsize to pixel measurement
pt = 1.33 # size of 1 pt in pixels
fontsize = fontsize * pt

inch = 96 # size of 1 inch in pixels
dpi = dpi / inch

function colorscheme_alpha(cscheme::ColorScheme, alpha=0.5)
    ncolors = length(cscheme)

    return ColorScheme([RGBA(get(cscheme, k), alpha) for k in range(0, 1, length=ncolors)])
end

function skipmissing_median(x)
    new_median = Vector{Union{Missing,Float64}}(missing, size(x, 1))
    for (i, row) in enumerate(eachrow(x))
        if any(ismissing.(row))
            continue
        end
        new_median[i] = median(row)
    end

    return new_median
end

"""
    label_lines(label::String; l_length=25)

Separate a sentence label into multiple lines at l_length point, or the nearest space
between words.
"""
function label_lines(label::String; l_length=25)
    if length(label) > l_length
        label_new = replace(label, Regex("(.{$(l_length)} )") => s"\1\n")

        return label_new
    end

    return label
end

function _axis_size(gdf, x_fig_size, y_fig_size, n_col; x_gap=1.5, y_gap=1.2)
    xsize = x_fig_size / (n_col * x_gap)
    n_fig_row = first(fldmod1(length(gdf), n_col))
    ysize = y_fig_size / (n_fig_row * y_gap)

    return xsize, ysize
end

function _setup_grouped_figure(dataframe, grouping; x_fig_size=2130, y_fig_size=1500, marker=LineElement, order="$(grouping)_average_latitude", multi_axis=true, fontsize=11pt, rev=true)
    dataframe = sort(dataframe, order; rev=rev)
    gdf = DataFrames.groupby(dataframe, grouping)

    fig = Figure(size=(x_fig_size, y_fig_size), fontsize=fontsize)
    plot_layout = nothing

    if multi_axis
        n_col = optimum_columns(size(unique(dataframe[:, grouping]), 1))
        plot_layout = figure_layouts(length(gdf), n_col)

        # Configure the legend for the figure
        last_figure = last(plot_layout)
        if last(last_figure) < n_col
            legend_position = (last_figure[1], last_figure[2]+1:n_col)
            n_banks = 1
        else
            legend_position = (last_figure[1] + 1, 1:n_col)
            n_banks = 3
        end

        legend_entries = []
        for (i, col) in enumerate([:green, :orange, :blue])
            LE = marker(; color=col, marker=:circle)
            push!(legend_entries, [LE])
        end

        Legend(
            fig[legend_position...],
            legend_entries,
            ["Low", "Medium", "High"],
            nbanks=n_banks,
            tellwidth=false,
            padding=(2.0, 2.0, 2.0, 2.0),             # shrink padding inside legend box
            labelsize=fontsize,    # smaller font
            framevisible=false,        # optional: remove box
            markerlabelgap=3,          # reduce space between marker and label
            rowgap=0,                  # reduce vertical spacing between items
            colgap=4,                  # reduce horizontal spacing
            patchsize=(5, 5)
        )
    end

    return fig, gdf, plot_layout
end

function _setup_grouped_axes(fig, plot_layout_xi, xticks; ylabel="", xlabel="", title="", xsize=220, ysize=150, background_color=:white, xticklabelrotation=0.0, spinewidth=1.0)

    ax = Axis(
        fig[plot_layout_xi...];
        backgroundcolor=background_color,
        #xticklabelsize=16,
        #yticklabelsize=16,
        xticks=xticks,
        xticksvisible=false,
        title=title,
        width=xsize,
        height=ysize,
        spinewidth=spinewidth,
        xticklabelrotation=xticklabelrotation
    )
    return ax
end

function timeseries_xticks(length_t, years)
    length_range = first(length_t):10:last(length_t)
    years_length = collect(years)[length_range]

    return (length_range, string.(years_length))
end

function bootstrap_mean_ts(timeseries)
    ts_mean = Vector{Union{Missing,Float64}}(missing, size(timeseries, 1))
    ts_lb = Vector{Union{Missing,Float64}}(missing, size(timeseries, 1))
    ts_ub = Vector{Union{Missing,Float64}}(missing, size(timeseries, 1))

    for t in 1:size(timeseries, 1)
        if all(.!ismissing.(timeseries[t, :]))
            bootstrap_t = bootstrap(mean, timeseries[t, :], BasicSampling(1000))
            ts_mean[t] = first(bootstrap_t.t0)
            ts_lb[t] = getindex(getindex(confint(bootstrap_t, NormalConfInt(0.95)), 1), 2)
            ts_ub[t] = getindex(getindex(confint(bootstrap_t, NormalConfInt(0.95)), 1), 3)
        end
    end

    return ts_mean, ts_lb, ts_ub
end

function figure_layouts(n_bioregions, n_col)
    return [fldmod1(x, n_col) for x in 1:n_bioregions]
end

function optimum_columns(n_bioregions)
    if n_bioregions < 5
        return 2
    elseif n_bioregions < 10
        return 3
    elseif n_bioregions < 20
        return 4
    elseif n_bioregions < 30
        return 5
    end

    return 6
end

function GCM_label(x, y, GCM)
    return text!(x, y; text=GCM, align=(:right, :center))
end

"""
    ecs_plot(
        ecs_values::Vector{Union{Float64, Int64}},
        low_conf_range::Vector{Float64},
        high_conf_range::Vector{Float64},
        GCM_labels::Vector{String};
        fig_sizes::Dict=fig_sizes,
        fontsize=fontsize
    )

Create methods publication figure that displays the ECS values for each GCM investigated
and display the likely/very-likely ranges of ECS values from the IPCC.
"""
function ecs_plot(
    ecs_values::Vector{Float64},
    low_conf_range::Vector{Float64},
    high_conf_range::Vector{Float64},
    GCM_labels::Vector{String};
    fig_sizes::Dict=fig_sizes,
    fontsize=fontsize
)
    fig_x_size = fig_sizes["ecs_width"]
    fig_y_size = fig_sizes["ecs_height"]

    low_min, low_max = extrema(low_conf_range)
    high_min, high_max = extrema(high_conf_range)

    fig = Figure(size=(fig_x_size, fig_y_size), fontsize=fontsize)
    ax = Axis(
        fig[1, 1],
        ylabel="Equilibrium Climate Sensitivity (\u00B0C)",
        height=10cm,
        width=8cm
    )
    hidexdecorations!(ax)
    poly!(
        ax,
        [(0.98, high_min), (1.02, high_min), (1.02, high_max), (0.98, high_max)];
        color=(:red, 0.2), linestyle=:dash, strokewidth=1.5, label="very likely"
    )
    poly!(
        ax,
        [(0.98, low_min), (1.02, low_min), (1.02, low_max), (0.98, low_max)];
        color=(:red, 0.6), strokewidth=1.5, label="likely"
    )
    scatter!(ax, ones(length(ecs_values)), ecs_values, markersize=15, color=:black)
    GCM_label.(fill(1.06, length(ecs_values)), ecs_values, GCM_labels)

    fig[2, 1] = Legend(fig, ax, "ECS assessed range", framevisible=false, nbanks=2)

    return fig
end

"""
    grouped_cluster_timeseries_plots(
        timeseries_array::YAXArray,
        dataframe::DataFrame,
        cluster_col::Union{Symbol, String},
        grouping::Union{Symbol, String},
        length_t::UnitRange;
        fig_sizes::Dict=fig_sizes,
        fontsize=fontsize,
        ytitle::String=""
    )

Plot grouped-scale triple-cluster timeseries for each group. Uses ADRIA.viz.clustered_scenarios!().
Designed to plot with 3 clusters per group.

# Arguments
- `timeseries_array` : YAXArray containing timeseries, with dims at least nrow(dataframe) * length(length_t).
- `dataframe` : Dataframe containing cluster assignment and grouping information for each location.
- `cluster_col` : Column in `dataframe` containing numerical cluster assignments.
- `grouping` : Column in `dataframe` containing the unique grouping levels for locations.
- `length_t` : Indices to plot timeseries from `timeseries_array`, e.g. 1:50 for timesteps 1:50.
- `fig_sizes` : Figure sizing dictionary containing timeseries_ width and height measurements.
- `ytitle` : Y-axis label for timeseries plots.
"""
function grouped_cluster_timeseries_plots(
    timeseries_array::YAXArray,
    dataframe::DataFrame,
    cluster_col::Union{Symbol,String},
    grouping::Union{Symbol,String},
    length_t::UnitRange;
    fig_sizes::Dict=fig_sizes,
    fontsize=fontsize,
    ytitle::String=""
)
    fig_x_size = fig_sizes["grouped_timeseries_width"]
    fig_y_size = fig_sizes["grouped_timeseries_height"]
    n_col = optimum_columns(length(unique(dataframe[:, grouping])))
    fig, gdf, plot_layout = _setup_grouped_figure(
        dataframe,
        grouping;
        x_fig_size=fig_x_size,
        y_fig_size=fig_y_size,
        fontsize=fontsize
    )

    labels = label_lines.(first(df[:, grouping]) for df in gdf; l_length=17)
    xsize, ysize = _axis_size(gdf, fig_x_size, fig_y_size, n_col; y_gap=1.2)

    for (xi, groupdf) in enumerate(gdf)
        plot_layout_xi = plot_layout[xi]
        groupdf = sort(groupdf, cluster_col)

        # Ensure that relative cover timeseries match the cluster allocations from groupdf
        group_timeseries = timeseries_array[length_t, timeseries_array.locations.∈[groupdf.UNIQUE_ID]]
        group_timeseries = group_timeseries[:, indexin(groupdf.UNIQUE_ID, String.(group_timeseries.locations))]
        timesteps = group_timeseries.timesteps

        clusters = Int64.(groupdf[:, cluster_col])


        ADRIA.viz.clustered_scenarios!(
            fig[plot_layout_xi...],
            group_timeseries,
            clusters;
            opts=Dict{Symbol,Any}(:legend => false),
            axis_opts=Dict(
                :title => labels[xi],
                :xticks => (
                    first(length_t):10:last(length_t),
                    string.(collect((first(timesteps):10:last(timesteps))))
                ),
                :spinewidth => 0.5,
                :ylabelpadding => 2,
                :xlabelpadding => 2,
                :xticksize => 2,
                :yticksize => 2
            )
        )
    end

    feat = timeseries_array.properties[:metric_feature]
    unit_text = timeseries_array.properties[:metric_unit]
    Label(
        fig[:, -1], "$feat [$unit_text]",
        rotation=π / 2,  # Rotate 90 degrees
        tellheight=false,
        tellwidth=true
    )

    if grouping != :gbr
        bottom_axes = [last(plot_layout[last.(plot_layout).==x]) for x in 1:n_col]
        not_bottom_axes = findall(plot_layout .∉ [bottom_axes])
        not_bottom_axes = filter(x -> x isa Axis, fig.content)[not_bottom_axes]
        # second_row = findfirst(last.(plot_layout) .== n_col) + 1
        # second_last_row = findfirst(x -> first(x) == maximum(first.(plot_layout)), plot_layout) - 1
        # middle_axes = filter(x -> x isa Axis, fig.content)[second_row:second_last_row]
        map(x -> hidexdecorations!(x; grid=false, ticks=false), not_bottom_axes)
    end

    axes = filter(x -> x isa Axis, fig.content)
    map(x -> hideydecorations!(x; ticks=false, ticklabels=false, grid=false), axes)

    if grouping == :bioregion
        map(x -> colsize!(fig.layout, x, Relative(1 / n_col)), 1:n_col)
        map(x -> rowsize!(fig.layout, x, Relative(1 / maximum(first.(plot_layout)))), 1:n_col)
    end

    last_figure = last(plot_layout)
    if last(last_figure) > n_col
        legend_position = (last_figure[1] + 1, 1:n_col)
        rowsize!(fig.layout, first(legend_position), Relative(0.03))
    end

    if grouping == :gbr
        colsize!(fig.layout, 2, Relative(0.1))
    end
    Makie.trim!(fig.layout)

    return fig
end

"""
    grouped_cluster_ridgeline_plot(
        dataframe::DataFrame,
        cluster_col::Union{Symbol, String},
        grouping::Union{Symbol, String},
        variable::Union{Symbol, String};
        title::String="",
        xlabel::String="Cluster",
        ylabel::String="",
        fig_sizes::Dict=fig_sizes,
        fontsize::Float64=fontsize,
        datalimits::Tuple=(-Inf,Inf),
        overlap::Union{Float64, Int64}=1
    )

Plot clustered distribution of `variable` values for each group in `grouping`. Designed to
plot with 3 clusters in each group. Distributions are plotted using `Makie.violin!()``.
X-axis is values from distributions, Y-axis is levels from `grouping`.

# Arguments
- `dataframe` : Contains values for `cluster_col`, `grouping`, `variable` columns.
- `cluster_col` : Column containing numerical cluster assignments for each location.
- `grouping` : Column containing location grouping levels used to separate locations.
Unique levels of `grouping` are used as yticklabels.
- `variable` : Column containing values to plot into density-distributions.
- `title`
- `xlabel`
- `ylabel`
- `fig_sizes` : Dictionary containing violin_ width and height measurements for figure size.
- `fontsize`
- `datalimits` : Tuple containing minimum and maximum values for x axis.
- `overlap` : Amount of distance between y-axis group ridgelines.
"""
function grouped_cluster_ridgeline_plot(
    dataframe::DataFrame,
    cluster_col::Union{Symbol,String},
    grouping::Union{Symbol,String},
    variable::Union{Symbol,String};
    title::String="",
    xlabel::String="Cluster",
    ylabel::String="",
    fig_sizes::Dict=fig_sizes,
    fontsize::Float64=fontsize,
    datalimits::Tuple=(-Inf, Inf),
    overlap::Union{Float64,Int64}=1
)

    fig_x_size = fig_sizes["violin_width"]
    fig_y_size = fig_sizes["violin_height"]
    # Prepare plot
    fig, gdf, plot_layout = _setup_grouped_figure(
        dataframe,
        grouping;
        x_fig_size=fig_x_size,
        y_fig_size=fig_y_size,
        fontsize=fontsize,
        marker=PolyElement,
        multi_axis=false,
        rev=false
    )

    #labels = label_lines.([first(df[:, grouping]) for df in gdf]; l_length=17)
    labels = [first(df[:, grouping]) for df in gdf]
    colors = [:green, :orange, :blue]
    ax = Axis(
        fig[1, 1],
        title=title,
        ylabel=ylabel,
        xlabel=xlabel,
        limits=(extrema(dataframe[:, variable]), nothing)
    )

    for (i, groupdf) in enumerate(gdf)

        groupdf = sort(groupdf, cluster_col)

        clusters = Int64.(groupdf[:, cluster_col])
        cluster_order = unique(clusters)

        if variable == :so_to_si
            vlines!(1; color=(:gray, 0.5), linewidth=2)
        elseif variable == :log_so_to_si
            vlines!(0; color=(:gray, 0.5), linewidth=2)
        end

        y_offset = (i - 1) * overlap

        for (j, cluster) in enumerate(cluster_order)
            d = groupdf[groupdf[:, cluster_col].==cluster, variable]
            cluster_color = colors[j]

            violin!(
                fill(y_offset, length(d)),
                d,
                color=(cluster_color, 0.2),
                orientation=:horizontal,
                show_median=true,
                mediancolor=cluster_color,
                side=:right,
                datalimits=datalimits
            )
            # rainclouds!(
            #     fill(y_offset, length(d)),
            #     d;
            #     color=(cluster_colors[cluster], 0.6),
            #     markersize=5,
            #     jitter_width=0.27,
            #     side_nudge=0.25,
            #     plot_boxplots=false,
            #     clouds=nothing,
            #     side = :left,
            #     orientation=:horizontal
            # )
        end
    end

    # Adjust y-axis labels to show group names
    y_ticks = [(i - 1) * overlap for i in eachindex(labels)]
    ax.yticks = (y_ticks, labels)
    hidespines!(ax)

    legend_entries = []
    for (i, col) in enumerate([:green, :orange, :blue])
        LE = PolyElement(; color=col, marker=:circle)
        push!(legend_entries, [LE])
    end

    Legend(
        fig[1, 2],
        legend_entries,
        ["Low", "Medium", "High"],
        nbanks=1
    )

    return fig
end

grouping_full_names = Dict(
    :bioregion => "Bioregion",
    :management_area => "Management Area",
    :gbr => ""
)

"""
    cluster_analysis_plots(
        GCM,
        analysis_layers,
        rel_cover,
        dhw_ts,
        grouping,
        fig_out_dir;
        grouping_full_names=grouping_full_names,
        dpi=dpi
    )

Plot all analysis plots for paper.
"""
function cluster_analysis_plots(
    GCM,
    analysis_layers,
    rel_cover,
    dhw_ts,
    grouping,
    fig_out_dir;
    grouping_full_names=grouping_full_names,
    dpi=dpi
)
    grouping_fn = grouping_full_names[grouping]

    overlap = 0.8
    # Filter out groups that don't have 3 clusters due to earlier filtering.
    groups_too_few_clusters = grouping_counts(grouping, analysis_layers, "$(GCM)_$(grouping)_clusters", 3, 3)
    analysis_layers = analysis_layers[analysis_layers[:, grouping].∉[groups_too_few_clusters], :]
    rel_cover = rel_cover = rel_cover[:, rel_cover.locations.∈[analysis_layers.UNIQUE_ID]]
    dhw_ts = dhw_ts[:, dhw_ts.locations.∈[rel_cover.locations]]

    mean_dhw_violin = grouped_cluster_ridgeline_plot(
        analysis_layers,
        Symbol("$(GCM)_$(grouping)_clusters"),
        grouping, Symbol("$(GCM)_mean_dhw");
        xlabel="mean DHW", ylabel="$(grouping_fn)", overlap=overlap
    )
    save(
        joinpath(fig_out_dir, "$(grouping)", "mean_dhw_$(grouping)_violin.png"),
        mean_dhw_violin,
        px_per_unit=dpi
    )

    dhw_tol_violin = grouped_cluster_ridgeline_plot(
        analysis_layers,
        Symbol("$(GCM)_$(grouping)_clusters"),
        grouping, Symbol("$(GCM)_mean_DHW_tol");
        xlabel="mean reef DHW tolerance", ylabel="$(grouping_fn)", overlap=overlap
    )
    save(
        joinpath(fig_out_dir, "$(grouping)", "dhw_tolerance_$(grouping)_violin.png"),
        dhw_tol_violin,
        px_per_unit=dpi
    )

    so_to_si_violin = grouped_cluster_ridgeline_plot(
        analysis_layers,
        Symbol("$(GCM)_$(grouping)_clusters"),
        grouping, :log_so_to_si;
        xlabel="Log source to sink ratio", ylabel="$(grouping_fn)", overlap=overlap
    )
    save(
        joinpath(fig_out_dir, "$(grouping)", "so_to_si_$(grouping)_violin.png"),
        so_to_si_violin,
        px_per_unit=dpi
    )

    total_strength_violin = grouped_cluster_ridgeline_plot(
        analysis_layers,
        Symbol("$(GCM)_$(grouping)_clusters"),
        grouping, :log_total_strength;
        xlabel="Log total connectivity strength", ylabel="$(grouping_fn)", overlap=overlap
    )
    save(
        joinpath(fig_out_dir, "$(grouping)", "log_total_strength_$(grouping)_violin.png"),
        total_strength_violin,
        px_per_unit=dpi
    )

    initial_proportion_violin = grouped_cluster_ridgeline_plot(
        analysis_layers,
        Symbol("$(GCM)_$(grouping)_clusters"),
        grouping, Symbol("initial_proportion");
        xlabel="Initial proportion of habitable area occupied", ylabel="$(grouping_fn)", overlap=overlap
    )
    save(
        joinpath(fig_out_dir, "$(grouping)", "initial_proportion_$(grouping)_violin.png"),
        initial_proportion_violin,
        px_per_unit=dpi
    )

    dhw_cover_cor_violin = grouped_cluster_ridgeline_plot(
        analysis_layers,
        Symbol("$(GCM)_$(grouping)_clusters"),
        grouping, Symbol("$(GCM)_dhw_cover_cor");
        xlabel="Total coral cover - DHW correlation", ylabel="$(grouping_fn)", overlap=overlap
    )
    save(
        joinpath(fig_out_dir, "$(grouping)", "dhw_cover_cor_$(grouping)_violin.png"),
        dhw_cover_cor_violin,
        px_per_unit=dpi
    )

    weighted_incom_conn = grouped_cluster_ridgeline_plot(
        analysis_layers,
        Symbol("$(GCM)_$(grouping)_clusters"),
        grouping, Symbol("$(GCM)_weighted_incoming_conn_log");
        xlabel="Log total weighted incoming connectivity", ylabel="$(grouping_fn)", overlap=overlap
    )
    save(
        joinpath(fig_out_dir, "$(grouping)", "weighted_incoming_conn_$(grouping)_violin.png"),
        weighted_incom_conn,
        px_per_unit=dpi
    )

    analysis_layers_depth = analysis_layers[analysis_layers.depth_mean.!=7, :]
    groups_too_few_clusters_depth = grouping_counts(
        grouping,
        analysis_layers_depth,
        "$(GCM)_$(grouping)_clusters",
        3,
        5
    )
    analysis_layers_depth = analysis_layers_depth[
        analysis_layers_depth[:, grouping].∉[groups_too_few_clusters_depth], :
    ]

    depth_median_violin = grouped_cluster_ridgeline_plot(
        analysis_layers_depth,
        Symbol("$(GCM)_$(grouping)_clusters"),
        grouping, Symbol("depth_med");
        xlabel="Median Depth [m]", ylabel="$(grouping_fn)", overlap=overlap
    )
    save(
        joinpath(fig_out_dir, "$(grouping)", "depth_$(grouping)_violin.png"),
        depth_median_violin,
        px_per_unit=dpi
    )

    bioregion_grouped_timeseries = grouped_cluster_timeseries_plots(
        rel_cover,
        analysis_layers,
        "$(GCM)_$(grouping)_clusters",
        grouping,
        1:50
    )
    save(
        joinpath(fig_out_dir, "$(grouping)", "$(grouping)_cover_timeseries.png"),
        bioregion_grouped_timeseries,
        px_per_unit=dpi
    )

    bior_reefs_dhw_plot = grouped_cluster_timeseries_plots(
        dhw_ts,
        analysis_layers,
        "$(GCM)_$(grouping)_clusters",
        grouping,
        1:50
    )
    save(
        joinpath(fig_out_dir, "$(grouping)", "grouping_dhw_timeseries.png"),
        bior_reefs_dhw_plot,
        px_per_unit=dpi
    )

    return nothing
end

"""
    grouped_GCM_cluster_timeseries_plots(
        timeseries_array::YAXArray,
        dataframe::DataFrame,
        cluster_col::Union{String, Symbol},
        grouping::Union{String, Symbol},
        length_t::UnitRange;
        fig_sizes::Dict=fig_sizes,
        fontsize::Float64=fontsize
    )

Plots clustered timeseries plots for each GCM, separated into different axes for `grouping` levels.
Each GCM is plotted vertically, with different GCMs in the x direction and different grouping levels
in the y direction in figure.

# Arguments
- `timeseries_array` : Timeseries values in array of nrow(dataframe) * length(length_t) * length(GCM-levels) dimensions.
- `dataframe` : Contains values for `cluster_col` and `grouping`.
- `cluster_col` : Column containing numerical cluster assignment values.
- `grouping` : Column containing grouping levels to separate timeseries in y direction.
- `length_t` : Timestep indices to plot from `timeseries_array`.
- `fig_sizes` : Dictionary containing timeseries_ width and height variables for figure sizing.
- `fontsize`
"""
function grouped_GCM_cluster_timeseries_plots(
    timeseries_array::YAXArray,
    dataframe::DataFrame,
    cluster_col::Union{String,Symbol},
    grouping::Union{Vector{String},Vector{Symbol}},
    length_t::UnitRange;
    fig_sizes::Dict=fig_sizes,
    fontsize::Float64=fontsize
)
    # Adjust width and height just for the aspect ratio
    # Doesn't seem to affect the word doc output
    fig_x_size = fig_sizes["gcm_timeseries_width"] * 1.2
    fig_y_size = fig_sizes["gcm_timeseries_height"] * 0.95
    n_col = optimum_columns(size(unique(dataframe[:, grouping]), 1))

    # Calculate layout scheme and grouped dataframe, but ignore figure.
    # It sets a legend using predetermined row numbers but we want more flexibility
    _, gdf, plot_layout = _setup_grouped_figure(
        dataframe,
        grouping;
        x_fig_size=fig_x_size,
        y_fig_size=fig_y_size,
        fontsize=fontsize,
        order=[:management_area_average_latitude, :GCM]
    )

    mgmt_area = [first(df.management_area) for df in gdf]
    mgmt_area = unique(replace.(mgmt_area, " Management Area" => ""))
    gcm_header = unique([first(df.GCM) for df in gdf])

    # Create custom figure
    fig = Figure(size=(fig_x_size, fig_y_size), fontsize=fontsize)

    current_row = 0
    current_col = 0
    for (xi, groupdf) in enumerate(gdf)
        plot_layout_xi = plot_layout[xi]
        groupdf = sort(groupdf, cluster_col)

        # Ensure that relative cover timeseries match the cluster allocations from groupdf
        gcm = first(groupdf.GCM)
        group_timeseries = timeseries_array[GCM=(timeseries_array.GCM .== gcm)][:, :, 1]
        group_timeseries = group_timeseries[locations=(group_timeseries.locations .∈ [groupdf.UNIQUE_ID])]
        group_timeseries = group_timeseries[length_t, indexin(groupdf.UNIQUE_ID, String.(group_timeseries.locations))]

        groupdf = groupdf[groupdf.UNIQUE_ID.∈[group_timeseries.locations], :]
        clusters = Int64.(groupdf[:, cluster_col])

        ADRIA.viz.clustered_scenarios!(
            fig[plot_layout_xi...],
            group_timeseries,
            clusters;
            opts=Dict{Symbol,Any}(:legend => false),
            axis_opts=Dict(
                :xticks => (
                    [1, 10, 20, 30, 40, 50],
                    string.([2025, 2035, 2045, 2055, 2065, 2075])
                ),
                :spinewidth => 0.5,
                :ylabelpadding => 2,
                :xlabelpadding => 2,
                :xticksize => 2,
                :yticksize => 2
            )
        )

        # Add shared row labels (management area names)
        if plot_layout_xi[1] != current_row
            current_row = plot_layout_xi[1]
            Label(
                fig[plot_layout_xi[1], 0],
                mgmt_area[plot_layout_xi[1]],
                tellheight=false,
                tellwidth=true,
                fontsize=6
            )
        end

        # Add shared column labels (GCM names)
        if current_col != plot_layout_xi[2]
            Label(
                fig[0, plot_layout_xi[2]],
                gcm_header[plot_layout_xi[2]],
                tellheight=true,
                tellwidth=false,
                fontsize=6
            )
            current_col = plot_layout_xi[2]
        end
    end

    # Tweak layout defaults
    # fig.layout.colgap = 5  # Default ~10, reduce to make plots tighter horizontally
    # fig.layout.rowgap = 5  # Same for vertical spacing
    # fig.layout.outerpadding = 5  # Reduce space around the entire figure

    linkaxes!(filter(x -> x isa Axis, fig.content)...)

    # Add figure y-axis label
    feat = timeseries_array.properties[:metric_feature]
    unit_text = timeseries_array.properties[:metric_unit]
    Label(
        fig[:, -1], "$feat [$unit_text]",
        rotation=π / 2,  # Rotate 90 degrees
        tellheight=false,
        tellwidth=true
    )

    n_row, n_col = size(fig.layout)

    # Add figure x-axis label
    Label(
        fig[n_row+1, 3], "Year",
        tellheight=false,
        tellwidth=false
    )

    # Add legend
    legend_entries = []
    for col in [:green, :orange, :blue]
        push!(legend_entries, LineElement(; color=col, marker=:circle))
    end
    Legend(
        fig[n_row+2, 3],
        legend_entries,
        ["Low", "Medium", "High"],
        nbanks=3,
        tellheight=true,
        tellwidth=false,
        padding=(1.0, 2.0, 2.0, 2.0),             # shrink padding inside legend box
        labelsize=fontsize,    # smaller font
        framevisible=false,        # optional: remove box
        markerlabelgap=3,          # reduce space between marker and label
        rowgap=0,                  # reduce vertical spacing between items
        colgap=4,                  # reduce horizontal spacing
        patchsize=(5, 5)
    )

    second_last_row = findfirst(x -> first(x) == maximum(first.(plot_layout)), plot_layout) - 1
    all_axes = filter(x -> x isa Axis, fig.content)
    middle_axes = all_axes[1:second_last_row]
    last_axes_row = all_axes[second_last_row+1:end]

    # Hide x/y axis details
    map(x -> hideydecorations!(x, ticks=false, ticklabels=false, grid=false), all_axes)
    map(x -> hidexdecorations!(x; grid=false, ticks=false), middle_axes)
    map(x -> hidexdecorations!(x; grid=false, ticks=false, ticklabels=false), last_axes_row)

    # Tweak height/width of rows/columns
    n_row, n_col = size(fig.layout)
    map(x -> colsize!(fig.layout, x, Auto(1.1)), -1:n_col-2)
    map(x -> rowsize!(fig.layout, x, Auto(1.1)), 0:n_row-1)
    rowsize!(fig.layout, n_row - 2, Auto(0.05))  # Year label
    rowsize!(fig.layout, n_row - 1, Auto(0.25))  # Legend

    # Clear out empty space
    Makie.trim!(fig.layout)

    return fig
end

"""
    carbonate_budget_variable_scatter(
        long_df::DataFrame,
        color_variable_col::Union{String, Symbol},
        year_col::Union{String, Symbol},
        carbonate_budget_col::Union{String, Symbol},
        year_variable_correlation::Dict;
        xlabel::String="Carbonate budget threshold [%] (correlation)",
        ylabel::String="Years above carbonate budget threshold",
        color_label::String="",
        fig_sizes::Dict=fig_sizes,
        fontsize::Float64=fontsize,
        alpha::Union{Float64, Int64}=0.5
    )

Create scatter plot displaying years reefs exceed carbonate budget thresholds (y axis), for
each unique carbonate budget threshold level in `carbonate_budget_col` (x axis). Scatter
points are coloured by `color_variable_col`.

# Arguments
- `long_df` : Longform DataFrame containing `color_variable_col` and `year_col` values under
each `carbonate_budget_col` threshold level.
- `color_variable_col` : Column containing Float64 values to colour the scatter points for each location.
- `year_col` : Column containing values for the number of years a reef maintains a positive carbonate budget.
- `carbonate_budget_col` : Column containing the carbonate budget threshold values used in analyses (Β).
- `year_variable_correlation` : Dictionary containing correlation values (or nothing values) for each unique value in `carbonate_budget_col`.
- `xlabel`
- `ylabel`
- `color_label`
- `fig_sizes` : Dictionary containing carb_ width and height variables for figure sizing.
- `fontsize`
- `alpha` : Transparency applied to scatter points for each location.
"""
function carbonate_budget_variable_scatter(
    long_df::DataFrame,
    color_variable_col::Union{String,Symbol},
    year_col::Union{String,Symbol},
    carbonate_budget_col::Union{String,Symbol},
    year_variable_correlation::Dict;
    xlabel::String="Carbonate budget threshold [%] (correlation)",
    ylabel::String="Years above carbonate budget threshold",
    color_label::String="",
    fig_sizes::Dict=fig_sizes,
    fontsize::Float64=fontsize,
    alpha_c::Union{Float64,Int64}=0.5
)
    x_fig_size = fig_sizes["carb_width"]
    y_fig_size = fig_sizes["carb_height"]

    labels = unique(long_df[:, carbonate_budget_col])
    labels = ["$(lab)\n($(round(year_variable_correlation[lab], digits=2)))" for lab in labels]

    # Get the default colormap (as a gradient with 256 colors)
    base_cmap = cgrad(:viridis, 256; rev=true)
    # Convert all colors to RGBA with alpha = 0.4 (adjust as needed)
    transparent_colors = [RGBA(c.r, c.g, c.b, alpha_c) for c in base_cmap.colors]
    # Wrap into a ColorScheme object
    transparent_cmap = ColorScheme(transparent_colors)

    theta_means = combine(groupby(long_df, carbonate_budget_col)) do sdf
        (;
            x=first(sdf[:, carbonate_budget_col]),
            y=mean(sdf[:, year_col])
        )
    end

    fig = Figure(size=(x_fig_size, y_fig_size), fontsize=fontsize)
    ax = Axis(
        fig[1, 1],
        ylabel=ylabel,
        xlabel=xlabel,
        xticks=(1:1:11, labels)
    )
    rain = rainclouds!(
        long_df[:, carbonate_budget_col],
        long_df[:, year_col];
        color=long_df[:, color_variable_col],
        plot_boxplots=false,
        clouds=nothing,
        show_median=false,
        markersize=6,
        jitter_width=1
    )
    rain.plots[1].colormap = transparent_cmap
    rainclouds!(
        ax,
        theta_means.x,
        theta_means.y;
        color=(:red, 0.8),
        markersize=8,
        show_median=false,
        clouds=nothing,
        plot_boxplots=false
    )

    Colorbar(fig[1, 2], rain, label=color_label)

    return fig
end

"""
    gcm_cluster_assignment_heatmap(
        dataframe::DataFrame,
        grouping::Union{String, Symbol},
        cluster_cols::Union{String, Symbol};
        fig_sizes::Dict=fig_sizes,
        title::String="",
        ylabel::String="",
        xlabel::String="",
        GCMs::Vector{String}=GCMs
    )

Plot a heatmap showing cluster assignments over the different GCMs for each location. Heatmaps
are grouped by `grouping` if there are more than 1 unique levels of `grouping`. Y axes show number
of locations in groups, x axes show the GCMs. Heatmaps are coloured with the cluster levels (low, medium, high).


# Arguments
- `dataframe` : Dataframe containing all of the required `cluster_cols` and `grouping` values.
- `grouping` : Column containing grouping to split heatmaps.
- `cluster_cols` : Vector of column names containing numerical cluster assignemnts for each GCM.
GCM cluster columns should be in the same order as labels in GCMs vector.
- `fig_sizes` : Dictionary containing timeseries_ width and height variables for sizing figure.
- `title` : Main plot title.
- `ylabel`
- `xlabel`
- `GCMs` : Vector containing GCM string labels (In the same order as GCM cluster columns in `cluster_cols`!).
"""
function gcm_cluster_assignment_heatmap(
    dataframe::DataFrame,
    grouping::Union{String,Symbol},
    cluster_cols::Union{Vector{String},Vector{Symbol}};
    fig_sizes::Dict=fig_sizes,
    title::String="",
    ylabel::String="",
    xlabel::String="",
    GCMs::Vector{String}=GCMs
)
    fig_x_size = fig_sizes["cluster_hm_width"]
    fig_y_size = fig_sizes["cluster_hm_height"]
    n_col = optimum_columns(length(unique(dataframe[:, grouping])))
    # Prepare plot
    fig, gdf, plot_layout = _setup_grouped_figure(
        dataframe,
        grouping;
        x_fig_size=fig_x_size,
        y_fig_size=fig_y_size,
        fontsize=9,
        marker=PolyElement,
        multi_axis=true
    )

    #labels = label_lines.([first(df[:, grouping]) for df in gdf]; l_length=17)
    # labels = [first(df[:, :bioregion]) for df in gdf]
    colors = [:green, :orange, :blue]
    labels = label_lines.(first(df[:, grouping]) for df in gdf; l_length=9)
    # xsize, ysize = _axis_size(gdf, fig_x_size, fig_y_size, n_col; y_gap=0.8)
    # xsize, ysize = (50, 50)

    for (xi, groupdf) in enumerate(gdf)
        plot_layout_xi = plot_layout[xi]
        ax = _setup_grouped_axes(
            fig,
            plot_layout_xi,
            ([1, 2, 3, 4, 5], GCMs);
            ylabel=ylabel,
            xlabel=xlabel,
            title=title,
            xsize=nothing,
            ysize=nothing,
            background_color=:white,
            xticklabelrotation=45.0,
            spinewidth=0.5
        )

        cluster_matrix = Matrix(groupdf[:, [cluster_cols...]])
        hm = heatmap!(ax, cluster_matrix', colormap=colors)
        # groupdf = sort(groupdf, cluster_col)

        # Ensure that relative cover timeseries match the cluster allocations from groupdf
        Label(
            fig[plot_layout_xi..., Top()],
            labels[xi],
            fontsize=8pt,
            font=:bold,
            padding=(0, 1, 2, 0),
            halign=:center
        )
    end

    if grouping != :gbr
        # second_row = findfirst(last.(plot_layout) .== n_col) + 1
        second_last_row = findfirst(x -> first(x) == maximum(first.(plot_layout)), plot_layout) - 1
        middle_axes = filter(x -> x isa Axis, fig.content)[1:second_last_row]
        # axes_after_1 = filter(x -> x isa Axis, fig.content)[2:end]

        map(x -> hidexdecorations!(x; grid=false, ticks=false), middle_axes)
        map(x -> hideydecorations!(x; ticks=false, ticklabels=false, grid=false), filter(x -> x isa Axis, fig.content))
    end

    if grouping == :bioregion
        n_row = maximum(first.(plot_layout))
        map(x -> rowsize!(fig.layout, x, Relative(1 / n_row)), 1:n_row)
        map(x -> colsize!(fig.layout, x, Relative(1 / n_col)), 1:n_col)
        rowgap!(fig.layout, 2)
    end

    return fig
end


function map_gbr_reefs(reef_df, color_col::Symbol, colormap, color_legend_label; management_region_fn="../data/GBRMPA_Management_Areas.gpkg", mainland_fn="../data/GBRMPA_Reef_Features.gpkg", fig_sizes=fig_sizes, fontsize=fontsize)
    regions = GDF.read(management_region_fn)
    regions.region_name = replace.(regions.AREA_DESCR, [" Management Area" => ""])
    qld = GDF.read(mainland_fn)
    qld = qld[qld.FEAT_NAME.=="Mainland", :SHAPE]

    map_width = fig_sizes["map_width"]
    map_height = fig_sizes["map_height"]
    region_col = tuple.([:blue, :green, :orange, :red], fill(0.2, 4)) # Manually set alpha value to 0.3

    if color_col == :bioregion
        ordered_reefs = sort(reef_df, :bioregion_average_latitude; rev=true)
    else
        ordered_reefs = reef_df
    end

    reef_colors = colormap
    reef_colors = tuple.(reef_colors, fill(0.5, length(reef_colors)))

    bioregion_colors_labels = Dict(
        color_col => unique(reef_df[:, color_col]),
        :ref => 1:length(unique(reef_df[:, color_col])),
        :color => reef_colors,
        :label => label_lines.((unique(reef_df[:, color_col])); l_length=17)
    )
    bioregion_colors_labels = DataFrame(bioregion_colors_labels)

    ordered_reefs = leftjoin(ordered_reefs, bioregion_colors_labels; on=color_col, order=:left)
    centroids = GO.centroid.(ordered_reefs.geometry)

    # GeoMakie currently has a bug where x/y labels never display.
    # We adjust map size to roughly align with the correct projection.
    bgcol = :gray90
    fig = Figure(size=(map_width + 75, map_height), fontsize=fontsize, backgroundcolor=bgcol)
    ax = Axis(
        fig[1, 1],
        width=map_width - 275,
        height=map_height - 120,
        xgridvisible=false,
        ygridvisible=false,
        limits=((142.5, 154.1), (-25, -10)),
        backgroundcolor=bgcol
    )
    poly!(qld; color=:darkgray)
    poly!(regions.SHAPE, color=region_col)

    scatter!(centroids; color=ordered_reefs.color, markersize=4)
    ax.ylabel = "Latitude"
    ax.xlabel = "Longitude"

    # Note key cities/towns for reference
    scatter!((145.754120, -16.925491); color=:black)
    text!((145.1, -16.925491); text="Cairns", align=(:center, :top))

    scatter!((146.8057, -19.2664); color=:black)
    text!((146.0, -19.2664); text="Townsville", align=(:center, :top))

    scatter!((149.182147, -21.142496); color=:black)
    text!((148.432147, -21.142496); text="Mackay", align=(:center, :top))

    scatter!((150.733333, -23.133333); color=:black)
    text!((150.0, -23.133333); text="Yeppoon", align=(:center, :top))

    # lines!([(146.25, -20.5), (147.0559, -19.2697)]; color=:black)
    # text!((146.25, -20.5); text="AIMS - Cape Cleveland", align=(:center, :top))

    Legend(
        fig[2, 1],
        [PolyElement(color=col, alpha=0.2) for col in region_col],
        regions.region_name,
        "Management regions",
        nbanks=2,
        patchsize=(10, 10),
        colgap=8,
        backgroundcolor=bgcol
    )
    Legend(
        fig[1, 2],
        [MarkerElement(color=bioregion_color, marker=:circle) for bioregion_color in unique(ordered_reefs.color)],
        titlecase.(unique(ordered_reefs.label)),
        color_legend_label,
        nbanks=2,
        colgap=6,
        rowgap=1,
        backgroundcolor=bgcol
    )

    return fig
end

function carbonate_budget_subplot!(
    fig, 
    long_df, 
    grid_layout_pos, 
    color_variable_col, 
    carbonate_budget_col, 
    year_col, 
    labels, 
    color_label, 
    theta_year_means,
    colormap
)


    if grid_layout_pos[2] == 1
        bar_pos = (1,0)
        flip_axis = false
        yaxis_pos = :right
    else
        bar_pos = (1,2)
        flip_axis = true
        yaxis_pos = :left
    end

    grid_x = GridLayout(fig[grid_layout_pos...])
    ax = Axis(
        grid_x[1, 1],
        xlabel="Carbonate budget threshold [%] (correlation)",
        ylabel="Years above carbonate budget threshold",
        xticks=(1:1:11, labels),
        xticklabelsize=fontsize-2,
        yaxisposition = yaxis_pos
    )
    rain = rainclouds!(
        ax,
        long_df[:, carbonate_budget_col],
        long_df[:, year_col];
        color=long_df[:, color_variable_col],
        plot_boxplots=false,
        clouds=nothing,
        show_median=false,
        markersize=6,
        jitter_width=0.7
    )
    rain.plots[1].colormap = colormap
    rainclouds!(
        ax,
        theta_year_means.x,
        theta_year_means.y;
        color = (:red, 0.8),
        markersize=8,
        show_median=false,
        clouds=nothing,
        plot_boxplots=false
    )

    Colorbar(grid_x[bar_pos...], rain, label=color_label, flip_vertical_label=false, size=6, flipaxis=flip_axis, spinewidth=0.0)

    return fig
end

"""
    carbonate_budget_variable_scatter(
        long_df::DataFrame,
        color_variable_col::Union{String, Symbol},
        year_col::Union{String, Symbol},
        carbonate_budget_col::Union{String, Symbol},
        year_variable_correlation::Dict;
        xlabel::String="Carbonate budget threshold [%] (correlation)",
        ylabel::String="Years above carbonate budget threshold",
        color_label::String="",
        fig_sizes::Dict=fig_sizes,
        fontsize::Float64=fontsize,
        alpha::Union{Float64, Int64}=0.5
    )

Create scatter plot displaying years reefs exceed carbonate budget thresholds (y axis), for
each unique carbonate budget threshold level in `carbonate_budget_col` (x axis). Scatter
points are coloured by `color_variable_col`.

# Arguments
- `long_df` : Longform DataFrame containing `color_variable_col` and `year_col` values under
each `carbonate_budget_col` threshold level.
- `color_variable_col` : Column containing Float64 values to colour the scatter points for each location.
- `year_col` : Column containing values for the number of years a reef maintains a positive carbonate budget.
- `carbonate_budget_col` : Column containing the carbonate budget threshold values used in analyses (Β).
- `year_variable_correlation` : Dictionary containing correlation values (or nothing values) for each unique value in `carbonate_budget_col`.
- `xlabel`
- `ylabel`
- `color_label`
- `fig_sizes` : Dictionary containing carb_ width and height variables for figure sizing.
- `fontsize`
- `alpha` : Transparency applied to scatter points for each location.
"""
function carbonate_budget_variable_scatter(
    long_df::DataFrame,
    year_col::Union{String,Symbol},
    carbonate_budget_col::Union{String,Symbol},
    year_depth_correlation::Dict,
    year_conn_correlation::Dict;
    depth_var_col=:depth_med,
    conn_var_col=:log_total_strength,
    depth_color_label="Reef median depth [m]",
    conn_color_label="Log total connectivity strength",
    fig_sizes::Dict=fig_sizes,
    fontsize::Float64=fontsize,
    alpha_c::Union{Float64,Int64}=0.5
)
    x_fig_size = fig_sizes["carb_width"]
    y_fig_size = fig_sizes["carb_height"]

    labels = unique(long_df[:, carbonate_budget_col])
    depth_labels = ["$(lab)\n($(round(year_depth_correlation[lab], digits=2)))" for lab in labels]
    conn_labels = ["$(lab)\n($(round(year_conn_correlation[lab], digits=2)))" for lab in labels]

    # Get the default colormap (as a gradient with 256 colors)
    base_cmap = cgrad(:viridis, 256; rev=true)
    # Convert all colors to RGBA with alpha = 0.4 (adjust as needed)
    transparent_colors = [RGBA(c.r, c.g, c.b, alpha_c) for c in base_cmap.colors]
    # Wrap into a ColorScheme object
    transparent_cmap = ColorScheme(transparent_colors)

    theta_means = combine(groupby(long_df, carbonate_budget_col)) do sdf
        (; 
            x = first(sdf[:, carbonate_budget_col]),
            y = mean(sdf[:, year_col])
        )
    end

    fig = Figure(size=(x_fig_size, y_fig_size), fontsize=fontsize)
    carbonate_budget_subplot!(fig, long_df, (1,1), depth_var_col, carbonate_budget_col, year_col, depth_labels, depth_color_label, theta_means, transparent_cmap)
    carbonate_budget_subplot!(fig, long_df, (1,2), conn_var_col, carbonate_budget_col, year_col, conn_labels, conn_color_label, theta_means, transparent_cmap)

    grid_layouts = filter(x -> x isa Axis, fig.content)
    hideydecorations!(grid_layouts[1]; grid=false, ticks=false)
    hideydecorations!(grid_layouts[2]; grid=false, ticks=false, ticklabels=false)
    map(x -> hidexdecorations!(x; grid=false, ticks=false, ticklabels=false), grid_layouts)
    Makie.trim!(fig.layout)
    colgap!(fig.layout, 3.5)

    Label(fig[0, 1:2], "Number of years exceeding positive carbonate budget threshold", tellwidth=false)
    Label(fig[2, 1:2, Top()], "Carbonate budget threshold [%] (correlation)", tellwidth=false)
    rowsize!(fig.layout, 2, Relative(0.025))

    return fig
end