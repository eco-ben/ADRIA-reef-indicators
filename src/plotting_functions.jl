"""
File includes the plotting functions for ADRIA-reef-indicators analysis.
"""

using Printf
using Colors
using Random
using TOML


config = TOML.parsefile("config.toml")

# Extract figure sizes and convert to pixel measurement
fig_sizes = config["fig_formats_cm"]
cm = 37.7952755906 # Size of 1cm in pixels
map!(x -> x*cm, values(fig_sizes))

# Extract figure fontsize and convert to pixel measurement
pt = 1.33 # size of 1 pt in pixels
misc_sizes = config["fig_misc"]
fontsize = misc_sizes["fontsize"] * pt

inch = 96 # size of 1 inch in pixels
dpi = misc_sizes["dpi"] / inch

function colorscheme_alpha(cscheme::ColorScheme, alpha = 0.5)
    ncolors = length(cscheme)

    return ColorScheme([RGBA(get(cscheme, k), alpha) for k in range(0, 1, length=ncolors)])
end

"""
    bellwether_reef_numbers_plot(df, variable_col, value_col, xlabel="lags", ylabel="number of bellwether reefs")

Intended to plot the numbers of bellwether reefs seen in each analysis level on a scatter plot.

# Arguments
- `df` : A long format dataframe with column lags containing Integers for each relevant lag,
a variable column with the analysis level names for each lag, a value column with the numbers of bellwether reefs.
- `variable_col` : name of column containing the analysis level strings for colouring points.
- `value_col` : name of column containing the numbers of bellwether reefs for each lag/analysis level.
"""
function bellwether_reef_numbers_plot(df, variable_col, value_col, xlabel="lags", ylabel="number of bellwether reefs")
    f = Figure()
    ax = Axis(f[1,1]; xlabel=xlabel, ylabel=ylabel)

    alph_palette = colorscheme_alpha(ColorSchemes.tableau_20, 0.7);
    palette = ColorSchemes.tableau_20.colors;

    color_indices = groupindices(DataFrames.groupby(df, variable_col))
    color_names = unique(DataFrame(indices=color_indices, names=df[:, variable_col]))

    # Create the unique legend entries for each level of variable_col
    unique_names = color_names.names
    legend_entries = []
    for name in eachrow(color_names)
        col = palette[name.indices]
        LE = PolyElement(; color=col)
        push!(legend_entries, [LE])
    end

    points = scatter!(df.lags, df[:, value_col], color=alph_palette[color_indices], markersize=10)

    Legend(f[2, 1], legend_entries, unique_names, nbanks=3, tellheight=true,
    tellwidth=false, orientation=:horizontal, labelsize=10)

    display(f)
end

function skipmissing_median(x)
    new_median = Vector{Union{Missing, Float64}}(missing, size(x, 1))
    for (i, row) in enumerate(eachrow(x))
        if any(ismissing.(row))
            continue
        end
        new_median[i] = median(row)
    end

    return new_median
end

function label_lines(label; l_length=25)
    if length(label) > l_length
        label_new = replace(label, Regex("(.{$(l_length)} )") => s"\1\n")

        return label_new
    end

    return label
end

function _axis_size(gdf, x_fig_size, y_fig_size, n_col; x_gap = 1.5, y_gap = 1.2)
    xsize = x_fig_size / (n_col*x_gap)
    n_fig_row = first(fldmod1(length(gdf), n_col))
    ysize = y_fig_size / (n_fig_row*y_gap)

    return xsize, ysize
end

function _extract_name_and_correlation(df,bellwether_reefs_col, correlation_col, grouping)
    bioregion = first(df[:, grouping])
    bellwether_cor = round(mean(df[df[:, bellwether_reefs_col] .== "bellwether", correlation_col]), sigdigits=2)
    non_bellwether_cor = round(mean(df[df[:, bellwether_reefs_col] .== "non-bellwether", correlation_col]), sigdigits=2)

    return "$(bioregion) ($(@sprintf("%.2f",bellwether_cor))) ($(@sprintf("%.2f",non_bellwether_cor)))"
end

function _setup_grouped_figure(dataframe, bellwether_reefs_col, grouping; x_fig_size=2130, y_fig_size=1500)
    dataframe = sort(dataframe, [bellwether_reefs_col, :management_area]; rev=true)
    categories = categorical(dataframe[:, bellwether_reefs_col])

    gdf = DataFrames.groupby(dataframe, grouping)

    fig = Figure(backgroundcolor = RGBf(0.98, 0.98, 0.98), size = (x_fig_size, y_fig_size))

    colors = [Makie.wong_colors(); Makie.wong_colors()];
    legend_entries = []
    for (i, col) in enumerate(unique(colors[indexin(categories, unique(categories))]))
        LE = MarkerElement(; color=col, marker=:circle)
        push!(legend_entries, [LE])
    end

    n_col = optimum_columns(length(unique(dataframe[:, grouping])))
    plot_layout = figure_layouts(length(gdf), n_col)
    last_figure = last(plot_layout)
    if last(last_figure) < n_col
        legend_position = (last_figure[1], last_figure[2]+1:n_col)
    else
        legend_position = (last_figure[1] + 1, 1:n_col)
    end

    Legend(
        fig[legend_position...],
        legend_entries,
        unique(categories),
        nbanks=1
    )
    return fig, gdf, plot_layout, colors, categories
end

function _setup_grouped_figure(dataframe, grouping; x_fig_size=2130, y_fig_size=1500, marker=LineElement, order=:bioregion_average_latitude, multi_axis=true, fontsize=11pt)
    dataframe = sort(dataframe, order; rev=true)
    
    gdf = DataFrames.groupby(dataframe, grouping)

    fig = Figure(size = (x_fig_size, y_fig_size), fontsize = fontsize)
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
            nbanks=n_banks
        )
    end

    return fig, gdf, plot_layout
end

function _setup_grouped_axes(fig, plot_layout_xi, xticks; ylabel="", xlabel="", title="", xsize=220, ysize=150, background_color=:white)

    ax = Axis(
        fig[plot_layout_xi...];
        backgroundcolor=background_color,
        #xticklabelsize=16,
        #yticklabelsize=16,
        xticks = xticks,
        xticksvisible=false,
        title=title,
        width=xsize,
        height=ysize
    )
    return ax
end

function timeseries_xticks(length_t, years)
    length_range = first(length_t):10:last(length_t)
    years_length = collect(years)[length_range]

    return (length_range, string.(years_length))
end

function bootstrap_mean_ts(timeseries)
    ts_mean = Vector{Union{Missing, Float64}}(missing, size(timeseries, 1))
    ts_lb = Vector{Union{Missing, Float64}}(missing, size(timeseries, 1))
    ts_ub = Vector{Union{Missing, Float64}}(missing, size(timeseries, 1))

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
    elseif n_bioregions <30
        return 5
    end

    return 6
end

function GCM_label(x,y,GCM)
    return text!(x, y; text = GCM, align = (:right, :center))
end

function ecs_plot(ecs_values, low_conf_range, high_conf_range, GCM_labels; fig_sizes=fig_sizes, fontsize=fontsize)
    fig_x_size = fig_sizes["ecs_width"]
    fig_y_size = fig_sizes["ecs_height"]

    low_min, low_max = extrema(low_conf_range)
    high_min, high_max = extrema(high_conf_range)

    fig = Figure(size = (fig_x_size, fig_y_size), fontsize = fontsize)
    ax = Axis(
        fig[1,1],
        ylabel = "Equilibrium Climate Sensitivity (°C)",
        height = 10cm,
        width = 8cm
    )
    hidexdecorations!(ax)
    poly!(
        ax, 
        [(0.98, high_min), (1.02, high_min), (1.02, high_max), (0.98, high_max)];
        color = (:red, 0.2), linestyle = :dash, strokewidth = 1.5, label = "likely"
    )
    poly!(
        ax, 
        [(0.98, low_min), (1.02, low_min), (1.02, low_max), (0.98, low_max)];
        color = (:red, 0.6), strokewidth = 1.5, label = "very likely"
    )
    scatter!(ax, ones(length(ecs_values)), ecs_values, markersize = 15, color=:black)
    GCM_label.(fill(1.06, length(ecs_values)), ecs_values, GCM_labels)

    fig[2,1] = Legend(fig, ax, "ECS assessed range", framevisible = false, nbanks = 2)

    display(fig)

    return fig    
end

function grouped_cluster_timeseries_plots(
    timeseries_array,
    dataframe,
    cluster_col,
    grouping,
    length_t;
    fig_sizes=fig_sizes,
    fontsize=fontsize,
    ytitle=""
)
    fig_x_size = fig_sizes["timeseries_width"]
    fig_y_size = fig_sizes["timeseries_height"]
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
        group_timeseries = timeseries_array[length_t, timeseries_array.locations .∈ [groupdf.UNIQUE_ID]]
        group_timeseries = group_timeseries[:, indexin(groupdf.UNIQUE_ID, String.(group_timeseries.locations))]
        timesteps = group_timeseries.timesteps

        clusters = Int64.(groupdf[:, cluster_col])

        
        ADRIA.viz.clustered_scenarios!(
            fig[plot_layout_xi...],
            group_timeseries,
            clusters;
            opts = Dict{Symbol, Any}(:legend => false),
            axis_opts = Dict(
                :title => labels[xi], 
                :xticks => (
                    first(length_t):10:last(length_t), 
                    string.(collect((first(timesteps):10:last(timesteps))))
                )
                # :height => ysize,
                # :width => xsize
            )
        )
    end

    axes_after_1 = filter(x -> x isa Axis, fig.content)[2:end]
    map(x -> hidexdecorations!(x; grid = false, ticks = false), axes_after_1)
    map(x -> hideydecorations!(x, ticks = false, ticklabels = false, grid = false, label = false), axes_after_1)

    if grouping == :bioregion
        rowgap!(fig.layout, 5)
        colgap!(fig.layout, 5)
        map(x -> colsize!(fig.layout, x, Relative(1 / n_col)), 1:n_col)
    end

    # resize_to_layout!(fig)
    display(fig)

    return fig
end

function grouped_cluster_violin_plots(
    dataframe,
    cluster_col,
    grouping,
    variable;
    title="",
    xlabel="Cluster",
    ylabel="",
    fig_sizes=fig_sizes,
    fontsize=fontsize,
    datalimits=(-Inf,Inf)
)
    fig_x_size = fig_sizes["violin_width"]
    fig_y_size = fig_sizes["violin_height"]
    n_col = optimum_columns(length(unique(dataframe[:, grouping])))
    fig, gdf, plot_layout = _setup_grouped_figure(
        dataframe,
        grouping;
        x_fig_size=fig_x_size,
        y_fig_size=fig_y_size,
        fontsize=fontsize,
        marker=PolyElement
    )

    labels = label_lines.(first(df[:, grouping]) for df in gdf; l_length = 10)
    colors = [:green, :orange, :blue]
    xsize, ysize = _axis_size(gdf, x_fig_size, y_fig_size, n_col; y_gap=1.2, x_gap=1.2)
    xticks = (1:3, ["Low", "Medium", "High"])

    for (xi, groupdf) in enumerate(gdf)
        #println("$(xi)")
        plot_layout_xi = plot_layout[xi]
        groupdf = sort(groupdf, cluster_col)

        clusters = Int64.(groupdf[:, cluster_col])     

        ax = _setup_grouped_axes(
            fig,
            plot_layout_xi,
            xticks;
            ylabel=ylabel,
            xlabel=xlabel,
            title=title,
            xsize=xsize,
            ysize=ysize,
            background_color=:white
        )

        f = violin!(
            ax,
            clusters,
            groupdf[:, variable];
            color=colors[indexin(clusters, unique(clusters))],
            show_median=true,
            datalimits=datalimits
        )
        f = rainclouds!(
            ax,
            clusters,
            groupdf[:, variable];
            color=:gray,
            markersize=5,
            jitter_width=0.27,
            side_nudge=0.001,
            plot_boxplots=false,
            clouds=nothing
        )

        if variable == :so_to_si
            hlines!(1; color=(:gray, 0.5), linewidth=4)
        elseif variable == :log_so_to_si
            hlines(0; color=(:gray, 0.5), linewidth=4)
        end

        Label(
            fig[plot_layout_xi..., Top()],
            labels[xi],
            fontsize = 8pt,
            font = :bold,
            padding = (0, 5, 5, 0),
            halign = :center
        )
    end
    n_fig_row = first(fldmod1(length(gdf), n_col))
    Label(
        fig[1:n_fig_row, 0],
        ylabel,
        rotation= pi/2,
        fontsize=10pt
    )

    #linkaxes!(filter(x -> x isa Axis, fig.content)...)
    #resize_to_layout!(fig)

    display(fig)

    return fig
end

function grouped_cluster_ridgeline_plot(
    dataframe,
    cluster_col,
    grouping,
    variable;
    title="",
    xlabel="Cluster",
    ylabel="",
    fig_sizes=fig_sizes,
    fontsize=fontsize,
    datalimits=(-Inf,Inf),
    overlap=1
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
        multi_axis=false
    )

    #labels = label_lines.([first(df[:, grouping]) for df in gdf]; l_length=17)
    labels = [first(df[:, grouping]) for df in gdf]
    colors = [:green, :orange, :blue]
    ax = Axis(
        fig[1,1],
        title=title,
        ylabel=ylabel,
        xlabel=xlabel,
        limits = (extrema(dataframe[:, variable]), nothing)
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
            d = groupdf[groupdf[:, cluster_col] .== cluster, variable]
            cluster_color = colors[j]

            violin!(
                fill(y_offset, length(d)), 
                d, 
                color=(cluster_color, 0.2), 
                orientation=:horizontal, 
                show_median=true,
                mediancolor=cluster_color,
                side = :right,
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
        fig[1,2],
        legend_entries,
        ["Low", "Medium", "High"],
        nbanks=1
    )


    display(fig)

    return fig
end

grouping_full_names = Dict(
    :bioregion => "Bioregion",
    :management_area => "Management Area",
    :gbr => ""
)

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
    groups_too_few_clusters = grouping_counts(grouping, analysis_layers, "$(GCM)_$(grouping)_clusters", 3, 5)
    analysis_layers = analysis_layers[analysis_layers[:, grouping] .∉ [groups_too_few_clusters], :]
    rel_cover = rel_cover = rel_cover[:, rel_cover.locations .∈ [analysis_layers.UNIQUE_ID]]
    dhw_ts = dhw_ts[:, dhw_ts.locations .∈ [rel_cover.locations]]

    mean_dhw_violin = grouped_cluster_ridgeline_plot(
        analysis_layers,
        Symbol("$(GCM)_$(grouping)_clusters"),
        grouping, Symbol("$(GCM)_mean_dhw");
        xlabel="mean DHW", ylabel="$(grouping_fn)", overlap = overlap
    );
    save(
        joinpath(fig_out_dir, "$(grouping)", "mean_dhw_$(grouping)_violin.png"), 
        mean_dhw_violin, 
        px_per_unit = dpi
    )

    # dhw_tol_violin = grouped_cluster_ridgeline_plot(
    #     analysis_layers,
    #     Symbol("$(GCM)_$(grouping)_clusters"),
    #     grouping, Symbol("$(GCM)_mean_DHW_tol");
    #     xlabel="mean reef DHW tolerance", ylabel="$(grouping_fn)", overlap = overlap
    # )
    # save(
    #     joinpath(fig_out_dir, "$(grouping)", "dhw_tolerance_$(grouping)_violin.png"),
    #     dhw_tol_violin,
    #     px_per_unit = dpi
    # )

    so_to_si_violin = grouped_cluster_ridgeline_plot(
        analysis_layers,
        Symbol("$(GCM)_$(grouping)_clusters"),
        grouping, :log_so_to_si;
        xlabel="Log source to sink ratio", ylabel="$(grouping_fn)", overlap = overlap
    );
    save(
        joinpath(fig_out_dir, "$(grouping)", "so_to_si_$(grouping)_violin.png"), 
        so_to_si_violin, 
        px_per_unit = dpi
    )

    total_strength_violin = grouped_cluster_ridgeline_plot(
        analysis_layers,
        Symbol("$(GCM)_$(grouping)_clusters"),
        grouping, :log_total_strength;
        xlabel="Log total connectivity strength", ylabel="$(grouping_fn)", overlap = overlap
    );
    save(
        joinpath(fig_out_dir, "$(grouping)", "log_total_strength_$(grouping)_violin.png"), 
        total_strength_violin,
        px_per_unit = dpi
    )

    initial_proportion_violin = grouped_cluster_ridgeline_plot(
        analysis_layers,
        Symbol("$(GCM)_$(grouping)_clusters"),
        grouping, Symbol("initial_proportion");
        xlabel="Initial proportion of habitable area occupied", ylabel="$(grouping_fn)", overlap = overlap
    );
    save(
        joinpath(fig_out_dir, "$(grouping)", "initial_proportion_$(grouping)_violin.png"), 
        initial_proportion_violin, 
        px_per_unit = dpi
    )

    dhw_cover_cor_violin = grouped_cluster_ridgeline_plot(
        analysis_layers,
        Symbol("$(GCM)_$(grouping)_clusters"),
        grouping, Symbol("$(GCM)_dhw_cover_cor");
        xlabel="Total coral cover - DHW correlation", ylabel="$(grouping_fn)", overlap = overlap
    );
    save(
        joinpath(fig_out_dir, "$(grouping)", "dhw_cover_cor_$(grouping)_violin.png"), 
        dhw_cover_cor_violin, 
        px_per_unit = dpi
    )

    analysis_layers_depth = analysis_layers[analysis_layers.depth_mean .!= 7, :]
    groups_too_few_clusters_depth = grouping_counts(
        grouping, 
        analysis_layers_depth, 
        "$(GCM)_$(grouping)_clusters", 
        3,
        5
    )
    analysis_layers_depth = analysis_layers_depth[
        analysis_layers_depth[:, grouping] .∉ [groups_too_few_clusters_depth], :
    ]

    depth_median_violin = grouped_cluster_ridgeline_plot(
        analysis_layers_depth,
        Symbol("$(GCM)_$(grouping)_clusters"),
        grouping, Symbol("depth_med");
        xlabel="Median Depth [m]", ylabel="$(grouping_fn)", overlap = overlap
    );
    save(
        joinpath(fig_out_dir, "$(grouping)", "depth_$(grouping)_violin.png"), 
        depth_median_violin, 
        px_per_unit = dpi
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
        px_per_unit = dpi
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
        px_per_unit = dpi
    )

    return nothing
end

function grouped_GCM_cluster_timeseries_plots(
    timeseries_array,
    dataframe,
    cluster_col,
    grouping,
    length_t;
    fig_sizes=fig_sizes,
    fontsize=fontsize
)
    fig_x_size = fig_sizes["timeseries_width"]
    fig_y_size = fig_sizes["timeseries_height"]
    n_col = optimum_columns(size(unique(dataframe[:, grouping]), 1))
    fig, gdf, plot_layout = _setup_grouped_figure(
        dataframe,
        grouping;
        x_fig_size=fig_x_size,
        y_fig_size=fig_y_size,
        fontsize=fontsize,
        order=[:management_area, :GCM]
    )

    labels = label_lines.("$(first(df.management_area)) - $(first(df.GCM))" for df in gdf; l_length=15)

    for (xi, groupdf) in enumerate(gdf)
        # println("$(xi)")
        plot_layout_xi = plot_layout[xi]
        groupdf = sort(groupdf, cluster_col)

        # Ensure that relative cover timeseries match the cluster allocations from groupdf
        gcm = first(groupdf.GCM)
        group_timeseries = timeseries_array[GCM = (timeseries_array.GCM .== gcm)][:, :, 1]
        group_timeseries = group_timeseries[locations = (group_timeseries.locations .∈ [groupdf.UNIQUE_ID])]
        group_timeseries = group_timeseries[length_t, indexin(groupdf.UNIQUE_ID, String.(group_timeseries.locations))]
        
        # group_timeseries_less_than_5 = [all(group_timeseries[:, i].data .< 5) for i in 1:size(group_timeseries, 2)]
        # group_timeseries_less_than_5_ind = findall(group_timeseries_less_than_5)
        # group_timeseries = group_timeseries[:, group_timeseries_less_than_5_ind]
        groupdf = groupdf[groupdf.UNIQUE_ID .∈ [group_timeseries.locations], :]

        clusters = Int64.(groupdf[:, cluster_col])

        ADRIA.viz.clustered_scenarios!(
            fig[plot_layout_xi...],
            group_timeseries,
            clusters;
            opts = Dict{Symbol, Any}(:legend => false),
            axis_opts = Dict(
                :title => labels[xi],
                # :yticks => [0,2,4],
                :xticks => ([1,10,20,30,40,50], string.([2025,2035,2045,2055,2065,2075]))
            )
        )
    end

    axes_after_1 = filter(x -> x isa Axis, fig.content)[2:end]
    map(x -> hidexdecorations!(x; grid = false, ticks = false), axes_after_1)
    map(x -> hideydecorations!(x, ticks = false, ticklabels = false, grid = false), axes_after_1)
    
    # Tweak layout defaults
    # fig.layout.colgap = 5  # Default ~10, reduce to make plots tighter horizontally
    # fig.layout.rowgap = 5  # Same for vertical spacing
    # fig.layout.outerpadding = 5  # Reduce space around the entire figure

    linkaxes!(filter(x -> x isa Axis, fig.content)...)
    display(fig)

    return fig
end

function carbonate_budget_variable_scatter(
    long_df,
    color_variable_col,
    year_col,
    carbonate_budget_col,
    year_variable_correlation;
    xlabel="Carbonate budget threshold [%] (correlation)",
    ylabel="Years above carbonate budget threshold",
    color_label="",
    fig_sizes=fig_sizes,
    fontsize=fontsize,
    alpha=0.5
)
    x_fig_size = fig_sizes["carb_width"]
    y_fig_size = fig_sizes["carb_height"]

    labels = unique(long_df[:, carbonate_budget_col])
    labels = ["$(lab)\n($(round(year_variable_correlation[lab], digits=2)))" for lab in labels]

    # Get the default colormap (as a gradient with 256 colors)
    base_cmap = cgrad(:viridis, 256; rev=true)
    # Convert all colors to RGBA with alpha = 0.4 (adjust as needed)
    transparent_colors = [RGBA(c.r, c.g, c.b, alpha) for c in base_cmap.colors]
    # Wrap into a ColorScheme object
    transparent_cmap = ColorScheme(transparent_colors)


    fig = Figure(size = (x_fig_size, y_fig_size), fontsize = fontsize)
    ax = Axis(
        fig[1,1],
        ylabel = ylabel,
        xlabel = xlabel,
        xticks = (1:1:11, labels)
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

    # scat = scatter!(
    #     ACCESS_ESM1_5_long.variable,
    #     ACCESS_ESM1_5_long.value,
    #     color=ACCESS_ESM1_5_long.depth_med,
    #     alpha = 0.5
    # )
    Colorbar(fig[1,2], rain, label = color_label)

    display(fig)

    return fig
end