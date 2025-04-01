"""
File includes the plotting functions for ADRIA-reef-indicators analysis.
"""

using Printf

function _convert_plottable(gdf::Union{DataFrame,DataFrameRow}, geom_col::Symbol)
    local plottable
    try
        if gdf isa DataFrame
            plottable = GeoMakie.geo2basic(AG.forceto.(gdf[!, geom_col], AG.wkbPolygon))
        else
            plottable = GeoMakie.geo2basic(AG.forceto(gdf[geom_col], AG.wkbPolygon))
        end
    catch
        # Column is already in a plottable form, or some unrelated error occurred
        if gdf isa DataFrame
            plottable = gdf[:, geom_col]
        else
            plottable = [gdf[geom_col]]
        end
    end

    return plottable
end

"""
    plot_map(gdf::DataFrame; geom_col::Symbol=:geometry, color_by::Symbol)

Convenience plot function.

# Arguments
- `gdf` : GeoDataFrame
- `color_by` : Column name holding factor to color reefs by (e.g. :management_area)
- `geom_col` : Column name holding geometries to plot
"""
function plot_map(gdf::Union{DataFrame,DataFrameRow}; geom_col::Symbol=:geometry)
    f = Figure(; size=(600, 900))
    ga = GeoAxis(
        f[1, 1];
        dest="+proj=latlong +datum=WGS84",
        xlabel="Longitude",
        ylabel="Latitude",
        xticklabelpad=15,
        yticklabelpad=40,
        xticklabelsize=10,
        yticklabelsize=10,
        aspect=AxisAspect(0.75),
        xgridwidth=0.5,
        ygridwidth=0.5,
    )

    plottable = _convert_plottable(gdf, geom_col)
    poly!(ga, plottable)

    display(f)

    return f, ga
end

function plot_map!(ga::GeoAxis, gdf::DataFrame; geom_col=:geometry, color=nothing)::Nothing

    plottable = _convert_plottable(gdf, geom_col)
    if !isnothing(color)
        poly!(ga, plottable, color=color)
    else
        poly!(ga, plottable)
    end

    # Set figure limits explicitly
    xlims!(ga)
    ylims!(ga)

    return nothing
end

function plot_map!(gdf::DataFrame; geom_col=:geometry, color=nothing)::Nothing
    return plot_map!(current_axis(), gdf; geom_col=geom_col, color=color)
end

function plot_map(gdf::Union{DataFrame,DataFrameRow}, color_by::Symbol; geom_col::Symbol=:geometry)
    f = Figure(; size=(600, 900))
    ga = GeoAxis(
        f[1, 1];
        dest="+proj=latlong +datum=WGS84",
        xlabel="Longitude",
        ylabel="Latitude",
        xticklabelpad=15,
        yticklabelpad=40,
        xticklabelsize=10,
        yticklabelsize=10,
        aspect=AxisAspect(0.75),
        xgridwidth=0.5,
        ygridwidth=0.5,
    )

    plottable = _convert_plottable(gdf, geom_col)

    # Define the unique colors and names for each level of factor color_by.
    # Use a different color palette for factors with high numbers of levels
    # (this palette is not as good for visualisation).
    if size(unique(gdf[:, color_by]),1) <= 20
        palette = ColorSchemes.tableau_20.colors
    else
        palette = ColorSchemes.flag_ec.colors
    end

    color_indices = groupindices(DataFrames.groupby(gdf, color_by))
    names = unique(DataFrame(indices=color_indices, names=gdf[:, color_by]))

    # Create the unique legend entries for each level of color_by
    unique_names = names.names
    legend_entries = []
    for name in eachrow(names)
        col = palette[name.indices]
        LE = PolyElement(; color=col)
        push!(legend_entries, [LE])
    end

    polys = poly!(ga, plottable, color=palette[color_indices])

    Legend(f[2, 1], legend_entries, unique_names, nbanks=3, tellheight=true,
    tellwidth=false, orientation=:horizontal, labelsize=10)

    display(f)
    return f, ga
end

function colorscheme_alpha(cscheme::ColorScheme, alpha = 0.5)
    ncolors = length(cscheme)

    return ColorScheme([RGBA(get(cscheme, k), alpha) for k in range(0, 1, length=ncolors)])
end

"""
    plot_lines(
    df::DataFrame,
    color_by::Symbol,
    y_start_col,
    y_end_col,
    x_lab::String,
    y_lab::String,
    alpha=1
    )

Convenience plot function for plotting timeseries lines with categorical colours.

# Arguments
- `df` : DataFrame with wide format (each unique line {reef} in a plot occupies it's own row in df)
- `color_by` : Column name holding factor to color reefs by (e.g. :management_area)
- `y_start_col` : Column name holding the first timepoint of the series (often 1 or 2)
- `y_end_col` : Column name holding the last timepoint of the series (often 79)
- `x_lab` : Label for X Axis
- `y_lab` : Label for Y Axis
- `alpha` : Alpha value for transparency (0-1)
"""
function plot_lines(
    df::DataFrame,
    color_by::Symbol,
    y_start_col,
    y_end_col,
    x_lab::String,
    y_lab::String,
    alpha=1
    )

    f = Figure()
    ax = Axis(f[1,1]; xlabel=x_lab, ylabel=y_lab)

    # Define the unique colors and names for each level of factor color_by.
    # Use a different color palette for factors with high numbers of levels
    # (this palette is not as good for visualisation).
    if size(unique(df[:, color_by]),1) == 2
        alph_palette = colorscheme_alpha(ColorSchemes.tableau_20, alpha)[2:3]
        palette = ColorSchemes.tableau_20.colors[2:3]
    elseif size(unique(df[:, color_by]),1) <= 20
        alph_palette = colorscheme_alpha(ColorSchemes.tableau_20, alpha)
        palette = ColorSchemes.tableau_20.colors
    else
        alph_palette = colorscheme_alpha(ColorSchemes.flag_ec, alpha)
        palette = ColorSchemes.flag_ec.colors
    end

    color_indices = groupindices(DataFrames.groupby(df, color_by))
    color_names = unique(DataFrame(indices=color_indices, names=df[:, color_by]))

    # Create the unique legend entries for each level of color_by
    unique_names = color_names.names
    legend_entries = []
    for name in eachrow(color_names)
        col = palette[name.indices]
        LE = PolyElement(; color=col)
        push!(legend_entries, [LE])
    end

    lines = series!(ax, Matrix(DataFrames.select(df, DataFrames.Between(Symbol(y_start_col), Symbol(y_end_col)))), color=alph_palette[color_indices])

    Legend(f[2, 1], legend_entries, unique_names, nbanks=3, tellheight=true,
    tellwidth=false, orientation=:horizontal, labelsize=10)

    display(f)
    return f, ax
end

# Checking how prominent the reefs with correlation at later lags are at earlier lags
using Colors, Random

# For horizintal and vertical text alignment:
justifyme(θ) = (0≤θ<π/2 || 3π/2<θ≤2π) ? :left : (π/2<θ<3π/2) ? :right : :center
justifymeV(θ) = π/4≤θ≤3π/4 ? :bottom : 5π/4<θ≤7π/4 ? :top : :center

# Radar plot function:
function radarplot(ax::Axis, v, val_labels; p_grid = maximum(v) * (1.0:length(val_labels)) / length(val_labels), title = "", labels = eachindex(v), labelsize = 1, points=true, spokeswidth= 1.5, spokescolor=:salmon, fillalpha=0.2, linewidth=1.5)
    # Axis attributes
    ax.xgridvisible = false
    ax.ygridvisible = false
    ax.xminorgridvisible = false
    ax.yminorgridvisible = false
    ax.leftspinevisible = false
    ax.rightspinevisible = false
    ax.bottomspinevisible = false
    ax.topspinevisible = false
    ax.xminorticksvisible = false
    ax.yminorticksvisible = false
    ax.xticksvisible = false
    ax.yticksvisible = false
    ax.xticklabelsvisible = false
    ax.yticklabelsvisible = false
    ax.aspect = DataAspect()
    ax.title = title
    #
    l = length(v)
    rad = (0:(l-1)) * 2π / l
    # Point coordinates
    x = v .* cos.(rad) .- 0.05
    y = v .* sin.(rad) .- 0.05
    # if p_grid != 0
        # Coordinates for radial grid
        xa = maximum(p_grid) * cos.(rad) * 1.1
        ya = maximum(p_grid) * sin.(rad) * 1.1
        # Coordinates for polar grid text
        radC = (rad[Int(round(l / 2))] + rad[1 + Int(round(l / 2))]) / 2.0
        xc = p_grid * cos(radC)
        yc = p_grid * sin(radC)
        for i in p_grid
            poly!(ax, Circle(Point2f(0, 0), i), color = :transparent, strokewidth=1, strokecolor=ax.xgridcolor)
        end
        text!(ax, xc, yc, text=val_labels, fontsize = 12, align = (:center, :baseline), color=ax.xlabelcolor)
        arrows!(ax, zeros(l), zeros(l), xa, ya, color=ax.xgridcolor, linestyle=:solid, arrowhead=' ')
        if length(labels) == l
            for i in eachindex(rad)
                text!(ax, xa[i], ya[i], text=string(labels[i]), fontsize = labelsize, markerspace = :data, align = (justifyme(rad[i]), justifymeV(rad[i])), color=ax.xlabelcolor)
            end
        elseif length(labels) > 1
            printstyled("WARNING! Labels omitted:  they don't match with the points ($(length(labels)) vs $l).\n", bold=true, color=:yellow)
        end
    # end
    pp = scatter!(ax, [(x[i], y[i]) for i in eachindex(x)], color=RGB{Float32}(0.4, 0.4, 0.4))
    cc = to_value(pp.color)
    m_color = RGBA{Float32}(comp1(cc), comp2(cc), comp3(cc), fillalpha)
    s_color = RGB{Float32}(comp1(cc), comp2(cc), comp3(cc))
    pp.color = m_color
    pp.strokecolor = s_color
    pp.strokewidth= linewidth
    arrows!(ax, zeros(l), zeros(l), x, y, color=spokescolor, linewidth=spokeswidth, arrowhead=' ')
    if points
        scatter!(ax, x, y)
    end
    ax
end

function radarplot!(ax, v)

    l = length(v)
    rad = (0:(l-1)) * 2π / l

    x = v .* cos.(rad) .- 0.15
    y = v .* sin.(rad) .- 0.15

    pp = scatter!(ax, [(x[i], y[i]) for i in eachindex(x)], color=RGB{Float32}(0.4, 0.4, 0.4))

    ax
end

"""
    radarplot_df(df, cols, val_labels, labels; spokeswidth = 0, labelsize = 1)

Plots a radar plot intended to visualise each of the bellwether reefs and check whether their
correlation occurs at all lags or only some lags.

# Arguments
- `df` : dataframe containing a row for each reef and columns for each relevant lag. Each
lag column has corresponding integer values if a reef is a bellwether reef at that lag, e.g.
if reef-1 is a bellwether reef at lag5 the value of df[reef-1, lag5] = 5.
- `cols` : vector of column names containing the lag values for each reef.
- `val_labels` : labels for the values of each ring in the plot.
- `labels` : Labels for each reef. (RME_UNIQUE_ID or similar identifier)
"""
function radarplot_df(df, cols, val_labels, labels; spokeswidth = 0, labelsize = 1)
    fig = Figure()
    ax = Axis(fig[1,1])
    max_vals = [maximum(df[:, col]) for col in cols]
    initial_column = df[:, cols[argmax(max_vals)]]

    f = radarplot(ax, initial_column, val_labels; labels = labels, spokeswidth = spokeswidth, labelsize = labelsize)

    for col in cols
        f = radarplot!(ax, df[:, col])
    end

    display(fig)

    return fig
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

function explore_regression_scatter_plots(x, y; color=:blue, xlab="", ylab="", title="")

    fig = Figure()
    ax = Axis(
        fig[1,1];
        xlabel = xlab,
        ylabel = ylab,
        title = title
    )
    f = scatter!(
        x, y; color=color
    )

    return fig
end

function label_lines(label)
    if length(label) > 25
        label_new = replace(label, r"(.{25} )" => s"\1\n")

        return label_new
    end

    return label
end

function _axis_size(gdf, x_fig_size, y_fig_size, n_col)
    xsize = x_fig_size / (n_col)*1.5
    n_fig_row = first(fldmod1(length(gdf), n_col))
    ysize = y_fig_size / (n_fig_row*1.5)

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

function _setup_grouped_figure(dataframe, grouping; x_fig_size=2130, y_fig_size=1500, marker=LineElement)
    dataframe = sort(dataframe, :management_area; rev=true)
    
    gdf = DataFrames.groupby(dataframe, grouping)

    fig = Figure(backgroundcolor = RGBf(0.98, 0.98, 0.98), size = (x_fig_size, y_fig_size))

    n_col = optimum_columns(length(unique(dataframe[:, grouping])))
    plot_layout = figure_layouts(length(gdf), n_col)

    # Configure the legend for the figure
    last_figure = last(plot_layout)
    if last(last_figure) < n_col
        legend_position = (last_figure[1], last_figure[2]+1:n_col)
    else
        legend_position = (last_figure[1] + 1, 1:n_col)
    end
    
    legend_entries = []
    for (i, col) in enumerate([:blue, :orange, :green])
        LE = marker(; color=col, marker=:circle)
        push!(legend_entries, [LE])
    end

    Legend(
        fig[legend_position...],
        legend_entries,
        ["High", "Medium", "Low"],
        nbanks=1
    )

    return fig, gdf, plot_layout
end

function _setup_grouped_axes(fig, plot_layout_xi, xticks; ylabel="", xlabel="", title="", xsize=220, ysize=150, background_color=:white)
    # if plot_layout_xi == (1,1)
    #     ax = Axis(
    #                 fig[plot_layout_xi...];
    #                 backgroundcolor=background_color,
    #                 xlabel = xlabel,
    #                 xlabelsize = 13,
    #                 ylabelsize = 13,
    #                 ylabel = ylabel,
    #                 xticks = xticks,
    #                 xticklabelsize=12,
    #                 yticklabelsize=12,
    #                 title=title,
    #                 width=xsize,
    #                 height=ysize
    #     )
    #     return ax
    # end

    ax = Axis(
        fig[plot_layout_xi...];
        backgroundcolor=background_color,
        xticklabelsize=16,
        yticklabelsize=16,
        xticks = xticks,
        xticksvisible=false,
        title=title,
        width=xsize,
        height=ysize
    )
    return ax
end

function grouped_violin_plots(
    dataframe,
    bellwether_reefs_col,
    grouping,
    variable;
    ylabel="", xlabel="", title="",
    x_fig_size=1500,
    y_fig_size=1000,
    datalimits=(-Inf, Inf)
)
    n_col = optimum_columns(length(unique(dataframe[:, grouping])))
    fig, gdf, plot_layout, colors, categories = _setup_grouped_figure(
        dataframe,
        bellwether_reefs_col,
        grouping;
        x_fig_size=x_fig_size,
        y_fig_size=y_fig_size
    )
    xsize, ysize = _axis_size(gdf, x_fig_size, y_fig_size, n_col)
    xticks = (1:2, levels(categories))

    labels = label_lines.([first(df.bioregion) for df in gdf])

    for (xi, groupdf) in enumerate(gdf)
        plot_layout_xi = plot_layout[xi]
        categories = categorical(groupdf[:, bellwether_reefs_col])

        bellwether_var_med = median(groupdf[groupdf[:, bellwether_reefs_col] .== "bellwether", variable])
        non_bellwether_var_75 = quantile(groupdf[groupdf[:, bellwether_reefs_col] .== "non-bellwether", variable], 0.75)
        non_bellwether_var_25 = quantile(groupdf[groupdf[:, bellwether_reefs_col] .== "non-bellwether", variable], 0.25)

        if bellwether_var_med < non_bellwether_var_25
            background_color = (:royalblue, 0.2)
        elseif bellwether_var_med > non_bellwether_var_75
            background_color = (:goldenrod1, 0.2)
        else
            background_color = :white
        end

        ax = _setup_grouped_axes(
            fig,
            plot_layout_xi,
            xticks;
            ylabel=ylabel,
            xlabel=xlabel,
            title=title,
            xsize=xsize,
            ysize=ysize,
            background_color=background_color
        )

        f = violin!(
            ax,
            categories.refs,
            groupdf[:, variable];
            color=colors[indexin(categories, unique(categories))],
            show_median=true,
            datalimits=datalimits
        )
        f = rainclouds!(
            ax,
            categories.refs,
            groupdf[:, variable];
            color=:black,
            markersize=5,
            jitter_width=0.27,
            side_nudge=0.001,
            plot_boxplots=false,
            clouds=nothing
        )

        if variable == :so_to_si
            hlines!(1; color=(:gray, 0.5), linewidth=3)
        end

        Label(
            fig[plot_layout_xi..., Top()],
            labels[xi],
            fontsize = 16,
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
        fontsize=18
    )

    #linkaxes!(filter(x -> x isa Axis, fig.content)...)
    resize_to_layout!(fig)

    display(fig)

    return fig
end

function timeseries_xticks(length_t, years)
    length_range = first(length_t):10:last(length_t)
    years_length = collect(years)[length_range]

    return (length_range, string.(years_length))
end


# function grouped_violin_plots(
#     dataframe,
#     bellwether_reefs_col,
#     grouping,
#     variable;
#     ylabel="", xlabel="", title="",
#     xticks=(1:length(unique(dataframe[:, bellwether_reefs_col])), unique(dataframe[:, bellwether_reefs_col])),
#     x_fig_size=2000,
#     y_fig_size=1700,
#     datalimits=extrema
# )
#     local violin_plot
#     try
#         violin_plot = _grouped_violin_plots(
#             dataframe,
#             bellwether_reefs_col,
#             grouping,
#             variable;
#             ylabel="", xlabel="", title="",
#             xticks=(1:length(unique(dataframe[:, bellwether_reefs_col])), unique(dataframe[:, bellwether_reefs_col])),
#             x_fig_size=2000,
#             y_fig_size=1700,
#             datalimits=(-Inf, Inf)
#         )
#     catch y
#         if isa(y, BoundsError) | isa()
#             violin_plot = _grouped_violin_plots(
#                 dataframe,
#                 bellwether_reefs_col,
#                 grouping,
#                 variable;
#                 ylabel="", xlabel="", title="",
#                 xticks=(1:length(unique(dataframe[:, bellwether_reefs_col])), unique(dataframe[:, bellwether_reefs_col])),
#                 x_fig_size=2000,
#                 y_fig_size=1700,
#                 datalimits=(-Inf, Inf)
#             )
#         else
#             rethrow(y)
#         end
#     end

#     return violin_plot
# end

function grouped_timeseries_plots(
    dataframe,
    bellwether_reefs_col,
    #correlation_col,
    grouping,
    length_t,
    lag;
    ylabel="", xlabel="", title="",
    x_fig_size=1500,
    y_fig_size=1000
)
    n_col = n_col = optimum_columns(length(unique(dataframe[:, grouping])))
    fig, gdf, plot_layout, colors = _setup_grouped_figure(
        dataframe,
        bellwether_reefs_col,
        grouping;
        x_fig_size=x_fig_size,
        y_fig_size=y_fig_size
    )
    xsize, ysize = _axis_size(gdf, x_fig_size, y_fig_size, n_col)
    xticks = timeseries_xticks(length_t, 2025:2099)
    #xticks = (1:10:55, ["$(year)" for year in 2025:10:2075])

    #labels = label_lines.([_extract_name_and_correlation(df, bellwether_reefs_col, correlation_col, grouping) for df in gdf])
    labels = label_lines.(first(df[:, grouping]) for df in gdf)

    for (xi, groupdf) in enumerate(gdf)
        plot_layout_xi = plot_layout[xi]
        categories = categorical(groupdf[:, bellwether_reefs_col])

        target_reefs = Matrix(
        DataFrames.select(
            groupdf[groupdf[:, bellwether_reefs_col] .== "bellwether", :],
            DataFrames.Between(Symbol(length_t[1]), Symbol(length_t[2]-lag))
            )
        )'
        target_reefs_median, target_reefs_lb, target_reefs_ub = bootstrap_mean_ts(target_reefs)

        non_target_reefs = Matrix(
            DataFrames.select(
                groupdf[groupdf[:, bellwether_reefs_col] .== "non-bellwether", :],
                DataFrames.Between(Symbol(length_t[1]+lag), Symbol(length_t[2]))
            )
        )'
        # if lag > 0
        #     non_target_reefs_buffer = fill(missing, lag, size(non_target_reefs, 2))
        #     non_target_reefs = vcat(non_target_reefs, non_target_reefs_buffer)
        # end
        non_target_reefs_median, non_target_reefs_lb, non_target_reefs_ub = bootstrap_mean_ts(non_target_reefs)

        ax = _setup_grouped_axes(
            fig,
            plot_layout_xi,
            xticks;
            ylabel=ylabel, xlabel=xlabel, title=title, xsize=xsize, ysize=ysize
        )

        #band_indices = first(length_t) + lag:last(length_t) + lag
        band_indices = first(length_t):(last(length_t)-lag)

        band!(band_indices, non_target_reefs_lb, non_target_reefs_ub; color=(colors[1], 0.2))
        #band!(first(length_t):last(length_t), non_target_reefs_lb, non_target_reefs_ub; color=(colors[1], 0.3))
        series!(non_target_reefs', solid_color=(colors[1], 0.2))
        series!(non_target_reefs_median', solid_color=:blue)

        
        #bands = band!(band_indices, target_reefs_lb[band_indices], target_reefs_ub[band_indices]; color=(colors[2], 0.3))
        band!(band_indices, target_reefs_lb, target_reefs_ub; color=(colors[2], 0.3))
        series!(ax, target_reefs', solid_color=(colors[2], 0.5))
        series!(target_reefs_median', solid_color=:red)
      

        Label(
            fig[plot_layout_xi..., Top()],
            labels[xi],
            fontsize = 16,
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
        fontsize=18
    )

    #linkaxes!(filter(x -> x isa Axis, fig.content)...)
    resize_to_layout!(fig)

    display(fig)

    return fig
end

function violin_limits(variable_vector)
    if stat_range(variable_vector) < 2
        max = ceil(maximum(variable_vector); digits=2)
        min = floor(minimum(variable_vector); digits=2)

        return (min, max)
    elseif (stat_range(variable_vector) < 0.5) & (minimum(variable_vector) < 0.000001)
        return (-0.005, ceil(maximum(variable_vector); digits=2))
    end

    return (-Inf, Inf)
end

# function violin_plots!(
#     fig,
#     gdf,
#     bellwether_reefs_col,
#     variable,
#     labels,
#     plot_layout,
#     xsize,
#     ysize,
#     colors;
#     xlabel=xlabel,
#     ylabel=ylabel,
#     title="",
#     xticks=unique(dataframe[:, bellwether_reefs_col])
# )
#     for (xi, groupdf) in enumerate(gdf)
#         categories = categorical(groupdf[:, bellwether_reefs_col])

#         bellwether_var_med = median(groupdf[groupdf[:, bellwether_reefs_col] .== "bellwether", variable])
#         non_bellwether_var_75 = quantile(groupdf[groupdf[:, bellwether_reefs_col] .== "non-bellwether", variable], 0.75)
#         non_bellwether_var_25 = quantile(groupdf[groupdf[:, bellwether_reefs_col] .== "non-bellwether", variable], 0.25)

#         if bellwether_var_med < non_bellwether_var_25
#             background_color = (:royalblue, 0.2)
#         elseif bellwether_var_med > non_bellwether_var_75
#             background_color = (:goldenrod1, 0.2)
#         else
#             background_color = :white
#         end

#         if xi == 1
#             ax = Axis(
#                 fig[plot_layout[xi]...];
#                 backgroundcolor=background_color,
#                 xlabel = xlabel,
#                 xlabelsize = 13,
#                 ylabelsize = 13,
#                 ylabel = ylabel,
#                 xticks = (unique(categories.refs), xticks),
#                 xticklabelsize=12,
#                 yticklabelsize=12,
#                 title=title,
#                 width=xsize,
#                 height=ysize
#             )
#         else
#             ax = Axis(
#                 fig[plot_layout[xi]...];
#                 backgroundcolor=background_color,
#                 xticklabelsize=12,
#                 yticklabelsize=12,
#                 xticks = (unique(categories.refs), ["",""]),
#                 xticksvisible=false,
#                 title=title,
#                 width=xsize,
#                 height=ysize
#             )
#         end

#         f = violin!(
#             ax,
#             categories.refs,
#             groupdf[:, variable];
#             color=colors[indexin(categories, unique(categories))],
#             show_median=true,
#             datalimits=extrema
#         )
#         f = rainclouds!(
#             ax,
#             categories.refs,
#             groupdf[:, variable];
#             color=:black,
#             markersize=5,
#             jitter_width=0.27,
#             side_nudge=0.001,
#             plot_boxplots=false,
#             clouds=nothing
#         )

#         if variable == :so_to_si
#             hlines!(1; color=(:gray, 0.5), linewidth=3)
#         end

#         Label(
#             fig[plot_layout[xi]..., Top()],
#             labels[xi],
#             fontsize = 11,
#             font = :bold,
#             padding = (0, 5, 5, 0),
#             halign = :center
#         )
#     end

#     if variable ∈ ["$(GCM)_mean_dhw", "$(GCM)_initial_coral_cover", :so_to_si]
#         linkaxes!(filter(x -> x isa Axis, fig.content)...)
#     end

#     return fig
# end

# function timeseries_plots!(
#     fig,
#     gdf,
#     bellwether_reefs_col,
#     plot_layout,
#     length_t,
#     lag,
#     labels,
#     xsize,
#     ysize;
#     xlabel=xlabel,
#     ylabel=ylabel,
#     title="",
#     xticks=unique(dataframe[:, bellwether_reefs_col])
# )
#     for (xi, groupdf) in enumerate(gdf)
#         categories = categorical(groupdf[:, bellwether_reefs_col])

#         target_reefs = Matrix(
#         DataFrames.select(
#             groupdf[groupdf[:, bellwether_reefs_col] .== "bellwether", :],
#             DataFrames.Between(Symbol(length_t[1]), Symbol(length_t[2]))
#             )
#         )'
#         if lag > 0
#             target_reefs_buffer = fill(missing, 4, size(target_reefs, 2))
#             target_reefs = vcat(target_reefs_buffer, target_reefs)
#         end
#         target_reefs_median, target_reefs_lb, target_reefs_ub = bootstrap_median_ts(target_reefs)

#         non_target_reefs = Matrix(
#             DataFrames.select(
#                 groupdf[groupdf[:, bellwether_reefs_col] .== "non-bellwether", :],
#                 DataFrames.Between(Symbol(length_t[1]), Symbol(length_t[2]))
#             )
#         )'
#         non_target_reefs_median, non_target_reefs_lb, non_target_reefs_ub = bootstrap_median_ts(non_target_reefs)

#         if xi == 1
#             ax = Axis(
#                 fig[plot_layout[xi]...];
#                 xlabel = xlabel,
#                 xlabelsize = 13,
#                 ylabelsize = 13,
#                 ylabel = ylabel,
#                 xticklabelsize=12,
#                 yticklabelsize=12,
#                 title=title,
#                 width=xsize,
#                 height=ysize
#             )
#         else
#             ax = Axis(
#                 fig[plot_layout[xi]...];
#                 xticklabelsize=12,
#                 yticklabelsize=12,
#                 title=title,
#                 width=xsize,
#                 height=ysize
#             )
#         end

#         band_indices = first(length_t) + lag:last(length_t) + lag
#         bands = band!(band_indices, target_reefs_lb[band_indices], target_reefs_ub[band_indices]; color=(colors[2], 0.3))
#         #series!(ax, target_reefs', solid_color=(colors[2], 0.3))
#         series!(target_reefs_median', solid_color=:red)

#         band!(first(length_t):last(length_t), non_target_reefs_lb, non_target_reefs_ub; color=(colors[1], 0.3))
#         #series!(non_target_reefs', solid_color=(colors[1], 0.2))
#         series!(non_target_reefs_median', solid_color=:blue)

#         Label(
#             fig[plot_layout[xi]..., Top()],
#             labels[xi],
#             fontsize = 11,
#             font = :bold,
#             padding = (0, 5, 5, 0),
#             halign = :center
#         )
#     end

#     return fig
# end

# function bioregion_grouped_violins(
#     dataframe, bellwether_reefs_col, grouping, variable, ncol;
#     xlabel="Bellwether Reefs",
#     ylabel="Value",
#     xticks=unique(dataframe[:, bellwether_reefs_col]),
#     title=""
# )
#     dataframe = sort(dataframe, [bellwether_reefs_col, :management_area]; rev=true)
#     categories = categorical(dataframe[:, bellwether_reefs_col])

#     gdf = DataFrames.groupby(dataframe, grouping)

#     # if length(gdf) < 27
#     #     labels = string.(collect('a':'z'))[1:length(gdf)]
#     # else
#     #     error("number of bioregions > 26 (number of letters for labelling)")
#     # end
#     labels = label_lines.([first(df.bioregion) for df in gdf])

#     x_fig_size, y_fig_size = (2130, 1500)
#     fig = Figure(backgroundcolor = RGBf(0.98, 0.98, 0.98), size = (x_fig_size, y_fig_size))

#     colors = [Makie.wong_colors(); Makie.wong_colors()];
#     legend_entries = []
#     for (i, col) in enumerate(unique(colors[indexin(categories, unique(categories))]))
#         LE = MarkerElement(; color=col, marker=:circle)
#         push!(legend_entries, [LE])
#     end

#     min_var, max_var = extrema(dataframe[:, variable])
#     # var_difference = max_var - min_var
#     # # if variable ∈ [:mean_dhw, :so_to_si]
#     # #     limits = (nothing, (floor(Int64, min_var - 2), ceil(Int64, max_var + 2)))
#     # #     yticks = round.(Int64, min_var:2:max_var)
#     # # elseif variable == :initial_coral_cover
#     # #     limits = (nothing, (floor(Int64, min_var - 5), ceil(Int64, max_var + 5)))
#     # #     yticks = round.(Int64, min_var:10:max_var)
#     # # elseif variable == :total_strength
#     # #     limits = (nothing, (round(min_var - 0.3, digits=2), round(max_var + 0.3, digits=2)))
#     # #     yticks = round.(min_var:0.5:max_var, digits=1)
#     # # elseif variable ∈ [:dhw_cover_cor, :dhw_evenness_cor]
#     # #     limits = (nothing, (-1.1,1.1))
#     # #     yticks = (-1:0.5:1)
#     # # elseif variable == :conn_score
#     # #     limits = (nothing, nothing)
#     # #     yticks = Makie.automatic
#     # # end

#     # if var_difference < 4
#     #     limits = (nothing, (min_var - (var_difference/5), max_var + (var_difference/5)))
#     #     yticks = round.(min_var:(var_difference / 4):max_var; digits=3)
#     # else
#     #     limits = (nothing, (floor(Int64, min_var - 1), ceil(Int64, max_var + 1)))
#     #     yticks = round.(min_var:(var_difference / 4):max_var; digits=1)
#     # end

#     # if length(gdf) < 10
#     #     xsize, ysize = 180, 120
#     #     label_size = 12
#     # else
#     #     xsize, ysize = 150, 150
#     #     label_size = 12
#     # end

#     xsize = x_fig_size / (ncol*2)
#     n_fig_row = first(fldmod1(length(gdf), ncol))
#     ysize = y_fig_size / (n_fig_row*1.5)

#     for (xi, groupdf) in enumerate(gdf)
#         categories = categorical(groupdf[:, bellwether_reefs_col])

#         bellwether_var_med = median(groupdf[groupdf[:, bellwether_reefs_col] .== "bellwether", variable])
#         non_bellwether_var_75 = quantile(groupdf[groupdf[:, bellwether_reefs_col] .== "non-bellwether", variable], 0.75)
#         non_bellwether_var_25 = quantile(groupdf[groupdf[:, bellwether_reefs_col] .== "non-bellwether", variable], 0.25)

#         if bellwether_var_med < non_bellwether_var_25
#             background_color = (:royalblue, 0.2)
#         elseif bellwether_var_med > non_bellwether_var_75
#             background_color = (:goldenrod1, 0.2)
#         else
#             background_color = :white
#         end

#         if xi == 1
#             ax = Axis(
#                 fig[fldmod1(xi, ncol)...];
#                 backgroundcolor=background_color,
#                 xlabel = xlabel,
#                 xlabelsize = label_size+1,
#                 ylabelsize = label_size+1,
#                 ylabel = ylabel,
#                 xticks = (unique(categories.refs), xticks),
#                 xticklabelsize=label_size,
#                 #limits= limits,
#                 #yticks = yticks,
#                 yticklabelsize=label_size,
#                 title=title,
#                 width=xsize,
#                 height=ysize
#             )
#         else
#             ax = Axis(
#                 fig[fldmod1(xi, ncol)...];
#                 backgroundcolor=background_color,
#                 xticklabelsize=label_size,
#                 yticklabelsize=label_size,
#                 xticks = (unique(categories.refs), ["",""]),
#                 xticksvisible=false,
#                 #limits= limits,
#                 #yticks = yticks,
#                 title=title,
#                 width=xsize,
#                 height=ysize
#             )
#         end

#         f = violin!(
#             ax,
#             categories.refs,
#             groupdf[:, variable];
#             color=colors[indexin(categories, unique(categories))],
#             show_median=true
#         )
#         f = rainclouds!(
#             ax,
#             categories.refs,
#             groupdf[:, variable];
#             color=:black,
#             markersize=5,
#             jitter_width=0.27,
#             side_nudge=0.001,
#             plot_boxplots=false,
#             clouds=nothing
#         )

#         if variable == :so_to_si
#             hlines!(1; color=(:gray, 0.5), linewidth=3)
#         end

#         Label(
#             fig[fldmod1(xi, ncol)..., Top()],
#             labels[xi],
#             fontsize = 11,
#             font = :bold,
#             padding = (0, 5, 5, 0),
#             halign = :center
#         )
#     end

#     if variable ∈ ["$(GCM)_mean_dhw", "$(GCM)_initial_coral_cover", :so_to_si]
#         linkaxes!(filter(x -> x isa Axis, fig.content)...)
#     end

#     last_figure = fldmod1(length(gdf), ncol)
#     if last_figure[2] < ncol
#         Legend(fig[last_figure[1], last_figure[2]+1:ncol], legend_entries, unique(categories), nbanks=1)
#     else
#         Legend(fig[last_figure[1] + 1, 1:ncol], legend_entries, unique(categories), nbanks=1)
#     end
#     resize_to_layout!(fig)

#     display(fig)

#     return fig
# end

function bootstrap_median_ts(timeseries)
    ts_median = Vector{Union{Missing, Float64}}(missing, size(timeseries, 1))
    ts_lb = Vector{Union{Missing, Float64}}(missing, size(timeseries, 1))
    ts_ub = Vector{Union{Missing, Float64}}(missing, size(timeseries, 1))

    for t in 1:size(timeseries, 1)
        if all(.!ismissing.(timeseries[t, :]))
            bootstrap_t = bootstrap(median, timeseries[t, :], BasicSampling(1000))
            ts_median[t] = first(bootstrap_t.t0)
            ts_lb[t] = getindex(getindex(confint(bootstrap_t, NormalConfInt(0.95)), 1), 2)
            ts_ub[t] = getindex(getindex(confint(bootstrap_t, NormalConfInt(0.95)), 1), 3)
        end
    end

    return ts_median, ts_lb, ts_ub
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

# function bioregion_grouped_lagged_timeseries(
#     dataframe, bellwether_reefs_col, length_t, grouping, lag, ncol;
#     xlabel="Years",
#     ylabel="Value",
#     title=""
# )
#     dataframe = sort(dataframe, [bellwether_reefs_col, :management_area]; rev=true)
#     categories = categorical(dataframe[:, bellwether_reefs_col])

#     gdf = DataFrames.groupby(dataframe, grouping)

#     labels = label_lines.([first(df.bioregion) for df in gdf])

#     x_fig_size, y_fig_size = (2130, 1500)
#     fig = Figure(backgroundcolor = RGBf(0.98, 0.98, 0.98), size = (x_fig_size, y_fig_size))

#     colors = [Makie.wong_colors(); Makie.wong_colors()];
#     legend_entries = [
#         [PolyElement(;color = colors[1]), LineElement(;color=:blue)],
#         [PolyElement(;color = colors[2]), LineElement(;color=:red)]
#     ]

#     xsize = x_fig_size / (ncol*2)
#     n_fig_row = first(fldmod1(length(gdf), ncol))
#     ysize = y_fig_size / (n_fig_row*1.5)

#     for (xi, groupdf) in enumerate(gdf)
#         categories = categorical(groupdf[:, bellwether_reefs_col])

#         target_reefs = Matrix(
#         DataFrames.select(
#             groupdf[groupdf[:, bellwether_reefs_col] .== "bellwether", :],
#             DataFrames.Between(Symbol(length_t[1]), Symbol(length_t[2]))
#             )
#         )'
#         if lag > 0
#             target_reefs_buffer = fill(missing, 4, size(target_reefs, 2))
#             target_reefs = vcat(target_reefs_buffer, target_reefs)
#         end
#         target_reefs_median, target_reefs_lb, target_reefs_ub = bootstrap_median_ts(target_reefs)



#         non_target_reefs = Matrix(
#             DataFrames.select(
#                 groupdf[groupdf[:, bellwether_reefs_col] .== "non-bellwether", :],
#                 DataFrames.Between(Symbol(length_t[1]), Symbol(length_t[2]))
#             )
#         )'
#         non_target_reefs_median, non_target_reefs_lb, non_target_reefs_ub = bootstrap_median_ts(non_target_reefs)

#         if xi == 1
#             ax = Axis(
#                 fig[fldmod1(xi, ncol)...];
#                 xlabel = xlabel,
#                 xlabelsize = label_size+1,
#                 ylabelsize = label_size+1,
#                 ylabel = ylabel,
#                 xticklabelsize=label_size,
#                 yticklabelsize=label_size,
#                 title=title,
#                 width=xsize,
#                 height=ysize
#             )
#         else
#             ax = Axis(
#                 fig[fldmod1(xi, ncol)...];
#                 xticklabelsize=label_size,
#                 yticklabelsize=label_size,
#                 title=title,
#                 width=xsize,
#                 height=ysize
#             )
#         end

#         band_indices = first(length_t) + lag:last(length_t) + lag
#         bands = band!(band_indices, target_reefs_lb[band_indices], target_reefs_ub[band_indices]; color=(colors[2], 0.3))
#         #series!(ax, target_reefs', solid_color=(colors[2], 0.3))
#         series!(target_reefs_median', solid_color=:red)

#         band!(first(length_t):last(length_t), non_target_reefs_lb, non_target_reefs_ub; color=(colors[1], 0.3))
#         #series!(non_target_reefs', solid_color=(colors[1], 0.2))
#         series!(non_target_reefs_median', solid_color=:blue)

#         Label(
#             fig[fldmod1(xi, ncol)..., Top()],
#             labels[xi],
#             fontsize = 11,
#             font = :bold,
#             padding = (0, 5, 5, 0),
#             halign = :center
#         )
#     end

#     linkaxes!(filter(x -> x isa Axis, fig.content)...)

#     last_figure = fldmod1(length(gdf), ncol)
#     if last_figure[2] < ncol
#         Legend(fig[last_figure[1], last_figure[2]+1:ncol], legend_entries, unique(categories), nbanks=1)
#     else
#         Legend(fig[last_figure[1] + 1, 1:ncol], legend_entries, unique(categories), nbanks=1)
#     end
#     resize_to_layout!(fig)

#     display(fig)

#     return fig, gdf.keymap
# end

function figure_layouts(n_bioregions, n_col)
    return [fldmod1(x, n_col) for x in 1:n_bioregions]
end

function plot_map_bellwether_reefs(gdf::Union{DataFrame,DataFrameRow}, color_by::Symbol, bellwether_reefs_col::Symbol; geom_col::Symbol=:geometry)
    if first(unique(gdf.management_area_short)) ∈ ["FarNorthern", "Cairns-Cooktown"]
        f = Figure(; size=(700, 1000))
    else
        f = Figure(; size=(1000, 800))
    end
    ga = GeoAxis(
        f[1, 1];
        dest="+proj=latlong +datum=WGS84",
        xlabel="Longitude",
        ylabel="Latitude",
        xticklabelpad=15,
        yticklabelpad=40,
        xticklabelsize=10,
        yticklabelsize=10,
    )

    plottable = _convert_plottable(gdf, geom_col)

    # Define the unique colors and names for each level of factor color_by.
    # Use a different color palette for factors with high numbers of levels
    # (this palette is not as good for visualisation).
    if size(unique(gdf[:, color_by]),1) <= 20
        palette = ColorSchemes.tableau_20.colors
        alph_palette = colorscheme_alpha(ColorSchemes.tableau_20, 0.6)
    else
        palette = ColorSchemes.flag_ec.colors
        alph_palette = colorscheme_alpha(ColorSchemes.flag_ec, 0.6)
    end

    color_indices = groupindices(DataFrames.groupby(gdf, color_by))
    names = unique(DataFrame(indices=color_indices, names=gdf[:, color_by]))

    # Create the unique legend entries for each level of color_by
    unique_names = names.names
    legend_entries = []
    for name in eachrow(names)
        col = palette[name.indices]
        LE = PolyElement(; color=col)
        push!(legend_entries, [LE])
    end
    bellwether_indices = gdf[:, bellwether_reefs_col] .== "bellwether"
    non_bellwether_indices = gdf[:, bellwether_reefs_col] .== "non-bellwether"

    non_bell_polys = poly!(ga, plottable[non_bellwether_indices], color=alph_palette[color_indices[non_bellwether_indices]])
    bell_polys = poly!(ga, plottable[bellwether_indices], color=palette[color_indices[bellwether_indices]], strokewidth=0.5)

    Legend(f[2, 1], legend_entries, unique_names, nbanks=3, tellheight=true,
    tellwidth=false, orientation=:horizontal, labelsize=9, patchsize=(11,11))
    resize_to_layout!(f)

    display(f)
    return f, ga
end

function optimum_columns(n_bioregions)
    if n_bioregions < 5
        return 2
    elseif n_bioregions < 10
        return 3
    elseif n_bioregions < 20
        return 4
    end

    return 5
end

function ecs_plot(ecs_values, low_conf_range, high_conf_range, impact_labels)
    fig = Figure()
    ax = Axis(
        ylab = "Equilibrium Climate Sensitivity ()"
    )
end

function grouped_cluster_timeseries_plots(
    timeseries_array,
    dataframe,
    cluster_col,
    grouping,
    length_t;
    x_fig_size=1500,
    y_fig_size=1000
)
    n_col = optimum_columns(length(unique(dataframe[:, grouping])))
    fig, gdf, plot_layout = _setup_grouped_figure(
        dataframe,
        grouping;
        x_fig_size=x_fig_size,
        y_fig_size=y_fig_size
    )

    labels = label_lines.(first(df[:, grouping]) for df in gdf)

    for (xi, groupdf) in enumerate(gdf)
        plot_layout_xi = plot_layout[xi]
        groupdf = sort(groupdf, cluster_col)

        # Ensure that relative cover timeseries match the cluster allocations from groupdf
        group_timeseries = timeseries_array[length_t, timeseries_array.locations .∈ [groupdf.UNIQUE_ID]]
        group_timeseries = group_timeseries[:, indexin(groupdf.UNIQUE_ID, String.(group_timeseries.locations))]

        clusters = Int64.(groupdf[:, cluster_col])

        ADRIA.viz.clustered_scenarios!(
            fig[plot_layout_xi...],
            group_timeseries,
            clusters;
            opts = Dict{Symbol, Any}(:legend => false),
            axis_opts = Dict(:title => labels[xi])
        )
    end

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
    x_fig_size=1500,
    y_fig_size=1000
)
    n_col = optimum_columns(length(unique(dataframe[:, grouping])))
    fig, gdf, plot_layout = _setup_grouped_figure(
        dataframe,
        grouping;
        x_fig_size=x_fig_size,
        y_fig_size=y_fig_size,
        marker=PolyElement
    )

    labels = label_lines.(first(df[:, grouping]) for df in gdf)
    colors = [:green, :orange, :blue]
    xsize, ysize = _axis_size(gdf, x_fig_size, y_fig_size, n_col)

    for (xi, groupdf) in enumerate(gdf)
        #println("$(xi)")
        plot_layout_xi = plot_layout[xi]
        groupdf = sort(groupdf, cluster_col)

        clusters = Int64.(groupdf[:, cluster_col])     

        # bellwether_var_med = median(groupdf[groupdf[:, bellwether_reefs_col] .== "bellwether", variable])
        # non_bellwether_var_75 = quantile(groupdf[groupdf[:, bellwether_reefs_col] .== "non-bellwether", variable], 0.75)
        # non_bellwether_var_25 = quantile(groupdf[groupdf[:, bellwether_reefs_col] .== "non-bellwether", variable], 0.25)

        # if bellwether_var_med < non_bellwether_var_25
        #     background_color = (:royalblue, 0.2)
        # elseif bellwether_var_med > non_bellwether_var_75
        #     background_color = (:goldenrod1, 0.2)
        # else
        #     background_color = :white
        # end

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
            hlines!(1; color=(:gray, 0.5), linewidth=3)
        end

        Label(
            fig[plot_layout_xi..., Top()],
            labels[xi],
            fontsize = 16,
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
        fontsize=18
    )

    #linkaxes!(filter(x -> x isa Axis, fig.content)...)
    resize_to_layout!(fig)

    display(fig)

    return fig
end

function cluster_analysis_plots(analysis_layers, rel_cover, dhw_ts, grouping, fig_out_dir)

    mean_dhw_violin = grouped_cluster_violin_plots(
        analysis_layers,
        Symbol("$(GCM)_$(grouping)_clusters"),
        grouping, Symbol("$(GCM)_mean_dhw");
        ylabel="mean DHW", xlabel="Clusters"
    );
    save(
        joinpath(fig_out_dir, "$(grouping)", "mean_dhw_$(grouping)_violin.png"), 
        mean_dhw_violin
    )

    so_to_si_violin = grouped_cluster_violin_plots(
        analysis_layers,
        Symbol("$(GCM)_$(grouping)_clusters"),
        grouping, :so_to_si;
        ylabel="Source to sink ratio", xlabel="Clusters"
    );
    save(
        joinpath(fig_out_dir, "$(grouping)", "so_to_si_$(grouping)_violin.png"), 
        so_to_si_violin
    )

    total_strength_violin = grouped_cluster_violin_plots(
        analysis_layers,
        Symbol("$(GCM)_$(grouping)_clusters"),
        grouping, :total_strength;
        ylabel="Total connectivity strength", xlabel="Clusters"
    );
    save(
        joinpath(fig_out_dir, "$(grouping)", "total_strength_$(grouping)_violin.png"), 
        total_strength_violin
    )

    initial_cover_violin = grouped_cluster_violin_plots(
        analysis_layers,
        Symbol("$(GCM)_$(grouping)_clusters"),
        grouping, Symbol("initial_coral_cover");
        ylabel="Initial coral cover", xlabel="Clusters"
    );
    save(
        joinpath(fig_out_dir, "$(grouping)", "initial_cover_$(grouping)_violin.png"), 
        initial_cover_violin
    )

    initial_proportion_violin = grouped_cluster_violin_plots(
        analysis_layers,
        Symbol("$(GCM)_$(grouping)_clusters"),
        grouping, Symbol("initial_proportion");
        ylabel="Initial proportion of habitable area occupied", xlabel="Clusters"
    );
    save(
        joinpath(fig_out_dir, "$(grouping)", "initial_proportion_$(grouping)_violin.png"), 
        initial_proportion_violin
    )

    dhw_cover_cor_violin = grouped_cluster_violin_plots(
        analysis_layers,
        Symbol("$(GCM)_$(grouping)_clusters"),
        grouping, Symbol("$(GCM)_dhw_cover_cor");
        ylabel="Total coral cover - DHW correlation", xlabel="Clusters"
    );
    save(
        joinpath(fig_out_dir, "$(grouping)", "dhw_cover_cor_$(grouping)_violin.png"), 
        dhw_cover_cor_violin
    )

    analysis_layers_depth = analysis_layers[analysis_layers.depth_mean .!= 7, :]
    depth_median_violin = grouped_cluster_violin_plots(
        analysis_layers_depth,
        Symbol("$(GCM)_$(grouping)_clusters"),
        grouping, Symbol("depth_med");
        ylabel="Median Depth (m)", xlabel="Clusters"
    );
    save(
        joinpath(fig_out_dir, "$(grouping)", "depth_$(grouping)_violin.png"), 
        depth_median_violin
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
        bioregion_grouped_timeseries
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
        bior_reefs_dhw_plot
    )

    return nothing
end