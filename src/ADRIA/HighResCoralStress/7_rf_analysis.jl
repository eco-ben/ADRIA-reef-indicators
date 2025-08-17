using CSV
using CategoricalArrays
using MLJ, MLJDecisionTreeInterface, StatsBase
using Distributed

include("../../common.jl")

context_layers = GDF.read(joinpath(output_path, "analysis_context_layers_carbonate.gpkg"))
context_layers.log_so_to_si = log10.(context_layers.so_to_si)
context_layers.log_total_strength = log10.(context_layers.total_strength)
context_layers.log_depth_med = log10.(context_layers.depth_med)

dhw_scenarios = open_dataset(joinpath(gbr_domain_path, "DHWs/dhwRCP45.nc"))
GCMs = dhw_scenarios.dhw.properties["members"]

# Prepare long-form of data to support multinomial analysis.
# Each row should be unique attribute combination of reef properties:
# reef, depth, connectivity, dhw, cluster and GCM
gcm_dhw_cols = [Symbol("$(GCM)_mean_dhw") for GCM in GCMs]
gcm_bioregion_cluster_cols = [Symbol("$(GCM)_bioregion_clusters") for GCM in GCMs]

collated_reef_properties = Vector{DataFrame}(undef, length(GCMs))
target_cols = [:UNIQUE_ID, :GBRMPA_ID, :depth_med, :total_strength, :so_to_si, :bioregion, :bioregion_average_latitude]
reef_properties = context_layers[:, target_cols]
reef_properties.abs_k_area = context_layers.area .* context_layers.k
for (i, GCM) in enumerate(GCMs)
    dhw_col = Symbol("$(GCM)_mean_dhw")
    bio_cluster_col = Symbol("$(GCM)_bioregion_clusters")

    bio_cluster_details = context_layers[:, [dhw_col, bio_cluster_col]]
    rename!(bio_cluster_details, dhw_col => :mean_dhw, bio_cluster_col => :cluster)
    bio_cluster_details.GCM .= GCM
    bio_cluster_details.cluster .= categorical(bio_cluster_details.cluster)

    collated_reef_properties[i] = hcat(
        copy(reef_properties),
        bio_cluster_details
    )
end

reef_properties = vcat(collated_reef_properties...)
reef_properties.GCM .= categorical(reef_properties.GCM)
reef_properties.bioregion .= categorical(reef_properties.bioregion)
reef_properties.abs_k_area ./= 1e6  # Convert to km^2 for clarity

target_cols = [:depth_med, :total_strength, :so_to_si, :bioregion, :mean_dhw, :abs_k_area]
Xs = reef_properties[:, target_cols]
y = vec(reef_properties[:, :cluster])
y = Int64.(getfield.(y, :ref))

Xs.bioregion .= categorical(Xs.bioregion)
bioregion_cats = Xs.bioregion
Xs.bioregion .= Int64.(getfield.(Xs.bioregion, :ref))

d = Dict(1 => "low", 2 => "med", 3 => "high")

X2 = deepcopy(Xs)
y2 = categorical(map(x -> d[x], y), levels=["low", "med", "high"], ordered=true)

# RandomForestClassifer expects only OrderedFactor (instead of multiclass/categoricals)
# So we have to coerce categoricals to this.
# https://stackoverflow.com/a/78139992
X2 = coerce(
    X2,
    :depth_med => Continuous,
    :total_strength => Continuous,
    :so_to_si => Continuous,
    :bioregion => OrderedFactor,
    :mean_dhw => Continuous,
    :abs_k_area => Continuous,
)

readable_names = OrderedDict(
    "mean_dhw" => "Mean DHW",
    "total_strength" => "Total Connectivity Strength",
    "depth_med" => "Median Depth",
    "bioregion" => "Bioregion",
    "abs_k_area" => "Carrying Capacity",
    "so_to_si" => "Outgoing to Incoming\nConnectivity Ratio",
)

# General (naive) overview
sp_cors = Dict(
    "mean_dhw" => round(corspearman(X2.mean_dhw, y); digits=2),
    "total_strength" => round(corspearman(X2.total_strength, y), digits=2),
    "depth_med" => round(corspearman(X2.depth_med, y); digits=2),
    "bioregion" => round(corspearman(levelcode.(X2.bioregion), y); digits=2),
    "abs_k_area" => round(corspearman(X2.abs_k_area, y); digits=2),
    "so_to_si" => round(corspearman(X2.so_to_si, y); digits=2),
)

rng = Random.seed!(76)
rf = @load RandomForestClassifier pkg = "DecisionTree"
model = rf(; rng=rng)

train_size = ceil(Int64, nrow(X2) * 0.6)
shuffle_set = shuffle(1:nrow(X2))

train_X = X2[shuffle_set[1:train_size], :]
test_X = X2[shuffle_set[train_size+1:end], :]
train_y = y2[shuffle_set[1:train_size]]
test_y = y2[shuffle_set[train_size+1:end]]

mach = MLJ.machine(model, train_X, train_y)
fit!(mach)

cluster_est = MLJ.predict(mach, train_X)
y_pred_mode = predict_mode(mach, train_X)
μ_accuracy = mean(y_pred_mode .== train_y)
@info μ_accuracy  # 0.99

cluster_est = MLJ.predict(mach, test_X)
y_pred_mode = predict_mode(mach, test_X)

μ_accuracy = mean(y_pred_mode .== test_y)
@info μ_accuracy  # 0.66

cm = confusion_matrix(y_pred_mode, test_y)
cm.mat

class_names = ["Low", "Med", "High"]
for i in 1:3
    class_name = class_names[i]

    # True Positives (TP) - the diagonal
    tp = cm.mat[i, i]

    # False Positives (FP)
    # Predicted as class i but actually other classes
    fp = sum(cm.mat[i, :]) - tp

    # False Negatives (FN)
    # Actually class i, but predicted as others
    fn = sum(cm.mat[:, i]) - tp

    # True Negatives (TN): everything else
    total_samples = sum(cm.mat)
    tn = total_samples - tp - fp - fn

    precision = tp / (tp + fp)
    recall = tp / (tp + fn)
    f1 = 2 * (precision * recall) / (precision + recall)

    println("""
    Class: $class_name
        TP (True Positives): $tp
        FP (False Positives): $fp
        FN (False Negatives): $fn
        TN (True Negatives): $tn
        Check: TP+FP+FN+TN = $(tp+fp+fn+tn), Total = $total_samples
        Precision: $(round(precision; digits=3))
        Recall: $(round(recall; digits=3))
        F1: $(round(f1; digits=3))
    """)
end


cv = CV(nfolds=3, rng=123)
evaluate!(mach, resampling=cv, measure=[MLJ.accuracy, MLJ.cross_entropy])
# Accuracy: 0.669
# LogLoss: 0.751

feat_imports = feature_importances(mach)
#       :depth_med => 0.2516920696138505
#        :mean_dhw => 0.24192024951943078
#      :abs_k_area => 0.14411108965272595
#  :total_strength => 0.13292360681207943
#        :so_to_si => 0.12930894620395628
#       :bioregion => 0.10004403819795683
# The precise values of the above are subject to change but the order of features is
# robust.

fi_dict = OrderedDict(feat_imports)
fi_names = string.(reverse(collect(keys(fi_dict))))
fi_vals = reverse(collect(values(fi_dict)))

# Reorder columns based on feature importance
X2 = X2[:, string.([first(fi) for fi in feat_imports])]

# Plot overview based on Spearman correlation
f = Figure(size=(800, 500))
r = 1
c = 1
n_rows = 2
n_cols = 3
for t_col in names(X2)
    d = X2[:, t_col]

    if d isa CategoricalArray
        d_order = sortperm(levelcode.(d))
        v = getfield.(d[d_order], :ref)
    else
        d_order = sortperm(d)
        v = d[d_order]
    end

    t_str = string(t_col)
    ax = Axis(f[r, c])
    scatter!(ax, v; color=y[d_order], alpha=0.2, markersize=6)
    ax.title = "$(readable_names[t_str])\n($(sp_cors[t_str]))"

    c += 1
    if c > n_cols
        c = 1
        r += 1
    end
end

Label(f[:, 0], "Feature value", rotation=π / 2)
Label(f[3, :], "Number of samples")
Legend(
    f[4, :],
    [PolyElement(color=col, alpha=0.8) for col in [:green, :orange, :blue]],
    ["Low", "Medium", "High"],
    "Cluster",
    nbanks=3
)
rowsize!(f.layout, 4, Auto(0.5))

save(
    joinpath(figs_path, "ts_cluster_feature_sp_corr.png"),
    f,
    px_per_unit=dpi
)


# Plot partial dependencies
function partial_dependence_multiclass(machine, data, feature_col; n_grid=50)
    feature_values = data[:, feature_col]

    is_cat_type = eltype(feature_values) <: Union{String,Symbol,CategoricalValue,AbstractString}
    if is_cat_type || isa(feature_values, CategoricalArray)
        # Skip categorical (uninformative)
        return zeros(n_grid)
        # grid = unique(feature_values)
    else
        grid = range(extrema(feature_values)..., length=n_grid)
    end

    # Get class labels from a sample prediction
    sample_pred = MLJ.predict(machine, data[1:1, :])
    class_labels = levels(sample_pred[1])

    # Initialize results for each class
    results = OrderedDict()
    for class_label in class_labels
        results[class_label] = Float64[]
    end

    data_modified = deepcopy(data)
    for grid_val in grid

        # Make feature constant
        data_modified[:, feature_col] .= grid_val

        # Prediction if feature is fixed
        prob_predictions = MLJ.predict(machine, data_modified)

        # Calculate average probability for each class
        for class_label in class_labels
            class_probs = pdf.(prob_predictions, Ref(class_label))
            push!(results[class_label], mean(class_probs))
        end
    end

    return (grid=collect(grid), class_values=results)
end

function plot_multiclass_pd(pdp_results::OrderedDict, matching_reef_df::DataFrame; fig=Figure(size=(800,500)))
    feature_names = keys(pdp_results)
    readable_names = Dict(
        "mean_dhw" => "Mean DHW",
        "total_strength" => "Total Connectivity Strength",
        "depth_med" => "Median Depth",
        "bioregion" => "Bioregion",
        "abs_k_area" => "Carrying Capacity",
        "so_to_si" => "Outgoing to Incoming\nConnectivity Ratio",
    )

    n_rows = 2
    n_cols = 3

    # Color palette for classes (low, med, high)
    colors = [(:green, 0.4), (:orange, 0.4), (:blue, 0.4)]

    r = 1
    c = 1
    for fn in feature_names
        if fn == "bioregion"
            # Skip bioregion - uninformative
            continue
        end

        grid, class_values = pdp_results[fn]
        ax = Axis(
            fig[r, c],
            xlabel=readable_names[fn],
            limits=(nothing, (0, 1))
        )

        class_labels = collect(keys(class_values))

        if eltype(grid) <: Union{Int,String,Symbol,CategoricalValue,AbstractString}
            # Categorical - grouped bar plot
            x_positions = sort(unique(grid))  # 1:length(grid)
            n_classes = length(class_labels)
            bar_width = 0.8 / n_classes

            for (i, class_label) in enumerate(class_labels)
                x_offset = x_positions .+ (i - (n_classes + 1) / 2) * bar_width
                barplot!(ax, x_offset, class_values[class_label],
                    width=bar_width,
                    color=colors[i],
                    alpha=0.7,
                    label=string(class_label))
            end

            ax.xticks = (x_positions, string.(grid))
            ax.xticklabelrotation = π / 4
        else
            # Continuous - multiple line plots
            for (i, class_label) in enumerate(class_labels)
                lines!(
                    ax, grid, class_values[class_label],
                    color=colors[i],
                    linewidth=3,
                    label=string(class_label)
                )
                scatter!(
                    ax, grid, class_values[class_label],
                    color=colors[i],
                    markersize=6,
                    alpha=0.7
                )
            end
            var_values = matching_reef_df[:, fn]
            y_dist_scatter = rand(length(var_values)) .* (0.9 - 0.8) .+ 0.8
            scatter!(
                ax, var_values, y_dist_scatter, color=(:gray, 0.3), markersize=6
            )
        end

        c += 1

        if c == n_cols + 1
            r += 1
            c = 1
        end
    end

    Legend(
        fig[n_rows+1, :],
        [PolyElement(color=col, alpha=0.8) for col in colors],
        ["Low", "Medium", "High"],
        orientation=:horizontal
    )

    # if isa(fig, GridLayout)
    #     rowsize!(fig, n_rows+1, Relative(0.15))
    #     for j in 1:3
    #         colsize!(fig, j, Relative(1/3))
    #     end
    # else
    #     rowsize!(fig.layout, n_rows+1, Relative(0.15))
    #     for j in 1:3
    #         colsize!(fig.layout, j, Relative(1/3))
    #     end
    # end

    # Label(fig[1:2, -1], "Probability", rotation=π / 2)
    Label(fig[0, 1:3], "Partial Dependence", font=:bold)

    return fig
end

pdp_vals = OrderedDict()
for n in names(X2)
    pdp_vals[n] = partial_dependence_multiclass(mach, test_X, n; n_grid=50)
end

# Plot feature importances
figure = Figure(
    size = (fig_sizes["cluster_hm_width"], fig_sizes["cluster_hm_height"] - 3centimetre), 
    fontsize=fontsize
)
gr1 = GridLayout(figure[1,1])
ax = Axis(
    gr1[1,1], 
    yticks=(1:6, map(x -> readable_names[x], fi_names)),
    title="Feature Importance"
)
bar = barplot!(
    ax,
    fi_vals,
    direction=:x,
)

gr2 = GridLayout(figure[2,1])
plot_multiclass_pd(pdp_vals, test_X; fig=gr2)
rowsize!(figure.layout, 1, Relative(0.25))
Label(figure[2,1, Left()], "Probability", rotation=π/2, tellwidth=false)
Label(figure[1,1, TopLeft()], "A", font=:bold, tellwidth=false, fontsize=fontsize+3)
Label(figure[2,1, TopLeft()], "B", font=:bold, tellwidth=false, fontsize=fontsize+3)

save(
    joinpath(figs_path, "ts_cluster_rf_feature_analysis.png"),
    figure,
    px_per_unit=dpi
)

### The below notes are from a previous version. This version was ammended to ensure consistent fontsize/dpi.
# Above two figures are manually joined together and given panel labels A and B.
# The file is saved as "ts_cluster_rf_feature_analysis.png" in the Figure directory.
# This manual two-part figure can be recreated by using the old code for panel A and the following code 
# for panel B (then adding the Probability label):
# f2 = plot_multiclass_pd(pdp_vals, test_X)


test_X.pred_y = y_pred_mode
test_X.y = test_y
test_X.bioregion = bioregion_cats[shuffle_set[train_size+1:end]]
test_X.bioregion_average_lat = reef_properties.bioregion_average_latitude[shuffle_set[train_size+1:end]]

perf_by_bio = combine(groupby(test_X, :bioregion)) do sdf
    (; n = nrow(sdf),
       acc = mean(sdf.y .== sdf.pred_y),
       size = mean(sdf.abs_k_area),
       depth = mean(sdf.depth_med),
       dhw = mean(sdf.mean_dhw),
       conn_str = mean(sdf.total_strength),
       so_to_si = mean(sdf.so_to_si),
       bior_ave_lat = mean(sdf.bioregion_average_lat)
    )
end
perf_by_bio = sort(perf_by_bio, :bior_ave_lat, rev=true)

bioregion_colors = distinguishable_colors(nrow(perf_by_bio))

bgcol=:gray90
fig = Figure(
    size = (fig_sizes["cluster_hm_width"], fig_sizes["cluster_hm_height"] + 3centimetre), 
    fontsize=fontsize,
    backgroundcolor=bgcol
)
ax = Axis(
    fig[1,1],
    ylabel = "Mean bioregion accuracy",
    xlabel = "Mean reef depth [m]",
    backgroundcolor=bgcol
)
Makie.scatter!(ax, perf_by_bio.depth, perf_by_bio.acc, color=bioregion_colors)
ax = Axis(
    fig[1,2],
    ylabel = "Mean bioregion accuracy",
    xlabel = "Mean carrying capacity [km²]",
    backgroundcolor=bgcol
)
Makie.scatter!(ax, perf_by_bio.size, perf_by_bio.acc, color=bioregion_colors)
ax = Axis(
    fig[2,1],
    ylabel = "Bioregion accuracy",
    xlabel = "Mean DHW [\u00B0C - Weeks]",
    backgroundcolor=bgcol
)
Makie.scatter!(ax, perf_by_bio.dhw, perf_by_bio.acc, color=bioregion_colors)
ax = Axis(
    fig[2,2],
    ylabel = "Bioregion accuracy",
    xlabel = "Mean total connectivity strength",
    backgroundcolor=bgcol
)
Makie.scatter!(ax, perf_by_bio.n, perf_by_bio.acc, color=bioregion_colors)
ax = Axis(
    fig[3,1],
    ylabel = "Bioregion accuracy",
    xlabel = "Mean source-to-sink ratio",
    backgroundcolor=bgcol
)
Makie.scatter!(ax, perf_by_bio.so_to_si, perf_by_bio.acc, color=bioregion_colors)
ax = Axis(
    fig[3,2],
    ylabel = "Bioregion accuracy",
    xlabel = "Number of reefs per bioregion",
    backgroundcolor=bgcol
)
Makie.scatter!(ax, perf_by_bio.n, perf_by_bio.acc, color=bioregion_colors)

Legend(
    fig[1:3, 0],
    [MarkerElement(; color=bio_col, marker=:circle) for bio_col in bioregion_colors],
    label_lines.(string.(perf_by_bio.bioregion); l_length=12),
    orientation=:vertical,
    nbanks=1,
    backgroundcolor=bgcol,
    framewidth=0.3
)

save(
    joinpath(figs_path, "rf_performance_across_predictors.png"),
    fig,
    px_per_unit=dpi
)