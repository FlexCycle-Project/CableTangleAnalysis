# ── PCA projection analysis ────────────────────────────────────────────────
#
# Runs PCA on the cleared catalog and produces 2D scatter plots of the
# principal components, coloured by knot type and sized by threading count.

using MultivariateStats
using DelimitedFiles
using Statistics

"""
    write_pca_results(model, data, header, metrics, scores, settings)

Write the PCA loadings/explained-variance table and the per-row PCA score
table to `settings.analysis_dir`.
"""
function write_pca_results(model, data, header, metrics, scores, settings::Settings)
    loadings = projection(model)
    variances = principalvars(model)
    total_variance = tvar(model)
    explained = total_variance > 0.0 ? variances ./ total_variance : zeros(length(variances))
    cumulative = cumsum(explained)
    ncomp = size(loadings, 2)

    mkpath(settings.analysis_dir)

    result_header = vcat(
        "component", "eigenvalue", "explained_variance",
        "explained_variance_pct", "cumulative_pct", metrics,
    )
    result_rows = Vector{Vector{String}}()
    for component in 1:ncomp
        row = String[
            "PC$(component)",
            string(round(variances[component], sigdigits = 8)),
            string(round(explained[component], sigdigits = 8)),
            string(round(explained[component] * 100.0, sigdigits = 6)),
            string(round(cumulative[component] * 100.0, sigdigits = 6)),
        ]
        append!(row, string.(round.(loadings[:, component], sigdigits = 8)))
        push!(result_rows, row)
    end

    result_path = joinpath(settings.analysis_dir, "pca_results.tsv")
    open(result_path, "w") do io
        println(io, join(result_header, '\t'))
        for row in result_rows
            println(io, join(row, '\t'))
        end
    end
    println("PCA results table saved to $result_path")

    header_names = string.(vec(header))
    seed_index = findfirst(==("seed"), header_names)
    threading_index = findfirst(==("totalThreadingCount"), header_names)
    knot_index = something(findfirst(==("knotType"), header_names), 1)

    score_header = vcat(
        "seed", "knotType", "totalThreadingCount",
        ["PC$(component)" for component in 1:ncomp]...,
    )
    score_path = joinpath(settings.analysis_dir, "pca_scores.tsv")
    open(score_path, "w") do io
        println(io, join(score_header, '\t'))
        for observation in 1:size(data, 1)
            seed = seed_index === nothing ? "" : string(data[observation, seed_index])
            threading = threading_index === nothing ? "" :
                string(data[observation, threading_index])
            row = String[seed, string(data[observation, knot_index]), threading]
            append!(row, string.(round.(scores[:, observation], sigdigits = 8)))
            println(io, join(row, '\t'))
        end
    end
    println("PCA score table saved to $score_path")
    return nothing
end

"""
    run_pca([settings::Settings]; datafile, outdir, metrics)

Load the cleared metrics file, run PCA on the selected metrics, and produce a
2D combined figure (PC1 vs PC2, PC1 vs PC3, PC2 vs PC3).

Returns the fitted PCA model and the projected data.
"""
function run_pca(settings::Settings = DEFAULT_SETTINGS;
        datafile = joinpath(settings.data_dir, "cleared_metrics_2.tsv"),
        outdir = joinpath(settings.figures_dir, "pca"),
        metrics = COMPLEXITY_FEATURES,
    )
    data, header = readdlm(datafile, '\t', header = true)

    metrics_range = Int[]
    for metric in metrics
        index = findfirst(==(metric), vec(header))
        index === nothing || push!(metrics_range, index)
    end
    @assert !isempty(metrics_range) "None of the requested metrics were found in the header."

    println("Using Metrics: ", header[metrics_range])

    # Standardise (z-score).  Center first, then scale by the standard
    # deviation of the centered columns (identical to the spread of the raw
    # columns, but keeps the historical arithmetic order).
    raw = data[:, metrics_range]'
    centered = raw .- mean(raw, dims = 2)
    x = centered ./ std(centered, dims = 2)

    model = MultivariateStats.fit(PCA, x, maxoutdim = 3, pratio = 0.999, mean = 0)
    pcavars = principalvars(model)
    explained = pcavars ./ tvar(model) .* 100   # % variance explained per PC
    scores = predict(model, x)

    metrics_found = string.(vec(header)[metrics_range])
    write_pca_results(model, data, header, metrics_found, scores, settings)

    # ── Knot-type colour mapping ──
    label_knot_type = string.(data[:, 1])
    unique_labels = [
        "Trivial Knot", "3_i", "4_i", "5_i", "6_i",
        "7_i", "8_i", "9_i", "10_i", "11a_i", "11n_i",
    ]
    color_map = Dict(label => index for (index, label) in enumerate(unique_labels))
    color_range = (minimum(values(color_map)), maximum(values(color_map)))

    # Threading count → marker size.
    threading = float.(data[:, findfirst(==("totalThreadingCount"), vec(header))])
    println("Mean NT: ", mean(threading))
    println("Var NT: ", var(threading))
    size_scale(value) = (value - mean(threading)) / var(threading) * 100 + 20

    colors = [get(color_map, label, 0) for label in label_knot_type]
    sizes = size_scale.(threading)

    mkpath(outdir)

    # Small datasets can yield fewer than three components (three samples
    # span at most two dimensions after centering).  Only draw the panels
    # that the fitted model actually supports.
    ncomp = size(scores, 1)
    panels = Tuple{Int, Int, String, Int}[]
    ncomp >= 2 && push!(panels, (1, 2, "PC1 vs PC2", 1))
    ncomp >= 3 && push!(panels, (1, 3, "PC1 vs PC3", 2))
    ncomp >= 3 && push!(panels, (2, 3, "PC2 vs PC3", 3))
    if isempty(panels)
        @warn "PCA produced $ncomp component(s); skipping 2D projection figure"
        return model, scores
    end

    fig = Figure(size = (2000, 700), fontsize = 25)
    axes = []
    for (x_component, y_component, title, column) in panels
        ax = Axis(fig[1, column];
            xlabel = "PC$(x_component) ($(round(explained[x_component], digits=1))%)",
            ylabel = "PC$(y_component) ($(round(explained[y_component], digits=1))%)",
            title = title,
            aspect = AxisAspect(1),
        )
        scatter!(ax, scores[x_component, :], scores[y_component, :];
            color = colors, colorrange = color_range, colormap = :managua10,
            markersize = sizes, alpha = 0.85,
            strokewidth = 0.4, strokecolor = (:black, 0.25))
        push!(axes, ax)
    end

    add_pca_legends!(fig, axes[end], unique_labels, color_map, color_range,
                     threading, size_scale, col = 3)

    save(joinpath(outdir, "pca2d_combined.pdf"), fig)
    return model, scores
end

"""
    add_pca_legends!(fig, ax, unique_labels, color_map, color_range,
                     threading, size_scale; col=2)

Add knot-type colour legend and threading-size legend to a PCA figure.
"""
function add_pca_legends!(fig, ax, unique_labels, color_map, color_range,
                          threading, size_scale; col = 2)
    legend_scatter = []
    legend_labels_color = []
    for label in unique_labels
        push!(legend_scatter, scatter!(ax, [NaN], [NaN];
            color = color_map[label], colorrange = color_range,
            colormap = :managua10, markersize = 25))
        push!(legend_labels_color, beautify_knot_labels([label])[1])
    end
    Legend(fig[2, 1:col-1], legend_scatter, legend_labels_color,
           "Knot Type", orientation = :horizontal)

    legend_msize = []
    legend_labels_size = []
    for t in round.(Int, LinRange(minimum(threading), maximum(threading), 5))
        push!(legend_msize, scatter!(ax, [NaN], [NaN];
            color = :lightgray, markersize = size_scale(t)))
        push!(legend_labels_size, rich("N", subscript("t"), ": $t"))
    end
    Legend(fig[2, col], legend_msize, legend_labels_size,
           "Threading Number", orientation = :horizontal)
end
