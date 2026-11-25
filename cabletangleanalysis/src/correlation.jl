# ── Correlation analysis ───────────────────────────────────────────────────
#
# Pairwise metric correlation with FDR-corrected p-values and a heatmap.
# Spearman's rank correlation is used because the metrics are not normally
# distributed: both columns are rank-transformed (average ranks for ties) and
# `HypothesisTests.CorrelationTest` is applied to the ranks, which is the
# standard t-test approximation for the Spearman coefficient.

using HypothesisTests
using MultipleTesting
using StatsBase: tiedrank
using CairoMakie

"""
    correlation([settings::Settings])

Correlation pipeline:
1. Load data and select metrics from `settings.correlation_metrics`.
2. Compute the Spearman rank correlation matrix with FDR-corrected p-values.
3. Print significant correlations and redundant pairs (|ρ| > 0.8).
4. Save the heatmap to `settings.figures_dir/correlation_heatmap.pdf`.

Returns `(corrMat, pMatCorrected, metricNames, fig)`.
"""
function correlation(settings::Settings = DEFAULT_SETTINGS)
    data, header = load_data(settings; defined_only = true)

    indices = [findfirst(==(name), vec(header)) for name in settings.correlation_metrics]
    filter!(!isnothing, indices)
    metric_data = float.(data[:, indices])
    metric_names = string.(header[indices])
    n = size(metric_data, 2)

    # Rank-transform once per metric; the diagonal stays at 1.0.
    ranked_data = tiedrank.(eachcol(metric_data))

    corr_matrix = zeros(n, n)
    p_values = zeros(n, n)
    for i in 1:n
        corr_matrix[i, i] = 1.0
        for j in (i + 1):n
            test = CorrelationTest(ranked_data[i], ranked_data[j])
            corr_matrix[i, j] = corr_matrix[j, i] = test.r
            p_values[i, j] = p_values[j, i] = pvalue(test)
        end
    end

    # Benjamini-Hochberg correction over the unique off-diagonal pairs.
    pair_pvalues = [p_values[i, j] for i in 1:n for j in (i + 1):n]
    corrected = isempty(pair_pvalues) ? Float64[] :
        adjust(pair_pvalues, BenjaminiHochberg())

    corrected_matrix = zeros(n, n)
    index = 1
    for i in 1:n, j in (i + 1):n
        corrected_matrix[i, j] = corrected_matrix[j, i] = corrected[index]
        index += 1
    end

    println("=== Spearman Correlation Matrix ===")
    display(round.(corr_matrix, digits = 3))
    println("\n=== P-Value Matrix (FDR-corrected) ===")
    display(round.(corrected_matrix, digits = 8))

    println("\n=== Significant Correlations (FDR-corrected p < 0.05) ===")
    for i in 1:n, j in (i + 1):n
        corrected_matrix[i, j] < 0.05 || continue
        println("$(metric_names[i]) ↔ $(metric_names[j]): " *
                "ρ = $(round(corr_matrix[i, j], digits=3)), " *
                "p = $(round(corrected_matrix[i, j], digits=4)) " *
                "[$(correlation_strength(corr_matrix[i, j]))]")
    end

    fig = plot_correlation_heatmap(corr_matrix, metric_names)
    outpath = joinpath(settings.figures_dir, "correlation_heatmap.pdf")
    mkpath(dirname(outpath))
    save(outpath, fig)

    find_redundant_metrics(corr_matrix, metric_names, 0.8)

    return corr_matrix, corrected_matrix, metric_names, fig
end

"""
    correlation_strength(r) -> String

Categorise a correlation coefficient by strength.
"""
function correlation_strength(r::Float64)
    absolute = abs(r)
    absolute < 0.2 && return "Very Weak"
    absolute < 0.4 && return "Weak"
    absolute < 0.6 && return "Moderate"
    absolute < 0.8 && return "Strong"
    return "Very Strong"
end

"""
    plot_correlation_heatmap(corr_matrix, metric_names)

Upper-triangle heatmap of the correlation matrix with text annotations.
"""
function plot_correlation_heatmap(corr_matrix::Matrix{Float64}, metric_names)
    n = size(corr_matrix, 1)

    fig = Figure(size = (900, 850), fontsize = 18)
    ax = Axis(fig[1, 1],
        yticks = (1:n, metric_names),
        xticks = (1:n, metric_names),
        xlabel = "Metrics", ylabel = "Metrics",
        title = "Cable Tangle Metrics — Spearman Correlation",
        xticklabelrotation = π / 4,
    )

    # Mask the lower triangle so the heatmap shows each pair once.
    plot_matrix = copy(corr_matrix)
    for i in 1:n, j in 1:n
        i > j && (plot_matrix[i, j] = NaN)
    end

    hm = heatmap!(ax, 1:n, 1:n, plot_matrix,
        colormap = :RdBu, colorrange = (-1.0, 1.0), nan_color = :white)
    Colorbar(fig[1, 2], hm, label = "Spearman rho", width = 20)

    for i in 1:n, j in 1:n
        i <= j || continue
        value = corr_matrix[i, j]
        text!(ax, i, j;
            text = "$(round(value, digits=2))",
            fontsize = 22, font = :bold, align = (:center, :center),
            color = abs(value) > 0.5 ? :white : :black)
    end

    colsize!(fig.layout, 1, Aspect(1, 1.0))
    colgap!(fig.layout, 1, 10)
    resize_to_layout!(fig)
    return fig
end

"""
    find_redundant_metrics(corrMat, metricNames, threshold=0.8)

Print metric pairs with |ρ| > threshold (suggesting redundancy).
"""
function find_redundant_metrics(corr_matrix::Matrix{Float64}, metric_names,
                                threshold::Float64 = 0.8)
    println("\n=== Redundant Metric Pairs (|ρ| > $threshold) ===")
    redundant = Tuple{String, String, Float64}[]

    n = size(corr_matrix, 1)
    for i in 1:n, j in (i + 1):n
        abs(corr_matrix[i, j]) > threshold || continue
        push!(redundant, (metric_names[i], metric_names[j], corr_matrix[i, j]))
        println("$(metric_names[i]) ↔ $(metric_names[j]): " *
                "ρ = $(round(corr_matrix[i, j], digits=3))")
    end

    if isempty(redundant)
        println("No highly redundant metric pairs found.")
    else
        println("\nRecommendation: Consider removing one metric from each pair above.")
    end
    return redundant
end
