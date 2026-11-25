# ── Knot-type distribution ─────────────────────────────────────────────────

using StatsBase
using CairoMakie

"""
    knot_distribution([settings])

Plot a bar chart of the knot-type distribution (relative frequencies) with
rich-text labels in Alexander–Briggs notation.  Every catalog row is counted,
including unidentified knot types, so the bars sum to the full catalog.
"""
function knot_distribution(settings::Settings = DEFAULT_SETTINGS)
    data, _ = load_data_and_group_knots(settings)
    knot_types = data[:, 1]

    freq = countmap(knot_types)
    labels = collect(keys(freq))
    counts = collect(values(freq))
    order = sortperm(counts, rev = true)
    labels, counts = labels[order], counts[order]

    fig = Figure(size = (600, 380))
    x = 1:2:2 * length(labels)

    ax = Axis(fig[1, 1],
        xlabel = "Knot Type (Alexander–Briggs Notation)",
        ylabel = "Relative Share [%]",
        xticklabelrotation = deg2rad(45),
        xticks = (x, beautify_knot_labels(labels)),
        xticklabelsize = 15.0,
    )

    colors = cgrad(:managua10, 10, categorical = true)
    shares = counts ./ sum(counts) * 100
    barplot!(ax, x, shares,
        color = (colors[end-2], 0.4),
        strokecolor = colors[end-2], strokewidth = 1)
    text!(ax, x, shares;
        text = string.(counts), align = (:center, :bottom), fontsize = 13)
    # Leave headroom above the tallest bar so its count label is not clipped.
    ylims!(ax, 0, maximum(shares) * 1.15)

    outpath = joinpath(settings.figures_dir, "knot_distribution.pdf")
    mkpath(dirname(outpath))
    save(outpath, fig)
    println("Plot saved to $outpath")
    return fig
end
