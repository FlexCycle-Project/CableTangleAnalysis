# ── Cut-evaluation analysis ────────────────────────────────────────────────
#
# These functions analyse how many cuts are needed to resolve cable tangles
# to zero and how many threadings a single cut removes.  The readers accept
# both the old two-strategy files and the newer multi-strategy files produced
# by CableSim: missing strategy columns are ignored.

using Statistics
using Polynomials
using CairoMakie

# ─── Strategy definitions ─────────────────────────────────────────────────

struct CutStrategySpec
    label::String
    cuts::Symbol
end

struct EvalStrategySpec
    label::String
    delta::Symbol
    std::Union{Nothing,Symbol}
end

const CUT_STRATEGY_SPECS = [
    CutStrategySpec("Heuristic",     :heuristic_cuts),
    CutStrategySpec("Random edge",   :random_cuts),
    CutStrategySpec("Random hoop",   :random_hoop_cuts),
    CutStrategySpec("Biggest hoop",  :biggest_hoop_cuts),
    CutStrategySpec("Smallest hoop", :smallest_hoop_cuts),
    CutStrategySpec("Max curvature", :max_curvature_cuts),
    CutStrategySpec("Oracle",        :oracle_cuts),
]

const EVAL_STRATEGY_SPECS = [
    EvalStrategySpec("Heuristic",     :delta_heuristic,      nothing),
    EvalStrategySpec("Random edge",   :delta_random,         :delta_random_std),
    EvalStrategySpec("Random hoop",   :delta_random_hoop,    :delta_random_hoop_std),
    EvalStrategySpec("Biggest hoop",  :delta_biggest_hoop,   nothing),
    EvalStrategySpec("Smallest hoop", :delta_smallest_hoop,  nothing),
    EvalStrategySpec("Max curvature", :delta_max_curvature,  nothing),
    EvalStrategySpec("Oracle",        :delta_oracle,         nothing),
]

# Canonical display order, shared by every cut-stage comparison.
const CUT_STRATEGY_LABELS = [spec.label for spec in CUT_STRATEGY_SPECS]

# Raw strategy names in the stage files → display labels.
const CUT_STRATEGY_LABEL_MAP = Dict(
    "heuristic" => "Heuristic",
    "random_edge" => "Random edge",
    "random_hoop" => "Random hoop",
    "biggest_hoop" => "Biggest hoop",
    "smallest_hoop" => "Smallest hoop",
    "max_curvature" => "Max curvature",
    "oracle" => "Oracle",
)

# ─── Shared data readers ──────────────────────────────────────────────────

"""
    _read_cuts_data(datafile)

Read a cuts-to-zero TSV and return the available strategy columns.  Rows with
`T_before == 0` or an unidentified `knotType` are removed.  Missing strategy
columns become `NaN` vectors.
"""
function _read_cuts_data(datafile)
    data, header = _read_table(datafile)

    t = _numeric_col(data, header, "T_before")
    keep = .!isnan.(t) .& (t .!= 0) .& _defined_knot_mask(data, header)
    col(name) = _numeric_col(data, header, name)[keep]

    return (
        T_before = t[keep],
        heuristic_cuts = col("heuristic_cuts"),
        random_cuts = col("random_cuts_mean"),
        random_hoop_cuts = col("random_hoop_cuts_mean"),
        biggest_hoop_cuts = col("biggest_hoop_cuts"),
        smallest_hoop_cuts = col("smallest_hoop_cuts"),
        max_curvature_cuts = col("max_curvature_cuts"),
        oracle_cuts = col("oracle_cuts"),
    )
end

"""
    _read_cut_eval_data(datafile)

Read a cut-evaluation TSV and return the available single-cut strategy columns.
Rows with `T_before == 0` or an unidentified `knotType` are removed.  Oracle
values are masked out when `oracle_evaluated` is not 1.
"""
function _read_cut_eval_data(datafile)
    data, header = _read_table(datafile)

    t = _numeric_col(data, header, "T_before")
    keep = .!isnan.(t) .& (t .!= 0) .& _defined_knot_mask(data, header)
    col(name) = _numeric_col(data, header, name)[keep]

    delta_oracle = col("delta_oracle")
    oracle_evaluated = col("oracle_evaluated")
    delta_oracle[oracle_evaluated .!= 1] .= NaN

    return (
        T_before = t[keep],
        delta_heuristic = col("delta_heuristic"),
        delta_random = col("delta_random_mean"),
        delta_random_std = col("delta_random_std"),
        delta_random_hoop = col("delta_random_hoop_mean"),
        delta_random_hoop_std = col("delta_random_hoop_std"),
        delta_biggest_hoop = col("delta_biggest_hoop"),
        delta_smallest_hoop = col("delta_smallest_hoop"),
        delta_max_curvature = col("delta_max_curvature"),
        delta_oracle = delta_oracle,
    )
end

function _count_zero_threading(datafile)
    data, header = _read_table(datafile)
    index = _col_index(header, "T_before")
    (index === nothing || isempty(data)) && return 0
    return count(value -> _to_f64(value) == 0.0, data[:, index])
end

# ─── Strategy accessors ───────────────────────────────────────────────────

function _available_cuts_strategies(d)
    strategies = NamedTuple[]

    for spec in CUT_STRATEGY_SPECS
        cuts = getproperty(d, spec.cuts)
        all(isnan.(cuts)) && continue
        push!(strategies, (label = spec.label, cuts = cuts))
    end

    return strategies
end

function _available_eval_strategies(d)
    n = length(d.T_before)
    strategies = NamedTuple[]

    for spec in EVAL_STRATEGY_SPECS
        delta = getproperty(d, spec.delta)
        all(isnan.(delta)) && continue
        std = spec.std === nothing ? fill(NaN, n) : getproperty(d, spec.std)
        push!(strategies, (label = spec.label, delta = delta, std = std))
    end

    return strategies
end

# ─── Default data file paths ──────────────────────────────────────────────

_cuts_datafile() = joinpath(@__DIR__, "..", "data", "knot_cuts_to_zero_v3.tsv")
_eval_datafile() = joinpath(@__DIR__, "..", "data", "knot_cut_eval_v3.tsv")

# ─── Histogram overlay ────────────────────────────────────────────────────

"""
    plot_cuts_histogram([datafile]; settings::Settings = DEFAULT_SETTINGS)

Overlaid histograms of cuts for the heuristic and every strategy whose paired
difference from it is significant.  Median lines are drawn without extra
legend entries; the terminal output lists means and medians.
"""
function plot_cuts_histogram(datafile = _cuts_datafile(); settings::Settings = DEFAULT_SETTINGS)
    d = _read_cuts_data(datafile)
    strategies = _significant_cuts_strategies(d)
    isempty(strategies) && error("No cut strategies found in $datafile")
    colors = _strategy_colors(length(strategies); palette_name = :tab10)

    # The other cut plots exclude T_before == 0 to focus on threaded tangles.
    # The histogram accounts for every catalog tangle, so add already-resolved
    # tangles back as zero-cut observations.
    zero_count = _count_zero_threading(datafile)

    all_values = Float64[]
    for strategy in strategies
        append!(all_values, strategy.cuts[isfinite.(strategy.cuts)])
        append!(all_values, zeros(zero_count))
    end
    isempty(all_values) && error("No finite cut values found in $datafile")

    min_value = floor(minimum(all_values))
    max_value = ceil(maximum(all_values))
    bins = (min_value - 0.5):1.0:(max_value + 0.5)

    fig = Figure(size = (800, 500), fontsize = 22)
    ax = Axis(fig[1, 1],
        xlabel = "Steps (cuts needed)", ylabel = "Count",
        xgridvisible = false, ygridvisible = false,
    )

    for (index, strategy) in enumerate(strategies)
        values = strategy.cuts[isfinite.(strategy.cuts)]
        isempty(values) && zero_count == 0 && continue
        values = vcat(values, zeros(zero_count))

        hist!(ax, values, bins = bins,
            color = (colors[index], 0.25),
            strokecolor = colors[index], strokewidth = 1.2,
            label = strategy.label,
        )
        vlines!(ax, [median(values)], color = colors[index],
            linestyle = :dash, linewidth = 2)

        println("$(strategy.label): mean = $(round(mean(values), digits=2)), " *
                "median = $(round(median(values), digits=2))")
    end

    ax.xticks = min_value:1:max_value
    ax.xtickformat = "{:.0f}"

    axislegend(ax, position = :rt)
    outpath = joinpath(settings.figures_dir, "cut_eval", "cuts_histogram.pdf")
    mkpath(dirname(outpath))
    save(outpath, fig)
    return fig
end

# ─── Single-cut comparison ────────────────────────────────────────────────

"""
    plot_cut_eval_comparison([datafile]; settings::Settings = DEFAULT_SETTINGS)

Scatter of each significant strategy's single-cut threading reduction against
the reference strategy (heuristic when available).  The identity line indicates
equal performance.  A linear fit is printed for each strategy.
"""
function plot_cut_eval_comparison(datafile = _eval_datafile(); settings::Settings = DEFAULT_SETTINGS)
    d = _read_cut_eval_data(datafile)
    strategies = _significant_eval_strategies(d)
    isempty(strategies) && error("No cut-evaluation strategies found in $datafile")
    reference = _reference_strategy(strategies)
    reference === nothing && error("No reference strategy available")
    colors = _strategy_colors(length(strategies))

    all_values = Float64[]
    for strategy in strategies
        append!(all_values, strategy.delta[isfinite.(strategy.delta)])
    end
    isempty(all_values) && error("No finite delta values found in $datafile")

    lo, hi = extrema(all_values)
    band = hi - lo
    band == 0 && (band = max(abs(hi), 1.0))
    lo -= 0.05 * band
    hi += 0.05 * band

    fig = Figure(size = (800, 600), fontsize = 22)
    ax = Axis(fig[1, 1],
        xlabel = "$(reference.label) delta",
        ylabel = "Strategy delta",
        aspect = DataAspect(),
        limits = (lo, hi, lo, hi),
    )

    lines!(ax, [lo, hi], [lo, hi], color = :gray, linestyle = :dash,
        linewidth = 1, label = "strategy = $(reference.label)")

    for (index, strategy) in enumerate(strategies)
        strategy.label == reference.label && continue

        mask = _finite_mask(reference.delta, strategy.delta)
        any(mask) || continue

        scatter!(ax, reference.delta[mask], strategy.delta[mask],
            color = colors[index], label = strategy.label, markersize = 6)

        if any(isfinite.(strategy.std))
            errmask = mask .& isfinite.(strategy.std)
            if any(errmask)
                errorbars!(ax, reference.delta[errmask], strategy.delta[errmask],
                    strategy.std[errmask], strategy.std[errmask],
                    color = (colors[index], 0.3), whiskerwidth = 3)
            end
        end

        x = float.(reference.delta[mask])
        y = float.(strategy.delta[mask])
        if length(x) >= 2 && length(unique(x)) >= 2
            coefficients = coeffs(Polynomials.fit(x, y, 1))
            println("[$(strategy.label)] delta = $(round(coefficients[1], digits=3)) + " *
                    "$(round(coefficients[2], digits=3)) · $(reference.label) delta")
        end
    end

    axislegend(ax, position = :lt)
    outpath = joinpath(settings.figures_dir, "cut_eval", "cut_eval_comparison.pdf")
    mkpath(dirname(outpath))
    save(outpath, fig)
    return fig
end

# ─── Strategy filtering for the overview plots ────────────────────────────

function _significant_cuts_strategies(d)
    strategies = _available_cuts_strategies(d)
    isempty(strategies) && return strategies
    rows = _paired_comparison_rows(
        strategies, "Heuristic", strategy -> strategy.cuts, false;
        metric = "cuts_to_zero", include_diagnostics = false)
    significant = Set(row.strategy for row in rows if row.significant)
    reference = _reference_strategy(strategies)
    return [strategy for strategy in strategies
            if strategy.label == reference.label || strategy.label in significant]
end

function _significant_eval_strategies(d)
    strategies = _available_eval_strategies(d)
    isempty(strategies) && return strategies
    rows = _paired_comparison_rows(
        strategies, "Heuristic", strategy -> strategy.delta, true;
        metric = "delta", include_diagnostics = false)
    significant = Set(row.strategy for row in rows if row.significant)
    reference = _reference_strategy(strategies)
    return [strategy for strategy in strategies
            if strategy.label == reference.label || strategy.label in significant]
end

# ─── Violin comparisons: threading reduction and cuts needed ───────────────
#
# Each metric is available in two variants: the first cut (single-cut stage
# data) and the complete cuts-to-zero sequence, restricted to resolved trials
# where that matters.

function _stage_columns(stage_file::AbstractString)
    isfile(stage_file) || return nothing
    data, header = _read_table(stage_file)
    isempty(data) && return nothing
    return (
        seed = _string_col(data, header, "seed"),
        strategy = _string_col(data, header, "strategy"),
        trial = _string_col(data, header, "trial"),
        knotType = _string_col(data, header, "knotType"),
        has_knot_type = _col_index(header, "knotType") !== nothing,
        stage = _numeric_col(data, header, "stage"),
        cuts = _numeric_col(data, header, "cuts"),
        threading = _numeric_col(data, header, "threading"),
        resolved = _numeric_col(data, header, "resolved"),
    )
end

_stage_label(raw) = get(CUT_STRATEGY_LABEL_MAP, raw, raw)

function _initial_threading(cols)
    initial = Dict{String, Float64}()
    for i in eachindex(cols.seed)
        cols.stage[i] == 0.0 || continue
        isfinite(cols.threading[i]) || continue
        initial[cols.seed[i]] = cols.threading[i]
    end
    return initial
end

# Seeds kept for cut/complexity analysis: they must start with threading and
# have a defined knot type.  Older files without a knotType column keep every
# threaded seed.
function _analysis_seed_set(cols)
    initial = _initial_threading(cols)
    isempty(initial) && return nothing

    defined = Dict{String, Bool}()
    if cols.has_knot_type
        for i in eachindex(cols.seed)
            seed = cols.seed[i]
            defined[seed] = get(defined, seed, false) ||
                is_defined_knot_type(cols.knotType[i])
        end
    end

    return Set(
        seed for (seed, value) in initial
        if isfinite(value) && value > 0 &&
           (!cols.has_knot_type || get(defined, seed, false))
    )
end

function _final_stage_rows(cols)
    final = Dict{Tuple{String, String, String}, NamedTuple}()
    for i in eachindex(cols.seed)
        key = (cols.seed[i], _stage_label(cols.strategy[i]), cols.trial[i])
        stage = cols.stage[i]
        if !haskey(final, key) || stage > final[key].stage
            final[key] = (
                stage = stage,
                cuts = cols.cuts[i],
                threading = cols.threading[i],
                resolved = cols.resolved[i],
            )
        end
    end
    return final
end

function _group_by_seed_strategy(grouped)
    seeds = sort(unique(first.(keys(grouped))))
    strategies = [label for label in CUT_STRATEGY_LABELS
                  if any(key -> key[2] == label, keys(grouped))]
    values = Dict(label => fill(NaN, length(seeds)) for label in strategies)
    for label in strategies
        for (index, seed) in enumerate(seeds)
            entries = get(grouped, (seed, label), nothing)
            entries === nothing && continue
            finite = entries[isfinite.(entries)]
            isempty(finite) || (values[label][index] = mean(finite))
        end
    end
    return (seeds = seeds, strategies = strategies, values = values)
end

"""
    first_cut_threading_reduction([stage_file])

Per-seed mean threading reduction after the first cut for every strategy.
Threading reduction is `T_before − T_after`; random strategies are averaged
over their trials.
"""
function first_cut_threading_reduction(
        stage_file::AbstractString = DEFAULT_SETTINGS.cut_eval_stage_file,
    )
    cols = _stage_columns(stage_file)
    cols === nothing && return nothing
    initial = _initial_threading(cols)
    analysis = _analysis_seed_set(cols)

    grouped = Dict{Tuple{String, String}, Vector{Float64}}()
    for i in eachindex(cols.seed)
        cols.stage[i] == 1.0 || continue
        isfinite(cols.threading[i]) || continue
        analysis === nothing || cols.seed[i] in analysis || continue
        before = get(initial, cols.seed[i], NaN)
        (isfinite(before) && before > 0) || continue
        label = _stage_label(cols.strategy[i])
        label in CUT_STRATEGY_LABELS || continue
        push!(get!(grouped, (cols.seed[i], label), Float64[]),
            before - cols.threading[i])
    end
    isempty(grouped) && return nothing
    return _group_by_seed_strategy(grouped)
end

"""
    resolved_threading_reduction([stage_file])

Per-seed mean threading reduction at the final cuts-to-zero stage, keeping
only trials that resolved threading to zero.
"""
function resolved_threading_reduction(
        stage_file::AbstractString = DEFAULT_SETTINGS.cuts_to_zero_stage_file,
    )
    cols = _stage_columns(stage_file)
    cols === nothing && return nothing
    initial = _initial_threading(cols)
    analysis = _analysis_seed_set(cols)
    final = _final_stage_rows(cols)

    grouped = Dict{Tuple{String, String}, Vector{Float64}}()
    for ((seed, label, _), row) in final
        row.resolved == 1.0 || continue
        label in CUT_STRATEGY_LABELS || continue
        analysis === nothing || seed in analysis || continue
        before = get(initial, seed, NaN)
        (isfinite(before) && before > 0) || continue
        isfinite(row.threading) || continue
        push!(get!(grouped, (seed, label), Float64[]), before - row.threading)
    end
    isempty(grouped) && return nothing
    return _group_by_seed_strategy(grouped)
end

"""
    cuts_needed_summary([stage_file]; resolved_only = false)

Per-seed mean number of cuts at the final cuts-to-zero stage.  With
`resolved_only = true` only trials that resolved threading to zero are kept.
"""
function cuts_needed_summary(
        stage_file::AbstractString = DEFAULT_SETTINGS.cuts_to_zero_stage_file;
        resolved_only::Bool = false,
    )
    cols = _stage_columns(stage_file)
    cols === nothing && return nothing
    final = _final_stage_rows(cols)
    analysis = _analysis_seed_set(cols)

    grouped = Dict{Tuple{String, String}, Vector{Float64}}()
    for ((seed, label, _), row) in final
        label in CUT_STRATEGY_LABELS || continue
        analysis === nothing || seed in analysis || continue
        resolved_only && row.resolved != 1.0 && continue
        isfinite(row.cuts) || continue
        push!(get!(grouped, (seed, label), Float64[]), row.cuts)
    end
    isempty(grouped) && return nothing
    return _group_by_seed_strategy(grouped)
end

# ── Split-violin view of the paired differences ───────────────────────────

# Per-seed paired differences (strategy − heuristic).  Seeds where either
# side is missing stay NaN so they are ignored consistently everywhere.
function _paired_differences(comparison)
    comparison === nothing && return nothing
    reference = get(comparison.values, "Heuristic", nothing)
    reference === nothing && return nothing

    diffs = Dict{String, Vector{Float64}}()
    for label in comparison.strategies
        label == "Heuristic" && continue
        values = comparison.values[label]
        diffs[label] = [
            isfinite(values[i]) && isfinite(reference[i]) ?
                values[i] - reference[i] : NaN
            for i in eachindex(values)
        ]
    end
    return diffs
end

_paired_conditional_path(outpath::String) = begin
    root, ext = splitext(outpath)
    root * "_paired" * ext
end

"""
    _plot_paired_conditional(comparison, outpath, title, metric; higher_is_better)

Split-violin view of the paired comparison.  For each strategy the left half
shows every paired difference (zeros included) and the right half shows only
the non-zero differences that enter the conditional signed-rank test.
"""
function _plot_paired_conditional(
        comparison,
        outpath::String,
        title::String,
        metric::String;
        higher_is_better::Bool = true,
    )
    comparison === nothing && return nothing
    isempty(comparison.strategies) && return nothing
    diffs = _paired_differences(comparison)
    diffs === nothing && return nothing

    strategies = [label for label in comparison.strategies if label != "Heuristic"]
    isempty(strategies) && return nothing

    # The differences are already relative to the heuristic, so the reference
    # is a zero vector and the finite mask selects the paired seeds.
    candidates = vcat(
        [(label = "Heuristic", delta = zeros(length(comparison.seeds)))],
        [(label = label, delta = diffs[label]) for label in strategies],
    )
    rows = _paired_comparison_rows(
        candidates, "Heuristic", strategy -> strategy.delta, higher_is_better;
        metric = metric, include_diagnostics = false,
    )
    row_map = Dict(row.strategy => row for row in rows)

    n = length(strategies)
    x = collect(1:n)
    palette = _strategy_colors(n)

    all_values = Float64[]
    for label in strategies
        append!(all_values, diffs[label][isfinite.(diffs[label])])
    end
    isempty(all_values) && return nothing

    lo = min(0.0, minimum(all_values))
    hi = max(0.0, maximum(all_values))
    if hi == lo
        # Every pair is tied (or there is no spread at all).  Give the panel a
        # readable range instead of collapsing onto a single value.
        lo, hi = -1.0, 1.0
    end
    span = hi - lo
    span == 0.0 && (span = 1.0)

    direction_note = higher_is_better ?
        "positive = strategy removes more threading (better)" :
        "negative = strategy needs fewer cuts (better)"

    fig = Figure(size = (1500, 700), fontsize = 17)
    ax = Axis(
        fig[1, 1],
        xlabel = "Strategy (paired against Heuristic)",
        ylabel = "Paired difference (strategy − heuristic)\n" * direction_note,
        title = title * "\n" *
            "left: all paired differences (zeros included) · " *
            "right: non-zero differences only",
        xticks = (x, strategies),
        xticklabelrotation = π / 4,
    )

    # Annotations must clear the density tails drawn for both violin halves.
    band_top = hi
    for label in strategies
        finite_diff = diffs[label][isfinite.(diffs[label])]
        isempty(finite_diff) && continue
        nonzero_diff = finite_diff[abs.(finite_diff) .> 1e-12]
        if isempty(nonzero_diff)
            band_top = max(band_top, maximum(finite_diff))
        else
            band_top = max(band_top, _violin_top(finite_diff), _violin_top(nonzero_diff))
        end
    end

    for (index, label) in enumerate(strategies)
        all_diff = diffs[label][isfinite.(diffs[label])]
        isempty(all_diff) && continue
        nonzero_diff = all_diff[abs.(all_diff) .> 1e-12]
        color = palette[index]

        if isempty(nonzero_diff)
            # All pairs are exact ties, so a density estimate would just be a
            # degenerate block.  Mark the tie mass with a flat bar instead.
            lines!(ax, [index - 0.34, index + 0.34], [0.0, 0.0];
                color = color, linewidth = 7)
        else
            violin!(
                ax, fill(Float64(index), length(all_diff)), all_diff;
                side = :left, width = 0.8, scale = :count,
                color = (color, 0.35), strokecolor = color, strokewidth = 1.1,
                show_median = true, mediancolor = :black, medianlinewidth = 1.8,
            )
            violin!(
                ax, fill(Float64(index), length(nonzero_diff)), nonzero_diff;
                side = :right, width = 0.8, scale = :count,
                color = (color, 0.75), strokecolor = color, strokewidth = 1.1,
                show_median = true, mediancolor = :black, medianlinewidth = 1.8,
            )
        end

        row = get(row_map, label, nothing)
        row === nothing && continue

        tied_pct = round(100 * row.tied_rate, digits = 1)
        discordant_pct = round(100 * row.n_discordant / row.n, digits = 1)
        label_text = "ties $(row.n_tied) ($(tied_pct)%) · " *
                     "disc $(row.n_discordant) ($(discordant_pct)%)\n" *
                     "$(row.n_improved) better · $(row.n_worse) worse"
        text!(ax, index, band_top + 0.03 * span;
            text = label_text, align = (:center, :bottom),
            fontsize = 12, color = :gray30)

        row.significant || continue
        marker_text = row.p_adjusted < 0.001 ? "***" :
            row.p_adjusted < 0.01 ? "**" : "*"
        marker_color = row.better == row.strategy ? :red :
            row.better == row.reference ? :green : :gray
        effect = _fmt_effect(row.effect_size)
        annotation = isempty(effect) ?
            marker_text : marker_text * "  r = " * effect
        text!(ax, index, band_top + 0.14 * span;
            text = annotation, align = (:center, :bottom),
            fontsize = 15, font = :bold, color = marker_color)
    end

    hlines!(ax, [0], color = :gray, linestyle = :dash, linewidth = 1)
    ylims!(ax, lo - 0.08 * span, band_top + 0.24 * span)

    if all(label -> all(abs.(diffs[label][isfinite.(diffs[label])]) .<= 1e-12), strategies)
        text!(ax, (1 + n) / 2, lo + 0.18 * span;
            text = "all paired differences are exactly zero — no discordant pairs, so no conditional test",
            align = (:center, :bottom), fontsize = 15,
            color = :gray30, font = :italic)
    end

    legend_elements = [
        PolyElement(color = (:gray, 0.35), strokecolor = :gray, strokewidth = 1),
        PolyElement(color = (:gray, 0.75), strokecolor = :gray, strokewidth = 1),
        MarkerElement(marker = '*', color = :red, markersize = 22),
        MarkerElement(marker = '*', color = :green, markersize = 22),
    ]
    Legend(
        fig[2, 1],
        legend_elements,
        [
            "all paired differences (zeros included)",
            "non-zero paired differences (conditional sample)",
            "strategy better than Heuristic",
            "Heuristic better than strategy",
        ],
        "Conditional Wilcoxon vs Heuristic (BH-adjusted):\n" *
        "* p<0.05, ** p<0.01, *** p<0.001\n" *
        "r = matched-pairs rank-biserial (strategy − Heuristic)";
        orientation = :horizontal, framevisible = false, tellwidth = false,
        labelsize = 16, titlesize = 16, titleposition = :top,
    )

    mkpath(dirname(outpath))
    save(outpath, fig)
    return fig
end

# ── Public violin entry points ────────────────────────────────────────────

function _plot_cuts_violin(
        comparison,
        outpath::String,
        title::String,
        ylabel::String,
        metric::String;
        higher_is_better::Bool = true,
        zero_line::Bool = true,
    )
    comparison === nothing && error("no stage metrics available for $metric")
    isempty(comparison.strategies) && error("no strategies available for $metric")
    return _plot_metric_violin(
        comparison.strategies, comparison.values;
        outpath = outpath, title = title, ylabel = ylabel, metric = metric,
        comparison = comparison,
        higher_is_better = higher_is_better,
        zero_line = zero_line,
    )
end

"""
    plot_threading_reduction_first_cut([stage_file]; settings)

Violin comparison of the threading reduction after the first cut for every
strategy.  Strategies whose paired Wilcoxon test against the heuristic is
significant (Benjamini-Hochberg corrected) are marked with red asterisks when
the strategy is better and green asterisks when the heuristic is better, with
the matched-pairs rank-biserial effect size printed above the marker.
"""
function plot_threading_reduction_first_cut(
        stage_file::AbstractString = DEFAULT_SETTINGS.cut_eval_stage_file;
        settings::Settings = DEFAULT_SETTINGS,
    )
    comparison = first_cut_threading_reduction(stage_file)
    outpath = joinpath(settings.figures_dir, "cut_eval", "threading_reduction_first_cut.pdf")
    fig = _plot_cuts_violin(
        comparison, outpath,
        "Threading reduction after the first cut",
        "Threading reduction (ΔT)", "threading_reduction_first_cut";
        higher_is_better = true, zero_line = true,
    )
    _plot_paired_conditional(
        comparison, _paired_conditional_path(outpath),
        "Threading reduction after the first cut", "threading_reduction_first_cut";
        higher_is_better = true,
    )
    return fig
end

"""
    plot_threading_reduction_resolved([stage_file]; settings)

Violin comparison of the total threading reduction after the cuts-to-zero
sequence has resolved threading to zero.  Only resolved trials are used.
"""
function plot_threading_reduction_resolved(
        stage_file::AbstractString = DEFAULT_SETTINGS.cuts_to_zero_stage_file;
        settings::Settings = DEFAULT_SETTINGS,
    )
    comparison = resolved_threading_reduction(stage_file)
    outpath = joinpath(settings.figures_dir, "cut_eval", "threading_reduction_resolved.pdf")
    fig = _plot_cuts_violin(
        comparison, outpath,
        "Threading reduction after complete resolution",
        "Threading reduction (ΔT)", "threading_reduction_resolved";
        higher_is_better = true, zero_line = true,
    )
    _plot_paired_conditional(
        comparison, _paired_conditional_path(outpath),
        "Threading reduction after complete resolution", "threading_reduction_resolved";
        higher_is_better = true,
    )
    return fig
end

"""
    plot_cuts_needed_all([stage_file]; settings)

Violin comparison of the number of cuts needed at the final cuts-to-zero
stage, including trials that did not resolve threading.  Lower is better.
"""
function plot_cuts_needed_all(
        stage_file::AbstractString = DEFAULT_SETTINGS.cuts_to_zero_stage_file;
        settings::Settings = DEFAULT_SETTINGS,
    )
    comparison = cuts_needed_summary(stage_file; resolved_only = false)
    outpath = joinpath(settings.figures_dir, "cut_eval", "cuts_needed_all.pdf")
    fig = _plot_cuts_violin(
        comparison, outpath,
        "Cuts needed for all trials", "Cuts needed", "cuts_needed_all";
        higher_is_better = false, zero_line = false,
    )
    _plot_paired_conditional(
        comparison, _paired_conditional_path(outpath),
        "Cuts needed for all trials", "cuts_needed_all";
        higher_is_better = false,
    )
    return fig
end

"""
    plot_cuts_needed_resolved([stage_file]; settings)

Violin comparison of the number of cuts needed for trials that resolved
threading to zero.  Lower is better.
"""
function plot_cuts_needed_resolved(
        stage_file::AbstractString = DEFAULT_SETTINGS.cuts_to_zero_stage_file;
        settings::Settings = DEFAULT_SETTINGS,
    )
    comparison = cuts_needed_summary(stage_file; resolved_only = true)
    outpath = joinpath(settings.figures_dir, "cut_eval", "cuts_needed_resolved.pdf")
    fig = _plot_cuts_violin(
        comparison, outpath,
        "Cuts needed for trials that resolved threading", "Cuts needed",
        "cuts_needed_resolved";
        higher_is_better = false, zero_line = false,
    )
    _plot_paired_conditional(
        comparison, _paired_conditional_path(outpath),
        "Cuts needed for trials that resolved threading", "cuts_needed_resolved";
        higher_is_better = false,
    )
    return fig
end

# ── Paired-test tables for the four violins ───────────────────────────────

function _write_threading_cuts_paired_tests(
        comparison,
        settings::Settings,
        filename::String,
        metric::String,
        title::String;
        higher_is_better::Bool,
    )
    comparison === nothing && return nothing
    rows = _paired_rows_for_metric(
        comparison, comparison.values, metric;
        higher_is_better = higher_is_better,
    )
    if isempty(rows)
        @warn "not enough paired observations for $metric"
        return rows
    end
    outpath = joinpath(settings.analysis_dir, filename)
    _write_paired_report(rows, outpath)
    _print_paired_report(rows, title)
    println("\nSaved paired test table to $outpath")
    return rows
end

"""
    threading_reduction_first_cut_paired_tests([stage_file]; settings)

Paired Wilcoxon comparisons of first-cut threading reduction against the
heuristic strategy.
"""
function threading_reduction_first_cut_paired_tests(
        stage_file::AbstractString = DEFAULT_SETTINGS.cut_eval_stage_file;
        settings::Settings = DEFAULT_SETTINGS,
    )
    return _write_threading_cuts_paired_tests(
        first_cut_threading_reduction(stage_file), settings,
        "threading_reduction_first_cut_paired_tests.tsv",
        "threading_reduction_first_cut",
        "Threading reduction after first cut vs heuristic";
        higher_is_better = true,
    )
end

"""
    threading_reduction_resolved_paired_tests([stage_file]; settings)

Paired Wilcoxon comparisons of threading reduction after complete resolution
against the heuristic strategy.
"""
function threading_reduction_resolved_paired_tests(
        stage_file::AbstractString = DEFAULT_SETTINGS.cuts_to_zero_stage_file;
        settings::Settings = DEFAULT_SETTINGS,
    )
    return _write_threading_cuts_paired_tests(
        resolved_threading_reduction(stage_file), settings,
        "threading_reduction_resolved_paired_tests.tsv",
        "threading_reduction_resolved",
        "Threading reduction after resolution vs heuristic";
        higher_is_better = true,
    )
end

"""
    cuts_needed_paired_tests([stage_file]; settings)

Paired Wilcoxon comparisons of the cuts needed for all trials against the
heuristic strategy.  Negative differences mean fewer cuts than the heuristic.
"""
function cuts_needed_paired_tests(
        stage_file::AbstractString = DEFAULT_SETTINGS.cuts_to_zero_stage_file;
        settings::Settings = DEFAULT_SETTINGS,
    )
    return _write_threading_cuts_paired_tests(
        cuts_needed_summary(stage_file; resolved_only = false), settings,
        "cuts_needed_paired_tests.tsv",
        "cuts_needed_all",
        "Cuts needed (all trials) vs heuristic";
        higher_is_better = false,
    )
end

"""
    cuts_needed_resolved_paired_tests([stage_file]; settings)

Paired Wilcoxon comparisons of the cuts needed for resolved trials against the
heuristic strategy.  Negative differences mean fewer cuts than the heuristic.
"""
function cuts_needed_resolved_paired_tests(
        stage_file::AbstractString = DEFAULT_SETTINGS.cuts_to_zero_stage_file;
        settings::Settings = DEFAULT_SETTINGS,
    )
    return _write_threading_cuts_paired_tests(
        cuts_needed_summary(stage_file; resolved_only = true), settings,
        "cuts_needed_resolved_paired_tests.tsv",
        "cuts_needed_resolved",
        "Cuts needed (resolved trials) vs heuristic";
        higher_is_better = false,
    )
end

"""
    cuts_violin_comparisons(settings::Settings = DEFAULT_SETTINGS)

Generate the four violin comparisons for threading reduction and cuts needed,
and write their paired-test tables.  Missing stage files are skipped.
"""
function cuts_violin_comparisons(settings::Settings = DEFAULT_SETTINGS)
    if isfile(settings.cut_eval_stage_file)
        plot_threading_reduction_first_cut(settings.cut_eval_stage_file; settings = settings)
        threading_reduction_first_cut_paired_tests(
            settings.cut_eval_stage_file; settings = settings)
    else
        @warn "no first-cut stage metrics found; skipping threading violin plot" settings.cut_eval_stage_file
    end

    if isfile(settings.cuts_to_zero_stage_file)
        for (plot_function, test_function) in (
            (plot_threading_reduction_resolved, threading_reduction_resolved_paired_tests),
            (plot_cuts_needed_all, cuts_needed_paired_tests),
            (plot_cuts_needed_resolved, cuts_needed_resolved_paired_tests),
        )
            plot_function(settings.cuts_to_zero_stage_file; settings = settings)
            test_function(settings.cuts_to_zero_stage_file; settings = settings)
        end
    else
        @warn "no cuts-to-zero stage metrics found; skipping resolved violins" settings.cuts_to_zero_stage_file
    end
    return nothing
end
