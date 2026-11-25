# ── Paired strategy statistics and violin plots ────────────────────────────
#
# Shared by the cut violins (cuts.jl) and the complexity violins
# (complexity.jl).  Every strategy is compared with the heuristic on the same
# seeds, so all tests operate on paired differences `strategy − heuristic`.

using Statistics
using Distributions
using HypothesisTests
using MultipleTesting
using Random
using CairoMakie

# ── Signed-rank helpers ────────────────────────────────────────────────────

function _average_ranks(values)
    n = length(values)
    order = sortperm(values)
    ranks = zeros(Float64, n)
    i = 1
    while i <= n
        j = i
        while j < n && values[order[j + 1]] == values[order[j]]
            j += 1
        end
        average = (i + j) / 2
        for k in i:j
            ranks[order[k]] = average
        end
        i = j + 1
    end
    return ranks
end

function _rank_biserial(diff)
    nonzero = abs.(diff) .> 1e-12
    any(nonzero) || return 0.0
    d = diff[nonzero]
    ranks = _average_ranks(abs.(d))
    positive = sum(ranks[d .> 0])
    negative = sum(ranks[d .< 0])
    denominator = positive + negative
    return denominator == 0 ? 0.0 : (positive - negative) / denominator
end

function _paired_pvalue(diff)
    nonzero = abs.(diff) .> 1e-12
    any(nonzero) || return 1.0
    try
        return pvalue(SignedRankTest(diff[nonzero]))
    catch
        return NaN
    end
end

function _sign_pvalue(diff)
    finite = diff[isfinite.(diff)]
    nonzero = finite[abs.(finite) .> 1e-12]
    length(nonzero) < 2 && return 1.0
    try
        return pvalue(SignTest(nonzero))
    catch
        return NaN
    end
end

# One-sided exact binomial test over all pairs, with ties counted as failures.
# This is the unconditional sensitivity check; it does not drive plot markers.
function _unconditional_better_pvalue(n_improved::Int, n::Int)
    (n <= 0 || n_improved <= 0) && return 1.0
    return ccdf(Binomial(n, 0.5), n_improved - 1)
end

function _ks_statistic(x, y)
    xs = sort(collect(x))
    ys = sort(collect(y))
    n = length(xs)
    m = length(ys)
    (n == 0 || m == 0) && return 0.0

    i = 1
    j = 1
    distance = 0.0
    while i <= n && j <= m
        if xs[i] <= ys[j]
            i += 1
        else
            j += 1
        end
        distance = max(distance, abs((i - 1) / n - (j - 1) / m))
    end
    return distance
end

"""
    _symmetry_pvalue(diff; nboot, seed)

Randomization diagnostic for symmetry of the paired differences about their
median.  The test statistic is the two-sample KS distance between the centred
differences and their reflection; the null distribution is obtained by random
sign flips.
"""
function _symmetry_pvalue(diff; nboot = 1000, seed = 20240505)
    finite = diff[isfinite.(diff)]
    length(finite) < 3 && return NaN

    centred = finite .- median(finite)
    observed = _ks_statistic(centred, -centred)
    observed <= 1e-12 && return 1.0

    rng = MersenneTwister(seed)
    more_extreme = 0
    for _ in 1:nboot
        flipped = ifelse.(rand(rng, Bool, length(centred)), centred, -centred)
        _ks_statistic(flipped, -flipped) >= observed && (more_extreme += 1)
    end
    return (more_extreme + 1) / (nboot + 1)
end

function _bootstrap_median_ci(diff; nboot = 10_000, seed = 20240505)
    n = length(diff)
    n == 0 && return (NaN, NaN)
    rng = MersenneTwister(seed)
    boot = Vector{Float64}(undef, nboot)
    for b in 1:nboot
        boot[b] = median(diff[rand(rng, 1:n, n)])
    end
    q = quantile(boot, [0.025, 0.975])
    return (q[1], q[2])
end

function _paired_effect_direction(diff, effect_size, mean_diff, median_nonzero_diff)
    # The median of all differences is a poor direction cue when many pairs are
    # exact ties, so use the rank-biserial effect size first, then the sign-test
    # majority, then the mean and the median of non-zero differences.
    tolerance = 1e-12
    if isfinite(effect_size) && abs(effect_size) > tolerance
        return effect_size > 0 ? :positive : :negative
    end

    n_positive = count(x -> x > tolerance, diff)
    n_negative = count(x -> x < -tolerance, diff)
    if n_positive != n_negative
        return n_positive > n_negative ? :positive : :negative
    end

    if isfinite(mean_diff) && abs(mean_diff) > tolerance
        return mean_diff > 0 ? :positive : :negative
    end

    if isfinite(median_nonzero_diff) && abs(median_nonzero_diff) > tolerance
        return median_nonzero_diff > 0 ? :positive : :negative
    end

    return :tie
end

# ── Paired comparison table ────────────────────────────────────────────────

function _reference_strategy(strategies; label = "Heuristic")
    isempty(strategies) && return nothing
    index = findfirst(strategy -> strategy.label == label, strategies)
    return index === nothing ? strategies[1] : strategies[index]
end

function _paired_comparison_rows(
        strategies,
        reference_label::String,
        value_getter,
        higher_is_better::Bool;
        metric::String,
        seed_base::Int = 20240505,
        include_diagnostics::Bool = true,
    )
    reference = _reference_strategy(strategies; label = reference_label)
    reference === nothing && return NamedTuple[]

    reference_values = value_getter(reference)
    records = NamedTuple[]
    pvalues = Float64[]
    unconditional_pvalues = Float64[]
    symmetry_pvalues = Float64[]
    sign_pvalues = Float64[]

    for (index, strategy) in enumerate(strategies)
        strategy.label == reference.label && continue
        values = value_getter(strategy)
        mask = isfinite.(values) .& isfinite.(reference_values)
        n = count(mask)
        n < 2 && continue

        diff = values[mask] .- reference_values[mask]
        median_diff = median(diff)
        mean_diff = mean(diff)
        nonzero = abs.(diff) .> 1e-12
        median_nonzero_diff = any(nonzero) ? median(diff[nonzero]) : NaN
        ci_low, ci_high = _bootstrap_median_ci(diff; seed = seed_base + index)
        p = _paired_pvalue(diff)
        effect_size = _rank_biserial(diff)
        symmetry_p = include_diagnostics ?
            _symmetry_pvalue(diff; seed = seed_base + index) : NaN
        sign_p = include_diagnostics ? _sign_pvalue(diff) : NaN

        tolerance = 1e-12
        if higher_is_better
            n_improved = count(x -> x > tolerance, diff)
            n_worse = count(x -> x < -tolerance, diff)
        else
            n_improved = count(x -> x < -tolerance, diff)
            n_worse = count(x -> x > tolerance, diff)
        end
        n_tied = n - n_improved - n_worse
        n_discordant = n_improved + n_worse
        probability_of_superiority = (n_improved + 0.5 * n_tied) / n
        net_improvement_rate = (n_improved - n_worse) / n
        unconditional_p = _unconditional_better_pvalue(n_improved, n)

        direction = _paired_effect_direction(
            diff, effect_size, mean_diff, median_nonzero_diff,
        )
        better = if direction === :tie
            "tie"
        elseif direction === :positive
            higher_is_better ? strategy.label : reference.label
        else
            higher_is_better ? reference.label : strategy.label
        end

        push!(records, (
            metric = metric,
            strategy = strategy.label,
            reference = reference.label,
            n = n,
            n_discordant = n_discordant,
            median_diff = median_diff,
            mean_diff = mean_diff,
            median_nonzero_diff = median_nonzero_diff,
            n_improved = n_improved,
            n_tied = n_tied,
            n_worse = n_worse,
            improved_rate = n_improved / n,
            tied_rate = n_tied / n,
            worse_rate = n_worse / n,
            probability_of_superiority = probability_of_superiority,
            net_improvement_rate = net_improvement_rate,
            ci_low = ci_low,
            ci_high = ci_high,
            p_value = p,
            unconditional_p = unconditional_p,
            effect_size = effect_size,
            better = better,
            symmetry_p = symmetry_p,
            sign_p = sign_p,
        ))
        push!(pvalues, p)
        push!(unconditional_pvalues, unconditional_p)
        push!(symmetry_pvalues, symmetry_p)
        push!(sign_pvalues, sign_p)
    end

    isempty(records) && return NamedTuple[]

    adjusted(values) = begin
        result = fill(NaN, length(values))
        valid = findall(isfinite.(values))
        isempty(valid) || (result[valid] = adjust(values[valid], BenjaminiHochberg()))
        result
    end

    p_adjusted = adjusted(pvalues)
    unconditional_adjusted = adjusted(unconditional_pvalues)
    symmetry_adjusted = adjusted(symmetry_pvalues)
    sign_adjusted = adjusted(sign_pvalues)

    return [
        (
            metric = row.metric,
            strategy = row.strategy,
            reference = row.reference,
            n = row.n,
            n_discordant = row.n_discordant,
            median_diff = row.median_diff,
            mean_diff = row.mean_diff,
            median_nonzero_diff = row.median_nonzero_diff,
            n_improved = row.n_improved,
            n_tied = row.n_tied,
            n_worse = row.n_worse,
            improved_rate = row.improved_rate,
            tied_rate = row.tied_rate,
            worse_rate = row.worse_rate,
            probability_of_superiority = row.probability_of_superiority,
            net_improvement_rate = row.net_improvement_rate,
            ci_low = row.ci_low,
            ci_high = row.ci_high,
            p_value = row.p_value,
            p_adjusted = p_adjusted[index],
            unconditional_p_value = row.unconditional_p,
            unconditional_p_adjusted = unconditional_adjusted[index],
            effect_size = row.effect_size,
            better = row.better,
            significant = isfinite(p_adjusted[index]) && p_adjusted[index] < 0.05,
            symmetry_p = row.symmetry_p,
            symmetry_p_adjusted = symmetry_adjusted[index],
            symmetry_ok = !isfinite(symmetry_adjusted[index]) ||
                symmetry_adjusted[index] >= 0.05,
            sign_p_value = row.sign_p,
            sign_p_adjusted = sign_adjusted[index],
        )
        for (index, row) in enumerate(records)
    ]
end

function _paired_rows_for_metric(
        comparison,
        values::Dict{String, Vector{Float64}},
        metric::String;
        higher_is_better::Bool = true,
    )
    strategies = [(label = label, delta = values[label]) for label in comparison.strategies]
    return _paired_comparison_rows(
        strategies, "Heuristic", strategy -> strategy.delta, higher_is_better;
        metric = metric,
    )
end

# ── Reports ────────────────────────────────────────────────────────────────

function _write_paired_report(rows, outpath::String)
    header = [
        "metric", "strategy", "reference", "n", "n_discordant",
        "median_diff", "mean_diff", "median_nonzero_diff",
        "n_improved", "n_tied", "n_worse",
        "improved_rate", "tied_rate", "worse_rate",
        "probability_of_superiority", "net_improvement_rate",
        "ci_low", "ci_high", "p_value", "p_adjusted",
        "unconditional_p_value", "unconditional_p_adjusted",
        "symmetry_p", "symmetry_p_adjusted", "symmetry_ok",
        "sign_p_value", "sign_p_adjusted",
        "effect_size", "better", "significant",
    ]
    rounded_6(value) = string(round(value, digits = 6))
    rounded_6_sig(value) = string(round(value, sigdigits = 6))

    rows_out = [
        [
            row.metric,
            row.strategy,
            row.reference,
            string(row.n),
            string(row.n_discordant),
            rounded_6(row.median_diff),
            rounded_6(row.mean_diff),
            rounded_6(row.median_nonzero_diff),
            string(row.n_improved),
            string(row.n_tied),
            string(row.n_worse),
            rounded_6(row.improved_rate),
            rounded_6(row.tied_rate),
            rounded_6(row.worse_rate),
            rounded_6(row.probability_of_superiority),
            rounded_6(row.net_improvement_rate),
            rounded_6(row.ci_low),
            rounded_6(row.ci_high),
            rounded_6_sig(row.p_value),
            rounded_6_sig(row.p_adjusted),
            rounded_6_sig(row.unconditional_p_value),
            rounded_6_sig(row.unconditional_p_adjusted),
            rounded_6_sig(row.symmetry_p),
            rounded_6_sig(row.symmetry_p_adjusted),
            string(row.symmetry_ok),
            rounded_6_sig(row.sign_p_value),
            rounded_6_sig(row.sign_p_adjusted),
            rounded_6(row.effect_size),
            row.better,
            string(row.significant),
        ]
        for row in rows
    ]
    _write_table(outpath, header, rows_out)
end

function _print_paired_report(rows, title::String)
    println("\n=== $title ===")
    println("median_diff = strategy − heuristic; effect_size = matched-pairs rank-biserial")
    println("CI = bootstrap 95% interval for the median difference")
    println("direction = effect size / sign-test majority (not the median of all differences)")
    println("symmetry = randomization KS diagnostic; sign = sign-test p-value")
    println("uncond_p = one-sided binomial over all pairs, ties as failures (sensitivity only)\n")
    for row in rows
        suffix = row.significant ? " *" : ""
        println("$(row.strategy) vs $(row.reference): " *
                "median Δ = $(round(row.median_diff, digits=2)), " *
                "median non-zero Δ = $(round(row.median_nonzero_diff, digits=2)), " *
                "improved/tied/worse = $(row.n_improved)/$(row.n_tied)/$(row.n_worse) " *
                "($(round(100 * row.improved_rate, digits=1))%/" *
                "$(round(100 * row.tied_rate, digits=1))%/" *
                "$(round(100 * row.worse_rate, digits=1))%), " *
                "P(superiority) = $(round(row.probability_of_superiority, digits=3)), " *
                "95% CI [$(round(row.ci_low, digits=2)), $(round(row.ci_high, digits=2))], " *
                "p_adj = $(round(row.p_adjusted, sigdigits=3)), " *
                "uncond_p_adj = $(round(row.unconditional_p_adjusted, sigdigits=3)), " *
                "symmetry_ok = $(row.symmetry_ok), " *
                "sign_p_adj = $(round(row.sign_p_adjusted, sigdigits=3)), " *
                "effect = $(round(row.effect_size, digits=3)), " *
                "better = $(row.better)$suffix")
    end
end

# ── Significance markers ───────────────────────────────────────────────────

function _significance_markers(
        comparison,
        values::Dict{String, Vector{Float64}},
        metric::String;
        higher_is_better::Bool = true,
    )
    rows = _paired_rows_for_metric(comparison, values, metric; higher_is_better = higher_is_better)
    markers = Dict{
        String,
        NamedTuple{(:text, :color, :effect), Tuple{String, Symbol, Float64}},
    }()
    for row in rows
        row.significant || continue
        text = row.p_adjusted < 0.001 ? "***" :
            row.p_adjusted < 0.01 ? "**" : "*"
        color = row.better == row.strategy ? :red :
            row.better == row.reference ? :green : :gray
        markers[row.strategy] = (text = text, color = color, effect = row.effect_size)
    end
    return markers
end

# ── Violin plot ────────────────────────────────────────────────────────────

function _plot_metric_violin(
        strategies::Vector{String},
        values::Dict{String, Vector{Float64}};
        outpath::String,
        title::String,
        ylabel::String,
        metric::String,
        comparison,
        higher_is_better::Bool = true,
        zero_line::Bool = true,
    )
    n = length(strategies)
    x = collect(1:n)
    palette = _strategy_colors(n)

    all_values = Float64[]
    for label in strategies
        append!(all_values, values[label][isfinite.(values[label])])
    end
    isempty(all_values) && error("no finite values for $metric")

    lo = minimum(all_values)
    hi = maximum(all_values)
    if zero_line
        lo = min(0.0, lo)
        hi = max(0.0, hi)
    end
    span = hi - lo
    span == 0.0 && (span = max(abs(hi), 1.0))

    markers = _significance_markers(
        comparison, values, metric; higher_is_better = higher_is_better,
    )

    # Annotations must clear the density tails, not just the largest value.
    density_top = hi
    for label in strategies
        top = _violin_top(values[label])
        isfinite(top) && (density_top = max(density_top, top))
    end

    fig = Figure(size = (1000, 440), fontsize = 18)
    ax = Axis(
        fig[1, 1],
        xlabel = "Cut strategy",
        ylabel = ylabel,
        title = title,
        xticks = (x, strategies),
        xticklabelrotation = π / 4,
    )

    println("  $(ylabel):")
    for (index, label) in enumerate(strategies)
        finite_values = values[label][isfinite.(values[label])]
        isempty(finite_values) && continue
        color = palette[index]
        violin!(
            ax,
            fill(Float64(index), length(finite_values)),
            finite_values;
            width = 0.78,
            color = (color, 0.55),
            strokecolor = color,
            strokewidth = 1.2,
            show_median = true,
            mediancolor = :black,
            medianlinewidth = 2.0,
            scale = :count,
        )

        suffix = haskey(markers, label) ?
            " $(markers[label].text) ($(markers[label].color))" : ""
        println("    $(label): n=$(length(finite_values)), " *
                "median=$(round(median(finite_values), digits=3))$suffix")

        haskey(markers, label) || continue
        marker = markers[label]
        # Anchor every marker on the same height above the tallest violin, so
        # the text never sits inside a density tail.
        effect = _fmt_effect(marker.effect)
        annotation = isempty(effect) ?
            marker.text : marker.text * "  r = " * effect
        text!(
            ax,
            index,
            density_top + 0.02 * span;
            text = annotation,
            align = (:center, :bottom),
            fontsize = 15,
            font = :bold,
            color = marker.color,
        )
    end

    zero_line && hlines!(ax, [0], color = :gray, linestyle = :dash, linewidth = 1)
    top_pad = isempty(markers) ? 0.06 : 0.16
    ylims!(ax, lo - 0.06 * span, max(hi, density_top) + top_pad * span)

    mkpath(dirname(outpath))
    save(outpath, fig)
    return fig
end
