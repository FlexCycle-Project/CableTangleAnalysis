# ── PCA-derived complexity score ───────────────────────────────────────────
#
# The complexity score is a PCA projection of five geometric metrics that are
# saved for every tangle.  CableSim writes these metrics for the initial tangle
# and for every cutting stage, so the same PCA model can score all stages.

using LinearAlgebra
using Statistics

const COMPLEXITY_FEATURES = [
    "segmentCount",
    "selfCollisionCount",
    "closePairCount",
    "totalCurvature",
    "outsideKnotCurvature",
]

const COMPLEXITY_UNIDENTIFIED = Set(["", "none", "None", "NONE", "??", "unknown", "Unknown"])
const COMPLEXITY_UNKNOTS = Set(["UN", "unknot", "Unknot"])

"""
    ComplexityModel

Fitted PC1 model used to score cable tangles.  `weights` is the weight vector
used for the score.  With `orientation_mode = :positive` the weights are the
absolute PC1 loadings so every metric increases complexity.  With `:signed`
the weights are the raw signed PC1 loadings.
"""
struct ComplexityModel
    features::Vector{String}
    mean::Vector{Float64}
    std::Vector{Float64}
    weights::Vector{Float64}
    signed_loadings::Vector{Float64}
    explained::Float64
    orientation::String
    orientation_mode::String
    raw_min::Float64
    raw_max::Float64
end

struct ComplexityScores
    pc1::Vector{Float64}
    pc1_signed::Vector{Float64}
    score::Vector{Float64}
    percentile::Vector{Float64}
end

# ── Catalog loading and filtering ──────────────────────────────────────────

function _catalog_records(catalog::AbstractString, features::Vector{String})
    data, header = _read_table(catalog)
    seed_index = _require_column(header, "seed")
    knot_index = _require_column(header, "knotType")
    feature_indices = [_require_column(header, feature) for feature in features]

    records = Dict{String, Any}()
    for row in 1:size(data, 1)
        seed = string(data[row, seed_index])
        isempty(seed) && continue

        values = _to_f64.(data[row, feature_indices])
        all(isfinite, values) || continue

        knot_type = string(data[row, knot_index])
        if haskey(records, seed)
            push!(records[seed]["knot_types"], knot_type)
        else
            records[seed] = Dict(
                "seed" => seed,
                "knot_types" => Set([knot_type]),
                "values" => values,
            )
        end
    end

    return collect(values(records))
end

_is_identified(types) = any(knot -> !(knot in COMPLEXITY_UNIDENTIFIED), types)

function _is_non_unknot(types)
    return any(
        knot -> !(knot in COMPLEXITY_UNIDENTIFIED) && !(knot in COMPLEXITY_UNKNOTS),
        types,
    )
end

function _passes_complexity_filter(types, mode::AbstractString)
    mode == "all" && return true
    mode == "identified" && return _is_identified(types)
    mode == "unidentified" && return !_is_identified(types)
    mode == "non-unknot" && return _is_non_unknot(types)
    error("unknown complexity knot filter '$mode'")
end

# ── PCA fit and scoring ────────────────────────────────────────────────────

"""
    fit_complexity_pca(catalogs; features, knot_filter, orientation)

Fit the PC1 complexity model on one or more catalog files.  Rows are
deduplicated by seed within each catalog.  The sign is fixed so the summed
loadings are positive; `orientation = :positive` then uses absolute loadings
as score weights.
"""
function fit_complexity_pca(
        catalogs::Vector{String};
        features::Vector{String} = COMPLEXITY_FEATURES,
        knot_filter::String = "identified",
        orientation::Symbol = :positive,
        min_training_rows::Int = 3,
    )
    isempty(catalogs) && error("complexity PCA needs at least one catalog")

    records = Any[]
    for catalog in catalogs
        append!(records, _catalog_records(catalog, features))
    end
    isempty(records) && error("no finite catalog rows for complexity PCA")

    training = filter(
        record -> _passes_complexity_filter(record["knot_types"], knot_filter),
        records,
    )
    if length(training) < min_training_rows
        @warn "only $(length(training)) row(s) matched knot filter '$knot_filter'; using all $(length(records)) row(s)"
        training = records
    end

    x = Matrix{Float64}(undef, length(training), length(features))
    for (row, record) in enumerate(training)
        x[row, :] = record["values"]
    end

    mean_vec = vec(mean(x, dims = 1))
    std_vec = vec(std(x, dims = 1; corrected = false))
    std_vec[std_vec .== 0.0] .= 1.0

    standardized = (x .- mean_vec') ./ std_vec'
    decomposition = svd(permutedims(standardized))
    loadings = decomposition.U[:, 1]
    variances = decomposition.S .^ 2
    explained = isempty(variances) ? 0.0 : variances[1] / sum(variances)

    signed_loadings = copy(loadings)
    orientation_name = "sum-positive"
    if sum(signed_loadings) < 0.0
        signed_loadings .*= -1.0
    elseif abs(sum(signed_loadings)) <= 1e-12
        first_nonzero = findfirst(value -> abs(value) > 1e-12, signed_loadings)
        if first_nonzero !== nothing && signed_loadings[first_nonzero] < 0.0
            signed_loadings .*= -1.0
            orientation_name = "first-nonzero-positive"
        end
    end

    weights = orientation == :signed ? copy(signed_loadings) : abs.(signed_loadings)

    placeholder = ComplexityModel(
        features, mean_vec, std_vec, weights, signed_loadings, explained,
        orientation_name, string(orientation), 0.0, 0.0,
    )

    # Fix the normalization range on the training rows only.  Unidentified
    # knot types do not contribute to the model or to its score scale.
    all_raw, _ = _complexity_components(placeholder, _record_matrix(training))
    raw_min = isempty(all_raw) ? 0.0 : minimum(all_raw)
    raw_max = isempty(all_raw) ? 0.0 : maximum(all_raw)

    model = ComplexityModel(
        features, mean_vec, std_vec, weights, signed_loadings, explained,
        orientation_name, string(orientation), raw_min, raw_max,
    )
    return model, records
end

function _record_matrix(records)
    isempty(records) && return Matrix{Float64}(undef, 0, 0)
    features = records[1]["values"]
    matrix = Matrix{Float64}(undef, length(records), length(features))
    for (row, record) in enumerate(records)
        matrix[row, :] = record["values"]
    end
    return matrix
end

function _complexity_components(model::ComplexityModel, x::AbstractMatrix)
    isempty(x) && return Float64[], Float64[]
    standardized = (x .- model.mean') ./ model.std'
    return standardized * model.weights, standardized * model.signed_loadings
end

function score_complexity(model::ComplexityModel, x::AbstractMatrix)
    raw, signed = _complexity_components(model, x)
    isempty(raw) && return ComplexityScores(Float64[], Float64[], Float64[], Float64[])

    span = model.raw_max - model.raw_min
    normalized = span > 0.0 ? (raw .- model.raw_min) ./ span : zeros(length(raw))

    order = sortperm(raw)
    ranks = zeros(Float64, length(raw))
    for (rank, index) in enumerate(order)
        ranks[index] = rank - 1
    end
    percentile = length(raw) > 1 ? ranks ./ (length(raw) - 1) : zeros(length(raw))

    return ComplexityScores(raw, signed, normalized, percentile)
end

# ── Output writers ─────────────────────────────────────────────────────────

function write_complexity_model(model::ComplexityModel, out_dir::AbstractString)
    mkpath(out_dir)
    rows = [
        [
            model.features[index],
            _fmt(model.mean[index]),
            _fmt(model.std[index]),
            _fmt(model.weights[index]),
            _fmt(model.signed_loadings[index]),
        ]
        for index in eachindex(model.features)
    ]
    _write_table(
        joinpath(out_dir, "complexity_model.tsv"),
        ["feature", "mean", "std", "pc1_loading", "pc1_loading_signed"],
        rows,
    )

    score_formula = model.orientation_mode == "signed" ?
        "pca_pc1_signed" : "weighted_sum_abs_pc1_loadings"
    summary = [
        "features=$(join(model.features, ','))",
        "pc1_explained_variance=$(_fmt(model.explained))",
        "orientation=$(model.orientation)",
        "orientation_mode=$(model.orientation_mode)",
        "score=$score_formula",
        "raw_min=$(_fmt(model.raw_min))",
        "raw_max=$(_fmt(model.raw_max))",
    ]
    open(joinpath(out_dir, "complexity_summary.txt"), "w") do io
        for line in summary
            println(io, line)
        end
    end
    return nothing
end

function _defined_complexity_records(records)
    return [record for record in records
            if any(is_defined_knot_type, record["knot_types"])]
end

function write_complexity_scores(
        model::ComplexityModel,
        records,
        out_dir::AbstractString;
        defined_only::Bool = true,
    )
    defined_only && (records = _defined_complexity_records(records))
    scores = score_complexity(model, _record_matrix(records))
    header = [
        "seed", "knotType", "complexity_pc1", "complexity_pc1_signed",
        "complexity_score", "complexity_percentile", model.features...,
    ]
    rows = [
        [
            record["seed"],
            join(sort(collect(record["knot_types"])), "|"),
            _fmt(scores.pc1[index]),
            _fmt(scores.pc1_signed[index]),
            _fmt(scores.score[index]),
            _fmt(scores.percentile[index]),
            [_fmt(value) for value in record["values"]]...,
        ]
        for (index, record) in enumerate(records)
    ]
    _write_table(joinpath(out_dir, "complexity_scores.tsv"), header, rows)
    return scores
end

# ── Augmenting existing tables ─────────────────────────────────────────────

function _seed_score_map(model::ComplexityModel, data, header)
    seed_index = _require_column(header, "seed")
    feature_indices = [_require_column(header, feature) for feature in model.features]

    seeds = String[]
    values = Vector{Float64}[]
    for row in 1:size(data, 1)
        seed = string(data[row, seed_index])
        isempty(seed) && continue
        x = _to_f64.(data[row, feature_indices])
        all(isfinite, x) || continue
        push!(seeds, seed)
        push!(values, x)
    end

    matrix = isempty(values) ? Matrix{Float64}(undef, 0, length(model.features)) :
        reduce(vcat, permutedims.(values))
    scores = score_complexity(model, matrix)
    return Dict(seed => index for (index, seed) in enumerate(seeds)), scores
end

function seed_score_map(model::ComplexityModel, records)
    scores = score_complexity(model, _record_matrix(records))
    return Dict(
        record["seed"] => (
            pc1 = scores.pc1[index],
            pc1_signed = scores.pc1_signed[index],
            score = scores.score[index],
            percentile = scores.percentile[index],
        )
        for (index, record) in enumerate(records)
    )
end

function augment_table_with_complexity(
        input_path::AbstractString,
        output_path::AbstractString,
        model::ComplexityModel;
        score_map::Union{Nothing, Dict} = nothing,
    )
    isfile(input_path) || return nothing
    data, header = _read_table(input_path)
    seed_index = _require_column(header, "seed")

    if score_map === nothing
        row_map, scores = _seed_score_map(model, data, header)
        score_map = Dict(
            seed => (
                pc1 = scores.pc1[index],
                pc1_signed = scores.pc1_signed[index],
                score = scores.score[index],
                percentile = scores.percentile[index],
            )
            for (seed, index) in row_map
        )
    end

    extra_columns = [
        "complexity_pc1", "complexity_pc1_signed",
        "complexity_score", "complexity_percentile",
    ]
    output_header = copy(header)
    for name in extra_columns
        name in output_header || push!(output_header, name)
    end

    rows = Vector{Vector{String}}()
    for row in 1:size(data, 1)
        entry = get(score_map, string(data[row, seed_index]), nothing)
        output_row = String[string(value) for value in data[row, :]]
        if entry === nothing
            append!(output_row, ["", "", "", ""])
        else
            append!(output_row, [
                _fmt(entry.pc1), _fmt(entry.pc1_signed),
                _fmt(entry.score), _fmt(entry.percentile),
            ])
        end
        push!(rows, output_row)
    end

    _write_table(output_path, output_header, rows)
    return output_path
end

# ── Stage scoring ──────────────────────────────────────────────────────────

"""
    score_stage_metrics(stage_file, model, output_file)

Score the rows in a cut-stage metrics file.  Stage 0 is the initial tangle;
later stages are scored with the same model.  The output adds
`complexity_score`, `delta_threading`, and `delta_complexity` (before − after).
Rows with an undefined `knotType` are skipped by default.
"""
function score_stage_metrics(
        stage_file::AbstractString,
        model::ComplexityModel,
        output_file::Union{Nothing, AbstractString} = nothing;
        defined_only::Bool = true,
    )
    isfile(stage_file) || return nothing
    data, header = _read_table(stage_file)
    isempty(data) && return nothing

    seed_index = _require_column(header, "seed")
    strategy_index = _require_column(header, "strategy")
    trial_index = _require_column(header, "trial")
    stage_index = _require_column(header, "stage")
    threading_index = _require_column(header, "threading")
    cuts_index = _require_column(header, "cuts")
    resolved_index = _require_column(header, "resolved")
    knot_index = _col_index(header, "knotType")
    feature_indices = [_require_column(header, feature) for feature in model.features]

    x = Matrix{Float64}(undef, size(data, 1), length(model.features))
    for row in 1:size(data, 1)
        x[row, :] = _to_f64.(data[row, feature_indices])
    end
    scores = score_complexity(model, x)

    keep_row(row) = !defined_only || knot_index === nothing ||
        is_defined_knot_type(data[row, knot_index])

    initial_score = Dict{String, Float64}()
    initial_score_signed = Dict{String, Float64}()
    initial_threading = Dict{String, Float64}()
    for row in 1:size(data, 1)
        keep_row(row) || continue
        _to_f64(data[row, stage_index]) == 0.0 || continue
        seed = string(data[row, seed_index])
        initial_score[seed] = scores.pc1[row]
        initial_score_signed[seed] = scores.pc1_signed[row]
        initial_threading[seed] = _to_f64(data[row, threading_index])
    end

    rows = Vector{Vector{String}}()
    for row in 1:size(data, 1)
        keep_row(row) || continue
        seed = string(data[row, seed_index])
        threading = _to_f64(data[row, threading_index])
        before_score = get(initial_score, seed, NaN)
        before_score_signed = get(initial_score_signed, seed, NaN)
        before_threading = get(initial_threading, seed, NaN)

        stage = _to_f64(data[row, stage_index])
        cuts = _to_f64(data[row, cuts_index])
        resolved = _to_f64(data[row, resolved_index])

        push!(rows, [
            seed,
            knot_index === nothing ? "" : string(data[row, knot_index]),
            string(data[row, strategy_index]),
            string(data[row, trial_index]),
            isfinite(stage) ? string(Int(round(stage))) : "",
            isfinite(cuts) ? string(Int(round(cuts))) : "",
            isfinite(threading) ? string(Int(round(threading))) : "",
            isfinite(resolved) ? string(Int(round(resolved))) : "",
            _fmt(scores.pc1[row]),
            _fmt(scores.pc1_signed[row]),
            _fmt(isfinite(before_threading) ? before_threading - threading : NaN),
            _fmt(isfinite(before_score) ? before_score - scores.pc1[row] : NaN),
            _fmt(isfinite(before_score_signed) ? before_score_signed - scores.pc1_signed[row] : NaN),
        ])
    end

    output_header = [
        "seed", "knotType", "strategy", "trial", "stage", "cuts", "threading",
        "resolved", "complexity_pc1", "complexity_pc1_signed", "delta_threading",
        "delta_complexity", "delta_complexity_signed",
    ]
    output_file === nothing || _write_table(output_file, output_header, rows)
    return (header = output_header, rows = rows)
end

# ── Comparison data for plots ──────────────────────────────────────────────

# Seeds kept for complexity comparisons: they must start with threading and
# have a defined knot type, matching the cut plots and paired tests.
function _analysis_stage_seeds(scored)
    seed_index = _col_index(scored.header, "seed")
    stage_index = _col_index(scored.header, "stage")
    threading_index = _col_index(scored.header, "threading")
    knot_index = _col_index(scored.header, "knotType")
    (seed_index === nothing || stage_index === nothing ||
        threading_index === nothing) && return nothing

    initial = Dict{String, Float64}()
    defined = Dict{String, Bool}()
    for row in scored.rows
        seed = row[seed_index]
        if knot_index !== nothing
            defined[seed] = get(defined, seed, false) ||
                is_defined_knot_type(row[knot_index])
        end
        _to_f64(row[stage_index]) == 0.0 || continue
        threading = _to_f64(row[threading_index])
        isfinite(threading) || continue
        initial[seed] = threading
    end
    isempty(initial) && return nothing
    return Set(
        seed for (seed, threading) in initial
        if threading > 0 && (knot_index === nothing || get(defined, seed, false))
    )
end

"""
    stage_comparison_deltas(stage_file, model)

Return per-seed, per-strategy single-cut complexity deltas.  Random strategies
are averaged over their trials.
"""
function stage_comparison_deltas(
        stage_file::AbstractString,
        model::ComplexityModel,
    )
    scored = score_stage_metrics(stage_file, model)
    scored === nothing && return nothing
    threaded = _analysis_stage_seeds(scored)

    seed_index = _col_index(scored.header, "seed")
    strategy_index = _col_index(scored.header, "strategy")
    stage_index = _col_index(scored.header, "stage")
    threading_index = _col_index(scored.header, "delta_threading")
    complexity_index = _col_index(scored.header, "delta_complexity")
    signed_index = _col_index(scored.header, "delta_complexity_signed")

    grouped = Dict{Tuple{String, String}, Vector{Tuple{Float64, Float64}}}()
    for row in scored.rows
        _to_f64(row[stage_index]) == 1.0 || continue
        seed = row[seed_index]
        threaded === nothing || seed in threaded || continue
        label = get(CUT_STRATEGY_LABEL_MAP, row[strategy_index], row[strategy_index])
        threading = _to_f64(row[threading_index])
        complexity = _to_f64(row[complexity_index])
        signed = _to_f64(row[signed_index])
        (isfinite(threading) && isfinite(complexity) && isfinite(signed)) || continue
        push!(get!(grouped, (seed, label), Tuple{Float64, Float64}[]),
            (complexity, signed))
    end

    seeds = sort(unique(first.(keys(grouped))))
    strategies = [label for label in CUT_STRATEGY_LABELS
                  if any(key -> key[2] == label, keys(grouped))]
    complexity, signed = _per_seed_means(grouped, seeds, strategies, 1, 2)

    return (
        seeds = seeds,
        strategies = strategies,
        complexity = complexity,
        complexity_signed = signed,
    )
end

"""
    resolved_complexity_reduction(stage_file, model)

Return per-seed, per-strategy complexity reduction at the final cuts-to-zero
stage, considering only trials that actually resolved threading to zero.
"""
function resolved_complexity_reduction(
        stage_file::AbstractString,
        model::ComplexityModel,
    )
    scored = score_stage_metrics(stage_file, model)
    scored === nothing && return nothing
    threaded = _analysis_stage_seeds(scored)

    seed_index = _require_column(scored.header, "seed")
    strategy_index = _require_column(scored.header, "strategy")
    trial_index = _require_column(scored.header, "trial")
    stage_index = _require_column(scored.header, "stage")
    resolved_index = _require_column(scored.header, "resolved")
    complexity_index = _require_column(scored.header, "delta_complexity")
    signed_index = _require_column(scored.header, "delta_complexity_signed")

    # Keep the last stage of each seed/strategy/trial sequence.
    final_stage = Dict{Tuple{String, String, String}, NamedTuple}()
    for row in scored.rows
        seed = row[seed_index]
        threaded === nothing || seed in threaded || continue
        label = get(CUT_STRATEGY_LABEL_MAP, row[strategy_index], row[strategy_index])
        key = (seed, label, row[trial_index])
        stage = _to_f64(row[stage_index])
        if !haskey(final_stage, key) || stage > final_stage[key].stage
            final_stage[key] = (
                stage = stage,
                resolved = _to_f64(row[resolved_index]),
                positive = _to_f64(row[complexity_index]),
                signed = _to_f64(row[signed_index]),
            )
        end
    end

    grouped = Dict{Tuple{String, String}, Vector{Tuple{Float64, Float64}}}()
    for ((seed, label, _), row) in final_stage
        row.resolved == 1.0 || continue
        (isfinite(row.positive) && isfinite(row.signed)) || continue
        push!(get!(grouped, (seed, label), Tuple{Float64, Float64}[]),
            (row.positive, row.signed))
    end

    seeds = sort(unique(first.(keys(grouped))))
    available = [label for label in CUT_STRATEGY_LABELS
                 if any(key -> key[2] == label, keys(grouped))]
    positive, signed = _per_seed_means(grouped, seeds, available, 1, 2)

    # Drop strategies without any resolved observation.
    keep = [label for label in available if any(isfinite, positive[label])]
    return (
        seeds = seeds,
        strategies = keep,
        complexity = Dict(label => positive[label] for label in keep),
        complexity_signed = Dict(label => signed[label] for label in keep),
    )
end

function _per_seed_means(grouped, seeds, strategies, positive_index::Int, signed_index::Int)
    positive = Dict(label => fill(NaN, length(seeds)) for label in strategies)
    signed = Dict(label => fill(NaN, length(seeds)) for label in strategies)
    for label in strategies
        for (index, seed) in enumerate(seeds)
            entries = get(grouped, (seed, label), nothing)
            entries === nothing && continue
            positive[label][index] = mean(getindex.(entries, positive_index))
            signed[label][index] = mean(getindex.(entries, signed_index))
        end
    end
    return positive, signed
end

# ── Complexity violins and paired tests ────────────────────────────────────

_ensure_model(settings::Settings, model) = model === nothing ?
    first(fit_complexity_pca(
        [settings.metrics_file];
        features = settings.complexity_features,
        knot_filter = settings.complexity_knot_filter,
        orientation = settings.complexity_orientation,
        min_training_rows = settings.complexity_min_training_rows,
    )) : model

"""
    plot_complexity_reduction_first_cut([stage_file]; settings, model)

Violin comparison of the signed PC1 complexity reduction after the first cut
for every strategy.  Significant paired Wilcoxon comparisons are marked with
red asterisks when the strategy is better and green asterisks when the
heuristic is better, with the matched-pairs rank-biserial effect size printed
above the marker.
"""
function plot_complexity_reduction_first_cut(
        stage_file::AbstractString = DEFAULT_SETTINGS.cut_eval_stage_file;
        settings::Settings = DEFAULT_SETTINGS,
        model::Union{Nothing, ComplexityModel} = nothing,
    )
    model = _ensure_model(settings, model)
    comparison = stage_comparison_deltas(stage_file, model)
    comparison === nothing && error("no stage metrics found in $stage_file")
    isempty(comparison.strategies) && error("no cut strategies found in $stage_file")

    outpath = joinpath(settings.figures_dir, "cut_eval", "cut_eval_comparison_complexity.pdf")
    return _plot_metric_violin(
        comparison.strategies, comparison.complexity_signed;
        outpath = outpath,
        title = "Complexity reduction after the first cut (signed PC1)",
        ylabel = "Signed PC1 complexity reduction",
        metric = "complexity_reduction_first_cut_signed",
        comparison = comparison,
    )
end

"""
    plot_complexity_reduction_resolved([stage_file]; settings, model)

Violin comparison of the signed PC1 complexity reduction after the
cuts-to-zero sequence has resolved threading to zero.
"""
function plot_complexity_reduction_resolved(
        stage_file::AbstractString = DEFAULT_SETTINGS.cuts_to_zero_stage_file;
        settings::Settings = DEFAULT_SETTINGS,
        model::Union{Nothing, ComplexityModel} = nothing,
    )
    model = _ensure_model(settings, model)
    comparison = resolved_complexity_reduction(stage_file, model)
    comparison === nothing && error("no stage metrics found in $stage_file")
    isempty(comparison.strategies) && error("no resolved cut strategies found in $stage_file")

    outpath = joinpath(settings.figures_dir, "cut_eval", "complexity_reduction_resolved.pdf")
    return _plot_metric_violin(
        comparison.strategies, comparison.complexity_signed;
        outpath = outpath,
        title = "Complexity reduction after complete threading resolution (signed PC1)",
        ylabel = "Signed PC1 complexity reduction",
        metric = "complexity_reduction_resolved_signed",
        comparison = comparison,
    )
end

"""
    complexity_reduction_first_cut_paired_tests([stage_file]; settings, model)

Paired Wilcoxon comparisons of first-cut complexity reduction against the
heuristic strategy.  Runs separately for the positive and signed PC1 scores.
"""
function complexity_reduction_first_cut_paired_tests(
        stage_file::AbstractString = DEFAULT_SETTINGS.cut_eval_stage_file;
        settings::Settings = DEFAULT_SETTINGS,
        model::Union{Nothing, ComplexityModel} = nothing,
    )
    isfile(stage_file) || return nothing
    model = _ensure_model(settings, model)
    comparison = stage_comparison_deltas(stage_file, model)
    comparison === nothing && return nothing

    rows = vcat(
        _paired_rows_for_metric(
            comparison, comparison.complexity,
            "complexity_reduction_first_cut_positive"),
        _paired_rows_for_metric(
            comparison, comparison.complexity_signed,
            "complexity_reduction_first_cut_signed"),
    )
    if isempty(rows)
        @warn "not enough paired first-cut complexity observations"
        return rows
    end

    outpath = joinpath(settings.analysis_dir, "complexity_reduction_first_cut_paired_tests.tsv")
    _write_paired_report(rows, outpath)
    _print_paired_report(rows, "Complexity reduction after first cut vs heuristic")
    println("\nSaved paired test table to $outpath")
    return rows
end

"""
    complexity_reduction_resolved_paired_tests([stage_file]; settings, model)

Paired Wilcoxon comparisons of complexity reduction after complete threading
resolution against the heuristic strategy.
"""
function complexity_reduction_resolved_paired_tests(
        stage_file::AbstractString = DEFAULT_SETTINGS.cuts_to_zero_stage_file;
        settings::Settings = DEFAULT_SETTINGS,
        model::Union{Nothing, ComplexityModel} = nothing,
    )
    isfile(stage_file) || return nothing
    model = _ensure_model(settings, model)
    comparison = resolved_complexity_reduction(stage_file, model)
    comparison === nothing && return nothing

    rows = vcat(
        _paired_rows_for_metric(
            comparison, comparison.complexity,
            "complexity_reduction_resolved_positive"),
        _paired_rows_for_metric(
            comparison, comparison.complexity_signed,
            "complexity_reduction_resolved_signed"),
    )
    if isempty(rows)
        @warn "not enough paired resolved complexity observations"
        return rows
    end

    outpath = joinpath(settings.analysis_dir, "complexity_reduction_resolved_paired_tests.tsv")
    _write_paired_report(rows, outpath)
    _print_paired_report(rows, "Complexity reduction after threading resolution vs heuristic")
    println("\nSaved paired test table to $outpath")
    return rows
end

# ── High-level entry point ─────────────────────────────────────────────────

"""
    analyze_complexity([settings]; catalogs, out_dir, augment)

Fit the PCA complexity model, write model/scores, augment catalog and
evaluation tables, and score the cut-stage metrics.  Returns the model.
"""
function analyze_complexity(
        settings::Settings = DEFAULT_SETTINGS;
        catalogs::Union{Nothing, Vector{String}} = nothing,
        out_dir::Union{Nothing, String} = nothing,
        augment::Bool = true,
    )
    catalog_paths = catalogs === nothing ? [settings.metrics_file] : catalogs
    output_dir = out_dir === nothing ? settings.complexity_dir : out_dir

    model, records = fit_complexity_pca(
        catalog_paths;
        features = settings.complexity_features,
        knot_filter = settings.complexity_knot_filter,
        orientation = settings.complexity_orientation,
        min_training_rows = settings.complexity_min_training_rows,
    )

    write_complexity_model(model, output_dir)
    scored_records = _defined_complexity_records(records)
    write_complexity_scores(model, scored_records, output_dir)
    score_map = seed_score_map(model, scored_records)

    if augment
        for catalog in catalog_paths
            base = splitext(catalog)
            augment_table_with_complexity(
                catalog,
                base[1] * "_with_complexity" * base[2],
                model;
                score_map = score_map,
            )
        end
        augment_table_with_complexity(
            settings.cut_eval_file,
            joinpath(dirname(settings.cut_eval_file), "cut_eval_with_complexity.tsv"),
            model;
            score_map = score_map,
        )
        augment_table_with_complexity(
            settings.cuts_file,
            joinpath(dirname(settings.cuts_file), "cuts_to_zero_with_complexity.tsv"),
            model;
            score_map = score_map,
        )
    end

    score_stage_metrics(
        settings.cut_eval_stage_file,
        model,
        joinpath(output_dir, "cut_eval_stage_complexity.tsv"),
    )
    score_stage_metrics(
        settings.cuts_to_zero_stage_file,
        model,
        joinpath(output_dir, "cuts_to_zero_stage_complexity.tsv"),
    )

    println("Complexity model written to $output_dir")
    println("PC1 explained variance: $(round(model.explained * 100.0, digits=2))%")
    for index in eachindex(model.features)
        println("  $(model.features[index]): $(round(model.weights[index], digits=6))")
    end
    return model
end
