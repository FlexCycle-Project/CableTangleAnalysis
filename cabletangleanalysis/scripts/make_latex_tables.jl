#!/usr/bin/env julia
#
# Generate concise booktabs LaTeX tables from CableTangleAnalysis outputs.
#
# Usage:
#   julia --project=<repo>/cabletangleanalysis \
#       <repo>/cabletangleanalysis/scripts/make_latex_tables.jl \
#       <analysis-dir> <latex-out-dir>
#
# <analysis-dir> must contain a tables/ subdirectory with pca_results.tsv and
# the complexity paired-test tables.

using DelimitedFiles

const STRATEGY_ORDER = [
    "Heuristic", "Random edge", "Random hoop", "Biggest hoop",
    "Smallest hoop", "Max curvature", "Oracle",
]

const METRIC_LABELS = Dict(
    "segmentCount" => "Segment count",
    "selfCollisionCount" => "Self-collision count",
    "closePairCount" => "Close-pair count",
    "totalCurvature" => "Total curvature",
    "outsideKnotCurvature" => "Outside-knot curvature",
)

function read_table(path::AbstractString)
    data, header = readdlm(path, '\t', header = true)
    return data, string.(vec(header))
end

function column_index(header, name)
    return findfirst(==(name), header)
end

function value_at(data, header, row, name)
    index = column_index(header, name)
    index === nothing && return ""
    return string(data[row, index])
end

function latex_escape(text::AbstractString)
    return replace(
        string(text),
        "\\" => "\\textbackslash{}",
        "_" => "\\_",
        "%" => "\\%",
        "&" => "\\&",
        "#" => "\\#",
    )
end

function fmt_number(value; digits = 3)
    parsed = tryparse(Float64, string(value))
    parsed === nothing && return "--"
    isfinite(parsed) || return "--"
    abs(parsed) < 0.5 * 10.0^(-digits) && (parsed = 0.0)
    return string(round(parsed, digits = digits))
end

function fmt_p(value)
    parsed = tryparse(Float64, string(value))
    parsed === nothing && return "--"
    isfinite(parsed) || return "--"
    parsed <= 0.0 && return "\$<10^{-300}\$"
    parsed >= 0.001 && return "\$" * string(round(parsed, digits = 3)) * "\$"
    exponent = floor(Int, log10(parsed))
    mantissa = parsed / 10.0^exponent
    return "\$$(round(mantissa, digits = 1))\\times10^{$(exponent)}\$"
end

function strategy_rank(strategy)
    index = findfirst(==(strategy), STRATEGY_ORDER)
    return index === nothing ? length(STRATEGY_ORDER) + 1 : index
end

function write_table(path::AbstractString, lines::Vector{String})
    mkpath(dirname(path))
    open(path, "w") do io
        for line in lines
            println(io, line)
        end
    end
    println("Wrote $path")
end

function write_pca_table(tables_dir::String, out_dir::String)
    input = joinpath(tables_dir, "pca_results.tsv")
    isfile(input) || (@warn "missing $input"; return)

    data, header = read_table(input)
    component_index = column_index(header, "component")
    explained_index = column_index(header, "explained_variance_pct")
    cumulative_index = column_index(header, "cumulative_pct")
    component_index === nothing && return

    metadata = Set([
        "component", "eigenvalue", "explained_variance",
        "explained_variance_pct", "cumulative_pct",
    ])
    metric_indices = [index for (index, name) in enumerate(header) if !(name in metadata)]
    metric_names = header[metric_indices]
    components = [string(data[row, component_index]) for row in 1:size(data, 1)]

    lines = String[
        "\\begin{table}[t]",
        "\\centering",
        "\\caption{PCA of the five cable-complexity metrics. Loadings are computed on standardized metrics; the sign of each component is arbitrary.}",
        "\\label{tab:pca}",
        "\\begin{tabular}{l" * repeat("r", length(components)) * "}",
        "\\toprule",
        "Metric & " * join(components, " & ") * " \\\\",
        "\\midrule",
    ]

    for (metric_column, metric_name) in zip(metric_indices, metric_names)
        label = get(METRIC_LABELS, metric_name, metric_name)
        values = [
            fmt_number(data[row, metric_column])
            for row in 1:size(data, 1)
        ]
        push!(lines, latex_escape(label) * " & " * join(values, " & ") * " \\\\")
    end

    push!(lines, "\\midrule")
    push!(
        lines,
        "Explained variance (\\%) & " *
        join([fmt_number(data[row, explained_index]) for row in 1:size(data, 1)], " & ") *
        " \\\\",
    )
    push!(
        lines,
        "Cumulative (\\%) & " *
        join([fmt_number(data[row, cumulative_index]) for row in 1:size(data, 1)], " & ") *
        " \\\\",
    )
    push!(lines, "\\bottomrule")
    push!(lines, "\\end{tabular}")
    push!(lines, "\\end{table}")

    write_table(joinpath(out_dir, "table_pca.tex"), lines)
end

function write_complexity_table(
        tables_dir::String,
        out_dir::String,
        filename::String,
        output_basename::String,
        metric_suffix::String,
        caption::String,
        label::String,
    )
    input = joinpath(tables_dir, filename)
    isfile(input) || (@warn "missing $input"; return)

    data, header = read_table(input)
    rows = NamedTuple[]

    for row in 1:size(data, 1)
        metric = value_at(data, header, row, "metric")
        endswith(metric, metric_suffix) || continue
        strategy = value_at(data, header, row, "strategy")
        push!(rows, (
            strategy = strategy,
            n = value_at(data, header, row, "n"),
            median = value_at(data, header, row, "median_diff"),
            ci_low = value_at(data, header, row, "ci_low"),
            ci_high = value_at(data, header, row, "ci_high"),
            p = value_at(data, header, row, "p_adjusted"),
            symmetry = value_at(data, header, row, "symmetry_ok"),
            better = value_at(data, header, row, "better"),
        ))
    end

    sort!(rows, by = row -> strategy_rank(row.strategy))

    lines = String[
        "\\begin{table}[t]",
        "\\centering",
        "\\caption{$caption Values are strategy \$-\$ heuristic; positive values mean more complexity reduction. \$p\$ is the Benjamini--Hochberg-adjusted Wilcoxon signed-rank \$p\$-value. Sym.\\ marks whether the symmetry diagnostic passed at \$\\alpha=0.05\$.}",
        "\\label{$label}",
        "\\begin{tabular}{lrrrcl}",
        "\\toprule",
        "Strategy & \$n\$ & Median \$\\Delta\$ [95\\% CI] & \$p\$ & Sym. & Better \\\\",
        "\\midrule",
    ]

    for row in rows
        symmetry = lowercase(row.symmetry)
        symmetry_text = symmetry == "true" ? "Y" : (symmetry == "false" ? "N" : "--")
        better = row.better == "tie" ? "Tie" : latex_escape(row.better)
        median_ci = "\$[$(fmt_number(row.ci_low)), $(fmt_number(row.ci_high))]\$"
        push!(
            lines,
            latex_escape(row.strategy) * " & " * row.n * " & " *
            "\$" * fmt_number(row.median) * "\$ " * median_ci * " & " *
            fmt_p(row.p) * " & " * symmetry_text * " & " * better * " \\\\",
        )
    end

    push!(lines, "\\bottomrule")
    push!(lines, "\\end{tabular}")
    push!(lines, "\\end{table}")

    write_table(joinpath(out_dir, output_basename * ".tex"), lines)
end

function main(args)
    length(args) >= 2 || error(
        "usage: make_latex_tables.jl <analysis-dir> <latex-out-dir>",
    )
    analysis_dir = abspath(args[1])
    out_dir = abspath(args[2])
    tables_dir = joinpath(analysis_dir, "tables")
    isdir(tables_dir) || error("missing tables directory: $tables_dir")

    mkpath(out_dir)
    write_pca_table(tables_dir, out_dir)
    write_complexity_table(
        tables_dir,
        out_dir,
        "complexity_reduction_first_cut_paired_tests.tsv",
        "table_complexity_first_cut",
        "signed",
        "Paired comparison of signed PC1 complexity reduction after the first cut.",
        "tab:complexity-first-cut",
    )
    write_complexity_table(
        tables_dir,
        out_dir,
        "complexity_reduction_resolved_paired_tests.tsv",
        "table_complexity_resolved",
        "signed",
        "Paired comparison of signed PC1 complexity reduction after complete threading resolution.",
        "tab:complexity-resolved",
    )

    write_table(
        joinpath(out_dir, "tables.tex"),
        String[
            "\\input{table_pca}",
            "\\input{table_complexity_first_cut}",
            "\\input{table_complexity_resolved}",
        ],
    )
end

main(ARGS)
