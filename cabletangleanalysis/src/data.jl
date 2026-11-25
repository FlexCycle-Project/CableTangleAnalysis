# ── Data loaders ───────────────────────────────────────────────────────────
#
# All paths come from a Settings struct (defaults to DEFAULT_SETTINGS).  See
# CableTangleAnalysis.jl for the struct definition.

using DelimitedFiles

# Knot labels that do not identify a knot type.  Tangles with these labels
# have no well-defined outside-knot geometry and are excluded from the metric
# and strategy analyses (but still appear in the knot distribution).
const UNDEFINED_KNOT_LABELS = Set(["", "none", "unknown", "??"])

is_defined_knot_type(label) =
    !(lowercase(strip(string(label))) in UNDEFINED_KNOT_LABELS)

# Shared grouping rule for knot labels:
#   undefined → "Undefined", UN → "Trivial Knot", composites → "Mult.",
#   "3_1" → "3_i" (grouped by crossing number).
function _grouped_knot_label(label)
    label = string(label)
    if !is_defined_knot_type(label)
        return "Undefined"
    elseif label == "UN"
        return "Trivial Knot"
    elseif occursin("#", label)
        return "Mult."
    elseif occursin('_', label)
        return string(split(label, '_')[1]) * "_i"
    end
    return label
end

"""
    load_data([settings::Settings]; defined_only = false)
        -> (data::Matrix, header::Vector{String})

Load the full catalog from `settings.metrics_file`.  With
`defined_only = true`, rows whose `knotType` is undefined (`none`, `unknown`,
`??`, or empty) are removed.
"""
function load_data(
        settings::Settings = DEFAULT_SETTINGS;
        defined_only::Bool = false,
    )
    data, header = readdlm(settings.metrics_file, header = true)
    if defined_only
        index = findfirst(==("knotType"), vec(header))
        if index !== nothing
            keep = [is_defined_knot_type(data[row, index]) for row in 1:size(data, 1)]
            data = data[keep, :]
        end
    end
    return data, header
end

"""
    load_data_and_group_knots([settings::Settings]; defined_only = false)
        -> (data::Matrix, header::Vector{String})

Like [`load_data`](@ref) but also replaces the first column with grouped knot
labels (see `_grouped_knot_label`).
"""
function load_data_and_group_knots(
        settings::Settings = DEFAULT_SETTINGS;
        defined_only::Bool = false,
    )
    data, header = load_data(settings; defined_only = defined_only)
    for row in 1:size(data, 1)
        data[row, 1] = _grouped_knot_label(data[row, 1])
    end
    return data, header
end

"""
    clear_data_of_unwanted_knots([settings::Settings])

Write the catalog without `"Undefined"` and `"Mult."` groups to
`data_dir/cleared_metrics_2.tsv`.  This file is the input for `run_pca()`.
"""
function clear_data_of_unwanted_knots(settings::Settings = DEFAULT_SETTINGS)
    data, header = load_data_and_group_knots(settings)

    keep = [string(data[row, 1]) ∉ ("Undefined", "Mult.", "none")
            for row in 1:size(data, 1)]
    cleaned = data[keep, :]

    outpath = joinpath(settings.data_dir, "cleared_metrics_2.tsv")
    open(outpath, "w") do io
        writedlm(io, header, '\t')     # header row (1×N from readdlm)
        writedlm(io, cleaned, '\t')
    end
    println("Cleaned data written to $outpath ($(size(cleaned, 1)) rows)")
    return nothing
end
