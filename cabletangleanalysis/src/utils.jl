# ── Shared table and column helpers ────────────────────────────────────────
#
# Small utilities used by every analysis file.  Keeping them in one place
# avoids the slightly different copies that used to live in cuts.jl and
# complexity.jl.

using DelimitedFiles

_col_index(header, name) = findfirst(==(name), vec(header))

function _require_column(header, name)
    index = _col_index(header, name)
    index === nothing && error("missing column '$name'")
    return index
end

_to_f64(value) = value isa Number ? Float64(value) : something(
    tryparse(Float64, string(value)),
    NaN,
)

function _read_table(path::AbstractString)
    data, header = readdlm(path, '\t', header = true)
    return data, string.(vec(header))
end

function _write_table(path::AbstractString, header, rows)
    directory = dirname(path)
    isempty(directory) || mkpath(directory)
    open(path, "w") do io
        println(io, join(string.(header), '\t'))
        for row in rows
            println(io, join(string.(row), '\t'))
        end
    end
end

function _numeric_col(data, header, name)
    index = _col_index(header, name)
    index === nothing && return fill(NaN, size(data, 1))
    return _to_f64.(data[:, index])
end

function _string_col(data, header, name)
    index = _col_index(header, name)
    index === nothing && return fill("", size(data, 1))
    return string.(data[:, index])
end

# Rows with an unidentified knot type have no well-defined outside-knot
# geometry.  If the column is absent (older files), keep every row.
function _defined_knot_mask(data, header)
    index = _col_index(header, "knotType")
    index === nothing && return trues(size(data, 1))
    return is_defined_knot_type.(_string_col(data, header, "knotType"))
end

function _finite_mask(arrays...)
    mask = trues(length(arrays[1]))
    for array in arrays
        mask .&= isfinite.(array)
    end
    return mask
end

function _fmt(value; digits = 8)
    isfinite(value) || return ""
    return string(round(value, digits = digits))
end

# Effect sizes are plotted with two decimals; clamp tiny negatives so the
# annotation never reads "-0.0".
function _fmt_effect(value; digits = 2)
    isfinite(value) || return ""
    rounded = round(value, digits = digits)
    abs(rounded) < 10.0^-digits && (rounded = 0.0)
    return string(rounded)
end
