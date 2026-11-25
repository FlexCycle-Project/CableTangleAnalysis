# ── Shared plotting utilities ──────────────────────────────────────────────

using CairoMakie
using Statistics

"""
    _violin_top(values) -> Float64

Upper edge of the density that `violin!` draws for `values`.  Makie evaluates
the KDE on a grid that extends four Silverman bandwidths past the data range,
so annotations must clear this edge rather than just the largest observation.
"""
function _violin_top(values)
    finite = values[isfinite.(values)]
    isempty(finite) && return NaN
    _, hi = extrema(finite)
    n = length(finite)
    bandwidth = if n <= 1
        0.9
    else
        var_width = std(finite)
        width = min(var_width, (quantile(finite, 0.75) - quantile(finite, 0.25)) / 1.34)
        if width == 0.0
            width = var_width == 0.0 ? 1.0 : var_width
        end
        0.9 * width * n^(-0.2)
    end
    return hi + 4.0 * bandwidth
end

"""
    beautify_knot_labels(labels) -> Vector{RichText}

Convert knot labels (e.g. `"3_1"`, `"10_124"`, `"11a_1"`) into
subscript-rendered rich text for use in CairoMakie axis ticks.
Replaces `'#'` with `'-'`.
"""
function beautify_knot_labels(labels)
    render_label(label) = begin
        parts = split(replace(label, "#" => "-"), "_")
        pieces = Any[rich(string(parts[1]))]
        for part in parts[2:end]
            m = match(r"^(\d+|i)(.*)$", part)
            if m === nothing
                push!(pieces, rich(string(part)))
                continue
            end
            push!(pieces, subscript(string(m.captures[1])))
            isempty(m.captures[2]) || push!(pieces, rich(string(m.captures[2])))
        end
        rich(pieces...)
    end

    return render_label.(labels)
end

"""
    _strategy_colors(n; palette_name = :managua10) -> Vector

One categorical colour per strategy, in strategy order.
"""
function _strategy_colors(n; palette_name = :managua10)
    palette = cgrad(palette_name, max(n, 1), categorical = true)
    return [palette[i] for i in 1:n]
end
