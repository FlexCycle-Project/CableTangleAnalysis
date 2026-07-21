import argparse
import csv
import html
import math
from collections import Counter
from pathlib import Path

import numpy as np


DEFAULT_FEATURES = [
    "totalHoopDensity",
    "threadedHoopDensity",
    "threadingDensity",
    "selfCollisionDensity",
    "closePairDensity",
    "totalCurvatureDensity",
    "outsideKnotCurvatureDensity",
]


def read_tsv(path):
    with open(path, newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def write_tsv(path, rows, fieldnames):
    with open(path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def parse_float(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return float("nan")


def load_matrix(rows, features):
    kept_rows = []
    values = []
    for row in rows:
        if row.get("knotType") == "none":
            continue

        feature_values = [parse_float(row.get(feature)) for feature in features]
        if not all(math.isfinite(value) for value in feature_values):
            continue

        kept_rows.append(row)
        values.append(feature_values)

    return kept_rows, np.asarray(values, dtype=float)


def run_pca(x):
    mean = x.mean(axis=0)
    std = x.std(axis=0)
    std[std == 0.0] = 1.0
    xs = (x - mean) / std

    _, singular_values, vt = np.linalg.svd(xs, full_matrices=False)
    variances = singular_values ** 2
    explained = variances / variances.sum()
    coords = xs @ vt.T
    return mean, std, vt, explained, coords


def write_coords(path, rows, coords, pc_count):
    fieldnames = ["knotType", "seed"] + [f"PC{i + 1}" for i in range(pc_count)]
    out_rows = []
    for row, coord in zip(rows, coords):
        out = {
            "knotType": row.get("knotType", ""),
            "seed": row.get("seed", ""),
        }
        for i in range(pc_count):
            out[f"PC{i + 1}"] = f"{coord[i]:.8f}"
        out_rows.append(out)
    write_tsv(path, out_rows, fieldnames)


def write_loadings(path, features, components, explained, pc_count):
    fieldnames = ["feature"] + [f"PC{i + 1}" for i in range(pc_count)]
    rows = []
    for feature_index, feature in enumerate(features):
        row = {"feature": feature}
        for pc_index in range(pc_count):
            row[f"PC{pc_index + 1}"] = f"{components[pc_index, feature_index]:.8f}"
        rows.append(row)

    explained_row = {"feature": "explainedVariance"}
    for pc_index in range(pc_count):
        explained_row[f"PC{pc_index + 1}"] = f"{explained[pc_index]:.8f}"
    rows.append(explained_row)
    write_tsv(path, rows, fieldnames)


def scale(value, src_min, src_max, dst_min, dst_max):
    if src_max <= src_min:
        return (dst_min + dst_max) * 0.5
    t = (value - src_min) / (src_max - src_min)
    return dst_min + t * (dst_max - dst_min)


def visual_spread(coords, strength):
    if strength <= 0.0:
        return coords

    spread = coords.copy()
    center = np.median(spread[:, :2], axis=0)
    xy = spread[:, :2] - center
    radius = np.linalg.norm(xy, axis=1)
    positive = radius[radius > 1e-12]
    if positive.size == 0:
        return coords

    reference = np.median(positive)
    reference = max(reference, 1e-12)
    # The factor is strongest near the dense center and fades out for far points.
    factor = 1.0 + strength * np.exp(-radius / reference)
    spread[:, :2] = center + xy * factor[:, None]
    return spread


def color_palette():
    return [
        "#2563eb", "#dc2626", "#16a34a", "#9333ea", "#ea580c",
        "#0891b2", "#be123c", "#4d7c0f", "#7c3aed", "#0f766e",
        "#ca8a04", "#c026d3", "#475569", "#65a30d", "#0284c7",
    ]


def write_html_plot(path, rows, coords, title, spread_strength):
    plot_coords = visual_spread(coords, spread_strength)
    pc1 = plot_coords[:, 0]
    pc2 = plot_coords[:, 1]
    pc1_min, pc1_max = float(pc1.min()), float(pc1.max())
    pc2_min, pc2_max = float(pc2.min()), float(pc2.max())

    knot_counts = Counter(row.get("knotType", "unknown") for row in rows)
    top_types = []
    for name, _ in knot_counts.most_common():
        if name == "??":
            continue
        top_types.append(name)
        if len(top_types) >= 12:
            break
    palette = color_palette()
    color_by_type = {name: palette[i % len(palette)] for i, name in enumerate(top_types)}
    color_by_type["UN"] = "#000000"
    other_color = "#94a3b8"

    width, height = 1200, 850
    margin_left, margin_right = 90, 260
    margin_top, margin_bottom = 60, 90
    plot_left = margin_left
    plot_right = width - margin_right
    plot_top = margin_top
    plot_bottom = height - margin_bottom

    def tick_values(vmin, vmax, count=5):
        if vmax <= vmin:
            return [vmin]
        return [vmin + (vmax - vmin) * i / (count - 1) for i in range(count)]

    svg = []
    svg.append(f'<svg viewBox="0 0 {width} {height}" width="{width}" height="{height}" role="img">')
    svg.append('<rect width="100%" height="100%" fill="#f8fafc"/>')
    svg.append(f'<rect x="{plot_left}" y="{plot_top}" width="{plot_right - plot_left}" height="{plot_bottom - plot_top}" fill="#ffffff" stroke="#cbd5e1"/>')

    for value in tick_values(pc1_min, pc1_max):
        x = scale(value, pc1_min, pc1_max, plot_left, plot_right)
        svg.append(f'<line x1="{x:.2f}" y1="{plot_top}" x2="{x:.2f}" y2="{plot_bottom}" stroke="#e2e8f0"/>')
        svg.append(f'<text x="{x:.2f}" y="{plot_bottom + 28}" font-size="14" text-anchor="middle" fill="#475569">{value:.2f}</text>')

    for value in tick_values(pc2_min, pc2_max):
        y = scale(value, pc2_min, pc2_max, plot_bottom, plot_top)
        svg.append(f'<line x1="{plot_left}" y1="{y:.2f}" x2="{plot_right}" y2="{y:.2f}" stroke="#e2e8f0"/>')
        svg.append(f'<text x="{plot_left - 14}" y="{y + 5:.2f}" font-size="14" text-anchor="end" fill="#475569">{value:.2f}</text>')

    svg.append(f'<line x1="{plot_left}" y1="{plot_bottom}" x2="{plot_right}" y2="{plot_bottom}" stroke="#334155" stroke-width="1.5"/>')
    svg.append(f'<line x1="{plot_left}" y1="{plot_top}" x2="{plot_left}" y2="{plot_bottom}" stroke="#334155" stroke-width="1.5"/>')
    svg.append(f'<text x="{(plot_left + plot_right) / 2:.2f}" y="{height - 28}" font-size="18" text-anchor="middle" fill="#0f172a">PC1</text>')
    svg.append(f'<text transform="translate(30 {(plot_top + plot_bottom) / 2:.2f}) rotate(-90)" font-size="18" text-anchor="middle" fill="#0f172a">PC2</text>')
    svg.append(f'<text x="{plot_left}" y="34" font-size="24" font-weight="700" fill="#0f172a">{html.escape(title)}</text>')
    svg.append(f'<text x="{plot_left}" y="55" font-size="14" fill="#64748b">Rows: {len(rows)}. Colored by KymoKnot type; uncommon types are gray.</text>')

    for row, coord, original_coord in zip(rows, plot_coords, coords):
        x = scale(float(coord[0]), pc1_min, pc1_max, plot_left, plot_right)
        y = scale(float(coord[1]), pc2_min, pc2_max, plot_bottom, plot_top)
        knot_type = row.get("knotType", "unknown")
        color = color_by_type.get(knot_type, other_color)
        point_title = (
            f"knotType={knot_type}, seed={row.get('seed', '')}, "
            f"PC1={original_coord[0]:.3f}, PC2={original_coord[1]:.3f}"
        )
        svg.append(f'<circle cx="{x:.2f}" cy="{y:.2f}" r="5.4" fill="{color}" fill-opacity="0.84" stroke="#0f172a" stroke-opacity="0.32" stroke-width="0.7"><title>{html.escape(point_title)}</title></circle>')

    legend_x = plot_right + 35
    legend_y = plot_top + 10
    svg.append(f'<text x="{legend_x}" y="{legend_y}" font-size="18" font-weight="700" fill="#0f172a">Knot type</text>')
    legend_y += 28
    for knot_type in top_types:
        color = color_by_type[knot_type]
        label = f"{knot_type} ({knot_counts[knot_type]})"
        svg.append(f'<circle cx="{legend_x + 8}" cy="{legend_y - 5}" r="6" fill="{color}" fill-opacity="0.8"/>')
        svg.append(f'<text x="{legend_x + 24}" y="{legend_y}" font-size="14" fill="#334155">{html.escape(label)}</text>')
        legend_y += 24

    other_count = sum(count for name, count in knot_counts.items() if name not in color_by_type)
    if other_count:
        svg.append(f'<circle cx="{legend_x + 8}" cy="{legend_y - 5}" r="6" fill="{other_color}" fill-opacity="0.8"/>')
        svg.append(f'<text x="{legend_x + 24}" y="{legend_y}" font-size="14" fill="#334155">other ({other_count})</text>')
    svg.append("</svg>")

    document = f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <title>{html.escape(title)}</title>
  <style>
    body {{ margin: 0; background: #e2e8f0; font-family: Georgia, 'Times New Roman', serif; }}
    main {{ max-width: {width}px; margin: 24px auto; box-shadow: 0 20px 60px rgba(15, 23, 42, 0.18); }}
  </style>
</head>
<body>
  <main>{''.join(svg)}</main>
</body>
</html>
"""

    with open(path, "w", encoding="utf-8") as f:
        f.write(document)


def default_output_prefix(input_path):
    path = Path(input_path)
    return str(path.with_suffix(""))


def main():
    parser = argparse.ArgumentParser(description="Run PCA and generate TSV/HTML outputs.")
    parser.add_argument("input", nargs="?", default="Demos/MyCableDemo/knot_seed_catalog_v3_pca_A_density.tsv")
    parser.add_argument("--features", nargs="+", default=DEFAULT_FEATURES)
    parser.add_argument("--pc-count", type=int, default=3)
    parser.add_argument("--output-prefix", default=None)
    parser.add_argument("--title", default="PCA A: density-normalized cable entanglement metrics")
    parser.add_argument(
        "--spread",
        type=float,
        default=0.0,
        help="Visual-only center spreading strength for the HTML plot. Use 0 for raw PCA coordinates.",
    )
    args = parser.parse_args()

    rows = read_tsv(args.input)
    kept_rows, x = load_matrix(rows, args.features)
    mean, std, components, explained, coords = run_pca(x)

    pc_count = min(args.pc_count, coords.shape[1])
    prefix = args.output_prefix or default_output_prefix(args.input)
    coords_path = f"{prefix}_pca_coords.tsv"
    loadings_path = f"{prefix}_pca_loadings.tsv"
    plot_path = f"{prefix}_pca_plot.html"

    write_coords(coords_path, kept_rows, coords, pc_count)
    write_loadings(loadings_path, args.features, components, explained, pc_count)
    write_html_plot(plot_path, kept_rows, coords, args.title, args.spread)

    print(f"Rows loaded: {len(kept_rows)}")
    print("Features:")
    for feature in args.features:
        print(f"  {feature}")
    print("Explained variance:")
    for i in range(pc_count):
        print(f"  PC{i + 1}: {explained[i] * 100.0:.2f}%")
    print(f"Wrote coordinates: {coords_path}")
    print(f"Wrote loadings:    {loadings_path}")
    print(f"Wrote plot:        {plot_path}")


if __name__ == "__main__":
    main()
