#!/usr/bin/env python3
"""Compute a PCA-derived cable complexity score and adapt evaluation outputs.

The PCA is trained on geometric metrics that do not directly depend on the
threading count:

    segmentCount, selfCollisionCount, closePairCount,
    totalCurvature, outsideKnotCurvature

PC1 loadings become the metric weights. The resulting PC1 projection is the
complexity score. Catalog and strategy-evaluation files are augmented with the
score, and strategy performance is summarized by complexity quartile.
"""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path
from typing import Any, Iterable

import numpy as np


DEFAULT_FEATURES = [
    "segmentCount",
    "selfCollisionCount",
    "closePairCount",
    "totalCurvature",
    "outsideKnotCurvature",
]

FEATURE_FALLBACKS = {
    "totalCurvature": ["totalCurvature", "totalCurvatureDensity"],
    "outsideKnotCurvature": ["outsideKnotCurvature", "outsideKnotCurvatureDensity"],
}

UNIDENTIFIED_TYPES = {"", "none", "None", "NONE", "??", "unknown", "Unknown"}
UNKNOT_TYPES = {"UN", "unknot", "Unknot"}

CUT_METRICS = [
    "T_before",
    "delta_heuristic",
    "delta_random_mean",
    "delta_random_hoop_mean",
    "delta_biggest_hoop",
    "delta_smallest_hoop",
    "delta_max_curvature",
    "delta_oracle",
]

CTZ_METRICS = [
    "heuristic_cuts",
    "random_cuts_mean",
    "random_hoop_cuts_mean",
    "biggest_hoop_cuts",
    "smallest_hoop_cuts",
    "max_curvature_cuts",
    "heuristic_resolved",
    "random_resolved_rate",
    "random_hoop_resolved_rate",
    "biggest_hoop_resolved",
    "smallest_hoop_resolved",
    "max_curvature_resolved",
]


def read_tsv(path: Path) -> tuple[list[dict[str, str]], list[str]]:
    with open(path, newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)
        fieldnames = list(reader.fieldnames or [])
    return rows, fieldnames


def write_tsv(path: Path, rows: Iterable[dict[str, Any]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({name: row.get(name, "") for name in fieldnames})


def parse_float(value: Any) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return float("nan")


def fmt(value: float, digits: int = 8) -> str:
    if not math.isfinite(value):
        return ""
    return f"{value:.{digits}f}"


def metric_value(row: dict[str, str], metric: str) -> float:
    if metric == "delta_oracle":
        evaluated = str(row.get("oracle_evaluated", "")).strip()
        if evaluated not in {"1", "true", "True"}:
            return float("nan")
    return parse_float(row.get(metric, ""))


def unique_variant_label(path: Path, used: set[str]) -> str:
    label = path.parent.name or path.stem
    candidate = label
    index = 2
    while candidate in used:
        candidate = f"{label}_{index}"
        index += 1
    used.add(candidate)
    return candidate


def resolve_features(fieldnames: list[str], requested: list[str]) -> dict[str, str]:
    resolved: dict[str, str] = {}
    for feature in requested:
        candidates = FEATURE_FALLBACKS.get(feature, [feature])
        match = next((name for name in candidates if name in fieldnames), None)
        if match is None:
            raise SystemExit(
                f"error: catalog is missing feature '{feature}' "
                f"(tried: {', '.join(candidates)})"
            )
        resolved[feature] = match
    return resolved


def is_identified(types: set[str]) -> bool:
    return any(knot_type not in UNIDENTIFIED_TYPES for knot_type in types)


def is_non_unknot(types: set[str]) -> bool:
    return any(
        knot_type not in UNIDENTIFIED_TYPES and knot_type not in UNKNOT_TYPES
        for knot_type in types
    )


def passes_filter(types: set[str], mode: str) -> bool:
    if mode == "all":
        return True
    if mode == "identified":
        return is_identified(types)
    if mode == "unidentified":
        return not is_identified(types)
    if mode == "non-unknot":
        return is_non_unknot(types)
    raise SystemExit(f"error: unknown knot filter '{mode}'")


def load_groups(
    catalogs: list[Path],
    requested_features: list[str],
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    """Return per-catalog metadata and one deduplicated record per (variant, seed)."""

    used_labels: set[str] = set()
    catalogs_meta: list[dict[str, Any]] = []
    groups: dict[tuple[str, str], dict[str, Any]] = {}

    for catalog in catalogs:
        rows, fieldnames = read_tsv(catalog)
        if "seed" not in fieldnames or "knotType" not in fieldnames:
            raise SystemExit(f"error: {catalog} must contain seed and knotType columns")

        resolved = resolve_features(fieldnames, requested_features)
        variant = unique_variant_label(catalog, used_labels)

        catalogs_meta.append(
            {
                "path": catalog,
                "variant": variant,
                "fieldnames": fieldnames,
                "resolved_features": resolved,
                "rows": rows,
            }
        )

        for row in rows:
            seed = str(row.get("seed", "")).strip()
            if not seed:
                continue
            key = (variant, seed)
            values = [parse_float(row.get(column, "")) for column in resolved.values()]
            if not all(math.isfinite(value) for value in values):
                continue

            knot_type = str(row.get("knotType", "")).strip()
            if key not in groups:
                groups[key] = {
                    "variant": variant,
                    "seed": seed,
                    "knot_types": {knot_type},
                    "values": values,
                }
            else:
                groups[key]["knot_types"].add(knot_type)

    return catalogs_meta, list(groups.values())


def select_training_groups(
    groups: list[dict[str, Any]],
    mode: str,
    min_rows: int,
) -> tuple[list[dict[str, Any]], str]:
    filtered = [group for group in groups if passes_filter(group["knot_types"], mode)]
    if len(filtered) >= min_rows:
        return filtered, mode

    print(
        f"warning: only {len(filtered)} row(s) matched knot filter '{mode}'; "
        f"falling back to all {len(groups)} row(s)."
    )
    return groups, "all"


def fit_pca(
    training_groups: list[dict[str, Any]],
    features: list[str],
) -> dict[str, Any]:
    if len(training_groups) < 2:
        raise SystemExit("error: PCA needs at least two catalog rows")

    x = np.asarray([group["values"] for group in training_groups], dtype=float)
    mean = x.mean(axis=0)
    std = x.std(axis=0)
    std[std == 0.0] = 1.0
    xs = (x - mean) / std

    _, singular_values, vt = np.linalg.svd(xs, full_matrices=False)
    variances = singular_values**2
    total_variance = float(variances.sum())
    explained = variances / total_variance if total_variance > 0.0 else np.zeros_like(variances)

    loadings = vt[0].copy()
    orientation = "sum-positive"
    if loadings.sum() < 0.0:
        loadings *= -1.0
    elif abs(float(loadings.sum())) <= 1e-12:
        first_nonzero = next((i for i, value in enumerate(loadings) if abs(value) > 1e-12), None)
        if first_nonzero is not None and loadings[first_nonzero] < 0.0:
            loadings *= -1.0
        orientation = "first-nonzero-positive"

    return {
        "features": features,
        "mean": mean,
        "std": std,
        "loadings": loadings,
        "explained": explained,
        "orientation": orientation,
    }


def score_groups(groups: list[dict[str, Any]], model: dict[str, Any]) -> list[dict[str, Any]]:
    x = np.asarray([group["values"] for group in groups], dtype=float)
    z = (x - model["mean"]) / model["std"]
    pc1 = z @ model["loadings"]

    minimum = float(pc1.min())
    maximum = float(pc1.max())
    if maximum > minimum:
        normalized = (pc1 - minimum) / (maximum - minimum)
    else:
        normalized = np.zeros_like(pc1)

    order = np.argsort(pc1)
    ranks = np.empty_like(order)
    ranks[order] = np.arange(len(pc1), dtype=float)
    if len(pc1) > 1:
        percentile = ranks / float(len(pc1) - 1)
    else:
        percentile = np.zeros_like(pc1)

    scored: list[dict[str, Any]] = []
    for index, group in enumerate(groups):
        scored.append(
            {
                **group,
                "pc1": float(pc1[index]),
                "score": float(normalized[index]),
                "percentile": float(percentile[index]),
            }
        )
    return scored


def score_map(scored: list[dict[str, Any]]) -> dict[tuple[str, str], dict[str, Any]]:
    return {(item["variant"], item["seed"]): item for item in scored}


def write_model(
    out_dir: Path,
    model: dict[str, Any],
    training_rows: int,
    total_rows: int,
    filter_requested: str,
    filter_used: str,
) -> None:
    rows = []
    for name, mean, std, loading in zip(
        model["features"], model["mean"], model["std"], model["loadings"]
    ):
        rows.append(
            {
                "feature": name,
                "mean": fmt(float(mean)),
                "std": fmt(float(std)),
                "pc1_loading": fmt(float(loading)),
            }
        )
    write_tsv(out_dir / "complexity_model.tsv", rows, ["feature", "mean", "std", "pc1_loading"])

    summary = [
        f"features={','.join(model['features'])}",
        f"knot_filter_requested={filter_requested}",
        f"knot_filter_used={filter_used}",
        f"training_rows={training_rows}",
        f"scored_rows={total_rows}",
        f"pc1_explained_variance={fmt(float(model['explained'][0]))}",
        f"orientation={model['orientation']}",
        "score=pca_pc1(standardized_features)",
    ]
    (out_dir / "complexity_summary.txt").write_text("\n".join(summary) + "\n", encoding="utf-8")


def write_scores(out_dir: Path, scored: list[dict[str, Any]], features: list[str]) -> None:
    fieldnames = [
        "variant",
        "seed",
        "knotType",
        "complexity_pc1",
        "complexity_score",
        "complexity_percentile",
        *features,
    ]
    rows = []
    for item in scored:
        rows.append(
            {
                "variant": item["variant"],
                "seed": item["seed"],
                "knotType": "|".join(sorted(item["knot_types"])),
                "complexity_pc1": fmt(item["pc1"]),
                "complexity_score": fmt(item["score"]),
                "complexity_percentile": fmt(item["percentile"]),
                **{
                    feature: fmt(float(value))
                    for feature, value in zip(features, item["values"])
                },
            }
        )
    write_tsv(out_dir / "complexity_scores.tsv", rows, fieldnames)


def augment_catalog(
    meta: dict[str, Any],
    scores: dict[tuple[str, str], dict[str, Any]],
    output_dir: Path | None = None,
) -> Path:
    catalog = meta["path"]
    variant = meta["variant"]
    destination = output_dir if output_dir is not None else catalog.parent
    output = destination / f"{catalog.stem}_with_complexity{catalog.suffix}"
    fieldnames = list(meta["fieldnames"])
    for name in ("complexity_pc1", "complexity_score", "complexity_percentile"):
        if name not in fieldnames:
            fieldnames.append(name)

    rows = []
    for row in meta["rows"]:
        seed = str(row.get("seed", "")).strip()
        score = scores.get((variant, seed))
        out_row = dict(row)
        out_row["complexity_pc1"] = fmt(score["pc1"]) if score else ""
        out_row["complexity_score"] = fmt(score["score"]) if score else ""
        out_row["complexity_percentile"] = fmt(score["percentile"]) if score else ""
        rows.append(out_row)

    write_tsv(output, rows, fieldnames)
    return output


def augment_eval_file(
    path: Path,
    variant: str,
    scores: dict[tuple[str, str], dict[str, Any]],
    output_dir: Path | None = None,
) -> Path | None:
    if not path.is_file():
        return None

    rows, fieldnames = read_tsv(path)
    if "seed" not in fieldnames:
        return None

    for name in ("complexity_pc1", "complexity_score", "complexity_percentile"):
        if name not in fieldnames:
            fieldnames.append(name)

    output_rows = []
    for row in rows:
        seed = str(row.get("seed", "")).strip()
        score = scores.get((variant, seed))
        out_row = dict(row)
        out_row["complexity_pc1"] = fmt(score["pc1"]) if score else ""
        out_row["complexity_score"] = fmt(score["score"]) if score else ""
        out_row["complexity_percentile"] = fmt(score["percentile"]) if score else ""
        output_rows.append(out_row)

    destination = output_dir if output_dir is not None else path.parent
    output = destination / f"{path.stem}_with_complexity{path.suffix}"
    write_tsv(output, output_rows, fieldnames)
    return output


def quartile_edges(values: list[float]) -> list[float]:
    finite = [value for value in values if math.isfinite(value)]
    if not finite:
        return [0.0, 0.0, 0.0, 0.0, 0.0]
    return [float(value) for value in np.quantile(finite, [0.0, 0.25, 0.5, 0.75, 1.0])]


def summarize_strategy_metrics(
    eval_name: str,
    rows_by_variant: dict[str, list[dict[str, str]]],
    metrics: list[str],
) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    all_rows = [row for rows in rows_by_variant.values() for row in rows]
    scores = [parse_float(row.get("complexity_score", "")) for row in all_rows]
    edges = quartile_edges(scores)

    summary_rows: list[dict[str, str]] = []
    correlation_rows: list[dict[str, str]] = []

    for variant, rows in rows_by_variant.items():
        for metric in metrics:
            if not rows or metric not in rows[0]:
                continue

            pairs = [
                (parse_float(row.get("complexity_score", "")), metric_value(row, metric))
                for row in rows
            ]
            pairs = [(x, y) for x, y in pairs if math.isfinite(x) and math.isfinite(y)]
            if len(pairs) >= 3:
                xs = np.asarray([pair[0] for pair in pairs], dtype=float)
                ys = np.asarray([pair[1] for pair in pairs], dtype=float)
                if xs.std() > 0.0 and ys.std() > 0.0:
                    correlation = float(np.corrcoef(xs, ys)[0, 1])
                else:
                    correlation = float("nan")
                correlation_rows.append(
                    {
                        "eval": eval_name,
                        "variant": variant,
                        "metric": metric,
                        "pearson_r": fmt(correlation),
                        "n": str(len(pairs)),
                    }
                )

        for metric in metrics:
            if not rows or metric not in rows[0]:
                continue
            for bin_index in range(4):
                low = edges[bin_index]
                high = edges[bin_index + 1]
                selected = []
                for row in rows:
                    score = parse_float(row.get("complexity_score", ""))
                    value = metric_value(row, metric)
                    if not math.isfinite(score) or not math.isfinite(value):
                        continue
                    if bin_index < 3:
                        in_bin = low <= score < high
                    else:
                        in_bin = low <= score <= high
                    if in_bin:
                        selected.append(value)
                if not selected:
                    continue
                values = np.asarray(selected, dtype=float)
                summary_rows.append(
                    {
                        "eval": eval_name,
                        "variant": variant,
                        "complexity_bin": str(bin_index + 1),
                        "complexity_min": fmt(low),
                        "complexity_max": fmt(high),
                        "metric": metric,
                        "n": str(len(selected)),
                        "mean": fmt(float(values.mean())),
                        "median": fmt(float(np.median(values))),
                    }
                )

    return summary_rows, correlation_rows


def analyze(
    catalogs: list[Path],
    out_dir: Path,
    features: list[str],
    knot_filter: str,
    min_training_rows: int,
    no_eval_join: bool,
    augment_dir: Path | None = None,
) -> list[dict[str, str]]:
    out_dir.mkdir(parents=True, exist_ok=True)

    catalogs_meta, groups = load_groups(catalogs, features)
    if not groups:
        raise SystemExit("error: no valid catalog rows found")

    training_groups, filter_used = select_training_groups(
        groups, knot_filter, min_training_rows
    )
    model = fit_pca(training_groups, features)
    scored = score_groups(groups, model)
    scores = score_map(scored)

    write_model(
        out_dir,
        model,
        len(training_groups),
        len(groups),
        knot_filter,
        filter_used,
    )
    write_scores(out_dir, scored, features)

    manifest_rows: list[dict[str, str]] = []
    eval_rows: dict[str, dict[str, list[dict[str, str]]]] = {
        "cut_eval": {},
        "cuts_to_zero": {},
    }

    for meta in catalogs_meta:
        variant = meta["variant"]
        variant_dir = meta["path"].parent
        augmented_catalog = augment_catalog(meta, scores, augment_dir)

        cut_joined = None
        ctz_joined = None
        if not no_eval_join:
            cut_joined = augment_eval_file(
                variant_dir / "cut_eval.tsv", variant, scores, augment_dir
            )
            ctz_joined = augment_eval_file(
                variant_dir / "cuts_to_zero.tsv", variant, scores, augment_dir
            )
            if cut_joined is not None:
                eval_rows["cut_eval"][variant] = read_tsv(cut_joined)[0]
            if ctz_joined is not None:
                eval_rows["cuts_to_zero"][variant] = read_tsv(ctz_joined)[0]

        manifest_rows.append(
            {
                "variant": variant,
                "catalog": str(meta["path"]),
                "catalog_with_complexity": str(augmented_catalog),
                "cut_eval_with_complexity": str(cut_joined) if cut_joined else "",
                "cuts_to_zero_with_complexity": str(ctz_joined) if ctz_joined else "",
            }
        )

    summary_rows, correlation_rows = summarize_strategy_metrics(
        "cut_eval", eval_rows["cut_eval"], CUT_METRICS
    )
    ctz_summary_rows, ctz_correlation_rows = summarize_strategy_metrics(
        "cuts_to_zero", eval_rows["cuts_to_zero"], CTZ_METRICS
    )
    summary_rows.extend(ctz_summary_rows)
    correlation_rows.extend(ctz_correlation_rows)

    write_tsv(
        out_dir / "complexity_manifest.tsv",
        manifest_rows,
        [
            "variant",
            "catalog",
            "catalog_with_complexity",
            "cut_eval_with_complexity",
            "cuts_to_zero_with_complexity",
        ],
    )
    write_tsv(
        out_dir / "strategy_by_complexity.tsv",
        summary_rows,
        [
            "eval",
            "variant",
            "complexity_bin",
            "complexity_min",
            "complexity_max",
            "metric",
            "n",
            "mean",
            "median",
        ],
    )
    write_tsv(
        out_dir / "strategy_correlations.tsv",
        correlation_rows,
        ["eval", "variant", "metric", "pearson_r", "n"],
    )

    print(f"PCA rows: {len(training_groups)} of {len(groups)}")
    print(f"PC1 explained variance: {model['explained'][0] * 100.0:.2f}%")
    print("PC1 weights:")
    for name, loading in zip(features, model["loadings"]):
        print(f"  {name}: {loading:.6f}")
    print(f"Wrote complexity outputs to: {out_dir}")

    return manifest_rows


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("catalogs", nargs="+", type=Path, help="Catalog TSV file(s)")
    parser.add_argument("--out-dir", type=Path, required=True, help="Output directory")
    parser.add_argument(
        "--features",
        nargs="+",
        default=DEFAULT_FEATURES,
        help="PCA feature columns (default: non-threading geometry metrics)",
    )
    parser.add_argument(
        "--knot-filter",
        choices=["all", "identified", "unidentified", "non-unknot"],
        default="identified",
        help="Rows used to fit PCA (default: identified knots)",
    )
    parser.add_argument(
        "--min-training-rows",
        type=int,
        default=3,
        help="Fall back to all rows if the filter yields fewer rows (default: 3)",
    )
    parser.add_argument(
        "--no-eval-join",
        action="store_true",
        help="Only augment the catalog; do not join complexity into evaluation TSVs",
    )
    parser.add_argument(
        "--per-variant",
        action="store_true",
        help="Fit a separate PCA model for each catalog and write one subdirectory per variant",
    )
    args = parser.parse_args()

    if args.per_variant:
        used_labels: set[str] = set()
        manifest_rows: list[dict[str, str]] = []
        for catalog in args.catalogs:
            label = unique_variant_label(catalog, used_labels)
            variant_out = args.out_dir / label
            print(f"\n=== PCA for variant '{label}' ===")
            rows = analyze(
                [catalog],
                variant_out,
                args.features,
                args.knot_filter,
                args.min_training_rows,
                args.no_eval_join,
                augment_dir=variant_out,
            )
            for row in rows:
                manifest_rows.append(
                    {
                        "variant": row["variant"],
                        "out_dir": str(variant_out),
                        "model": str(variant_out / "complexity_model.tsv"),
                        "scores": str(variant_out / "complexity_scores.tsv"),
                        "catalog_with_complexity": row["catalog_with_complexity"],
                        "cut_eval_with_complexity": row["cut_eval_with_complexity"],
                        "cuts_to_zero_with_complexity": row["cuts_to_zero_with_complexity"],
                    }
                )

        write_tsv(
            args.out_dir / "per_variant_manifest.tsv",
            manifest_rows,
            [
                "variant",
                "out_dir",
                "model",
                "scores",
                "catalog_with_complexity",
                "cut_eval_with_complexity",
                "cuts_to_zero_with_complexity",
            ],
        )
    else:
        analyze(
            args.catalogs,
            args.out_dir,
            args.features,
            args.knot_filter,
            args.min_training_rows,
            args.no_eval_join,
        )


if __name__ == "__main__":
    main()
