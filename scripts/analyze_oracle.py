#!/usr/bin/env python3
"""Analyze oracle cut performance against threading count and complexity score.

This is the second-pass oracle analysis: it consumes evaluation files that
already contain the PCA-derived complexity score from the first pass.

Outputs per creep variant:
  oracle_summary.tsv
  oracle_by_threading.tsv
  oracle_by_complexity.tsv
  oracle_by_threading_complexity.tsv
  oracle_by_threading.pdf
  oracle_by_complexity.pdf
  oracle_heatmap.pdf
"""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


THREADING_BINS = [
    (1.0, 5.0, "1-5"),
    (6.0, 10.0, "6-10"),
    (11.0, 20.0, "11-20"),
    (21.0, 40.0, "21-40"),
    (41.0, math.inf, "41+"),
]

COMPLEXITY_BINS = ["Q1", "Q2", "Q3", "Q4"]


def read_tsv(path: Path) -> list[dict[str, str]]:
    with open(path, newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str]) -> None:
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


def fmt(value: float, digits: int = 6) -> str:
    return "" if not math.isfinite(value) else f"{value:.{digits}f}"


def threading_bin(t_before: float) -> str:
    for low, high, label in THREADING_BINS:
        if low <= t_before <= high:
            return label
    return THREADING_BINS[-1][2]


def complexity_quartiles(scores: list[float]) -> list[float]:
    finite = [score for score in scores if math.isfinite(score)]
    if not finite:
        return [0.0, 0.0, 0.0, 0.0, 0.0]
    return [float(value) for value in np.quantile(finite, [0.0, 0.25, 0.5, 0.75, 1.0])]


def complexity_bin(score: float, edges: list[float]) -> str | None:
    if not math.isfinite(score):
        return None
    for index in range(4):
        low, high = edges[index], edges[index + 1]
        if index < 3:
            if low <= score < high:
                return COMPLEXITY_BINS[index]
        elif low <= score <= high:
            return COMPLEXITY_BINS[index]
    return None


def finite_values(rows: list[dict[str, Any]], key: str) -> list[float]:
    return [row[key] for row in rows if math.isfinite(row[key])]


def summary_stats(rows: list[dict[str, Any]]) -> dict[str, Any]:
    deltas = finite_values(rows, "delta_oracle")
    fractions = finite_values(rows, "fraction_oracle")
    t_before = finite_values(rows, "t_before")
    heuristic = finite_values(rows, "delta_heuristic")
    beats = [
        row["delta_oracle"] > row["delta_heuristic"]
        for row in rows
        if math.isfinite(row["delta_oracle"]) and math.isfinite(row["delta_heuristic"])
    ]
    return {
        "n": len(rows),
        "mean_delta_oracle": fmt(float(np.mean(deltas))) if deltas else "",
        "median_delta_oracle": fmt(float(np.median(deltas))) if deltas else "",
        "mean_fraction_oracle": fmt(float(np.mean(fractions))) if fractions else "",
        "median_fraction_oracle": fmt(float(np.median(fractions))) if fractions else "",
        "mean_t_before": fmt(float(np.mean(t_before))) if t_before else "",
        "mean_delta_heuristic": fmt(float(np.mean(heuristic))) if heuristic else "",
        "oracle_beats_heuristic_rate": fmt(float(np.mean(beats))) if beats else "",
    }


def load_oracle_rows(path: Path) -> list[dict[str, Any]]:
    rows = read_tsv(path)
    parsed: list[dict[str, Any]] = []
    for row in rows:
        evaluated = str(row.get("oracle_evaluated", "")).strip()
        if evaluated not in {"1", "true", "True"}:
            continue
        t_before = parse_float(row.get("T_before"))
        delta_oracle = parse_float(row.get("delta_oracle"))
        delta_heuristic = parse_float(row.get("delta_heuristic"))
        complexity = parse_float(row.get("complexity_score"))
        if (
            not math.isfinite(t_before)
            or t_before <= 0.0
            or not math.isfinite(delta_oracle)
            or delta_oracle < 0.0
            or not math.isfinite(complexity)
        ):
            continue
        parsed.append(
            {
                "seed": str(row.get("seed", "")),
                "knotType": str(row.get("knotType", "")),
                "t_before": t_before,
                "delta_oracle": delta_oracle,
                "delta_heuristic": delta_heuristic,
                "complexity_score": complexity,
                "fraction_oracle": max(0.0, min(1.0, delta_oracle / t_before)),
            }
        )
    return parsed


def summarize_by_threading(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    result = []
    for _, _, label in THREADING_BINS:
        selected = [row for row in rows if threading_bin(row["t_before"]) == label]
        if not selected:
            continue
        result.append({"threading_bin": label, **summary_stats(selected)})
    return result


def summarize_by_complexity(rows: list[dict[str, Any]], edges: list[float]) -> list[dict[str, Any]]:
    result = []
    for label in COMPLEXITY_BINS:
        selected = [row for row in rows if complexity_bin(row["complexity_score"], edges) == label]
        if not selected:
            continue
        result.append(
            {
                "complexity_bin": label,
                "complexity_min": fmt(edges[COMPLEXITY_BINS.index(label)]),
                "complexity_max": fmt(edges[COMPLEXITY_BINS.index(label) + 1]),
                **summary_stats(selected),
            }
        )
    return result


def summarize_2d(rows: list[dict[str, Any]], edges: list[float]) -> list[dict[str, Any]]:
    result = []
    for _, _, thread_label in THREADING_BINS:
        for complexity_label in COMPLEXITY_BINS:
            selected = [
                row
                for row in rows
                if threading_bin(row["t_before"]) == thread_label
                and complexity_bin(row["complexity_score"], edges) == complexity_label
            ]
            if not selected:
                continue
            result.append(
                {
                    "threading_bin": thread_label,
                    "complexity_bin": complexity_label,
                    **summary_stats(selected),
                }
            )
    return result


def plot_by_threading(rows: list[dict[str, Any]], outpath: Path) -> None:
    labels = [row["threading_bin"] for row in rows]
    x = np.arange(len(labels))
    mean_delta = [parse_float(row["mean_delta_oracle"]) for row in rows]
    mean_fraction = [parse_float(row["mean_fraction_oracle"]) for row in rows]

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2), constrained_layout=True)
    axes[0].plot(x, mean_delta, marker="o", color="#1f77b4", linewidth=2)
    axes[0].set_xticks(x, labels)
    axes[0].set_xlabel("Initial threading count T_before")
    axes[0].set_ylabel("Mean threading removed by oracle")
    axes[0].grid(alpha=0.25)

    axes[1].plot(x, mean_fraction, marker="o", color="#d62728", linewidth=2)
    axes[1].set_ylim(0.0, 1.02)
    axes[1].set_xticks(x, labels)
    axes[1].set_xlabel("Initial threading count T_before")
    axes[1].set_ylabel("Mean fraction removed by oracle")
    axes[1].grid(alpha=0.25)

    fig.suptitle("Oracle single-cut performance by threading number")
    fig.savefig(outpath)
    plt.close(fig)


def plot_by_complexity(rows: list[dict[str, Any]], outpath: Path) -> None:
    labels = [row["complexity_bin"] for row in rows]
    x = np.arange(len(labels))
    mean_delta = [parse_float(row["mean_delta_oracle"]) for row in rows]
    mean_fraction = [parse_float(row["mean_fraction_oracle"]) for row in rows]

    fig, axes = plt.subplots(1, 2, figsize=(9, 4.2), constrained_layout=True)
    axes[0].plot(x, mean_delta, marker="o", color="#2ca02c", linewidth=2)
    axes[0].set_xticks(x, labels)
    axes[0].set_xlabel("Complexity quartile")
    axes[0].set_ylabel("Mean threading removed by oracle")
    axes[0].grid(alpha=0.25)

    axes[1].plot(x, mean_fraction, marker="o", color="#9467bd", linewidth=2)
    axes[1].set_ylim(0.0, 1.02)
    axes[1].set_xticks(x, labels)
    axes[1].set_xlabel("Complexity quartile")
    axes[1].set_ylabel("Mean fraction removed by oracle")
    axes[1].grid(alpha=0.25)

    fig.suptitle("Oracle single-cut performance by complexity score")
    fig.savefig(outpath)
    plt.close(fig)


def plot_heatmap(rows: list[dict[str, Any]], outpath: Path) -> None:
    thread_labels = [label for _, _, label in THREADING_BINS]
    matrix = np.full((len(thread_labels), len(COMPLEXITY_BINS)), np.nan)
    counts = np.zeros_like(matrix)
    for row in rows:
        try:
            i = thread_labels.index(row["threading_bin"])
            j = COMPLEXITY_BINS.index(row["complexity_bin"])
        except ValueError:
            continue
        matrix[i, j] = parse_float(row["mean_fraction_oracle"])
        counts[i, j] = parse_float(row["n"])

    fig, ax = plt.subplots(figsize=(7.5, 5.2), constrained_layout=True)
    image = ax.imshow(matrix, aspect="auto", cmap="viridis", vmin=0.0, vmax=1.0)
    ax.set_xticks(np.arange(len(COMPLEXITY_BINS)), COMPLEXITY_BINS)
    ax.set_yticks(np.arange(len(thread_labels)), thread_labels)
    ax.set_xlabel("Complexity quartile")
    ax.set_ylabel("Initial threading count T_before")
    ax.set_title("Mean fraction of threading removed by oracle")
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            if math.isfinite(matrix[i, j]):
                ax.text(
                    j,
                    i,
                    f"{matrix[i, j]:.2f}\nn={int(counts[i, j])}",
                    ha="center",
                    va="center",
                    color="white" if matrix[i, j] < 0.55 else "black",
                    fontsize=9,
                )
    fig.colorbar(image, ax=ax, label="Mean fraction removed")
    fig.savefig(outpath)
    plt.close(fig)


def analyze_variant(cut_eval_path: Path, out_dir: Path) -> None:
    rows = load_oracle_rows(cut_eval_path)
    if not rows:
        print(f"warning: no oracle rows in {cut_eval_path}; skipping")
        return

    out_dir.mkdir(parents=True, exist_ok=True)
    edges = complexity_quartiles([row["complexity_score"] for row in rows])

    overall = summary_stats(rows)
    threading_rows = summarize_by_threading(rows)
    complexity_rows = summarize_by_complexity(rows, edges)
    combined_rows = summarize_2d(rows, edges)

    summary_fields = [
        "n",
        "mean_delta_oracle",
        "median_delta_oracle",
        "mean_fraction_oracle",
        "median_fraction_oracle",
        "mean_t_before",
        "mean_delta_heuristic",
        "oracle_beats_heuristic_rate",
    ]
    write_tsv(out_dir / "oracle_summary.tsv", [overall], summary_fields)
    write_tsv(
        out_dir / "oracle_by_threading.tsv",
        threading_rows,
        ["threading_bin", *summary_fields],
    )
    write_tsv(
        out_dir / "oracle_by_complexity.tsv",
        complexity_rows,
        ["complexity_bin", "complexity_min", "complexity_max", *summary_fields],
    )
    write_tsv(
        out_dir / "oracle_by_threading_complexity.tsv",
        combined_rows,
        ["threading_bin", "complexity_bin", *summary_fields],
    )

    plot_by_threading(threading_rows, out_dir / "oracle_by_threading.pdf")
    plot_by_complexity(complexity_rows, out_dir / "oracle_by_complexity.pdf")
    plot_heatmap(combined_rows, out_dir / "oracle_heatmap.pdf")

    print(
        f"{cut_eval_path.parent.name}: {len(rows)} oracle samples, "
        f"mean fraction removed={float(overall['mean_fraction_oracle']):.3f}"
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dataset-root", type=Path, required=True)
    parser.add_argument(
        "--complexity-root",
        type=Path,
        default=None,
        help="Per-variant complexity output root (default: <dataset-root>/complexity/per_variant)",
    )
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument(
        "--variants",
        nargs="+",
        default=["zero", "current", "high"],
        help="Creep variants to analyze",
    )
    args = parser.parse_args()

    complexity_root = args.complexity_root or (args.dataset_root / "complexity" / "per_variant")
    for variant in args.variants:
        cut_eval = complexity_root / variant / "cut_eval_with_complexity.tsv"
        if not cut_eval.is_file():
            fallback = args.dataset_root / variant / "cut_eval_with_complexity.tsv"
            cut_eval = fallback if fallback.is_file() else cut_eval
        if not cut_eval.is_file():
            print(f"warning: missing oracle/complexity eval for {variant}: {cut_eval}")
            continue
        analyze_variant(cut_eval, args.out_dir / variant)


if __name__ == "__main__":
    main()
