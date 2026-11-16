#!/usr/bin/env python3
"""Merge complete-oracle results into the cuts-to-zero evaluation tables.

The C++ oracle-only pass writes a compact oracle TSV and its own stage-metrics
TSV.  This script joins the oracle columns into cuts_to_zero.tsv by seed and
appends the oracle stage rows to cuts_to_zero_stage_metrics.tsv.

The merge is idempotent: existing oracle_* columns and existing oracle stage
rows are replaced.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

ORACLE_COLUMNS = [
    "oracle_cuts",
    "oracle_resolved",
    "oracle_exhaustive",
    "oracle_nodes",
    "oracle_sequence",
]


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with open(path, newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fieldnames = list(reader.fieldnames or [])
        rows = list(reader)
    if not fieldnames:
        raise SystemExit(f"missing header in {path}")
    return fieldnames, rows


def write_tsv(path: Path, fieldnames: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="ignore",
        )
        writer.writeheader()
        for row in rows:
            writer.writerow({name: row.get(name, "") for name in fieldnames})


def merge_cuts(
    cuts_path: Path,
    oracle_path: Path,
    output_path: Path,
) -> None:
    cuts_fields, cuts_rows = read_tsv(cuts_path)
    oracle_fields, oracle_rows = read_tsv(oracle_path)

    missing = [name for name in ("seed", *ORACLE_COLUMNS) if name not in oracle_fields]
    if missing:
        raise SystemExit(f"oracle TSV is missing columns: {', '.join(missing)}")

    oracle_by_seed = {row["seed"]: row for row in oracle_rows}

    for column in ORACLE_COLUMNS:
        if column not in cuts_fields:
            cuts_fields.append(column)

    for row in cuts_rows:
        oracle = oracle_by_seed.get(row.get("seed", ""))
        for column in ORACLE_COLUMNS:
            row[column] = oracle.get(column, "") if oracle is not None else ""

    write_tsv(output_path, cuts_fields, cuts_rows)


def merge_stage_metrics(
    stage_path: Path,
    oracle_stage_path: Path,
    output_path: Path,
) -> None:
    stage_fields, stage_rows = read_tsv(stage_path)
    oracle_fields, oracle_rows = read_tsv(oracle_stage_path)

    if stage_fields != oracle_fields:
        raise SystemExit(
            "stage-metrics headers differ between the base file and the oracle file"
        )

    kept = [row for row in stage_rows if row.get("strategy") != "oracle"]
    write_tsv(output_path, stage_fields, kept + oracle_rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cuts-to-zero", type=Path, required=True)
    parser.add_argument("--oracle", type=Path, required=True)
    parser.add_argument("--stage-metrics", type=Path, required=True)
    parser.add_argument("--oracle-stage-metrics", type=Path, required=True)
    parser.add_argument("--output-cuts", type=Path, required=True)
    parser.add_argument("--output-stage-metrics", type=Path, required=True)
    args = parser.parse_args()

    merge_cuts(args.cuts_to_zero, args.oracle, args.output_cuts)
    merge_stage_metrics(
        args.stage_metrics,
        args.oracle_stage_metrics,
        args.output_stage_metrics,
    )


if __name__ == "__main__":
    main()
