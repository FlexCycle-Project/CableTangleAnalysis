#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BIN="${BIN:-$ROOT/bin/MyCableDemo}"
PYTHON_BIN="${PYTHON_BIN:-python3}"
if [ "$PYTHON_BIN" = "python3" ] && [ -x "$ROOT/.venv/bin/python" ]; then
  PYTHON_BIN="$ROOT/.venv/bin/python"
fi

CATALOG=""
CUTS_TO_ZERO=""
STAGE_METRICS=""
OUT_DIR=""
WORKERS=8
START_SEED=""
END_SEED=""
MAX_CUTS=10
MAX_NODES=5000000
CREEP=0.02
SETTLE_STEPS=300
THRESHOLD=0.0001
FORCE=0
IN_PLACE=0

usage() {
  cat <<'EOF'
Usage: scripts/run_oracle_only.sh [options]

Runs the complete hoop-midpoint oracle without regenerating the catalog or
re-running the other cut strategies, then joins the oracle columns into
cuts_to_zero.tsv and appends the oracle stage rows.

Required:
  --catalog PATH            Existing catalog TSV (knotType + seed columns)
  --cuts-to-zero PATH       Existing cuts_to_zero.tsv
  --stage-metrics PATH      Existing cuts_to_zero_stage_metrics.tsv

Options:
  --out-dir DIR             Output directory for parts/logs (default: timestamped)
  --workers N               Parallel seed workers (default: 8)
  --start N                 First seed, inclusive
  --end N                   Last seed, inclusive
  --max-cuts N              Maximum cuts per seed (default: 10)
  --max-nodes N             Search node cap per seed (default: 5000000)
  --creep V                 Rest-creep value (default: 0.02)
  --creep-settle-steps N    Settle steps (default: 300)
  --creep-threshold V       Creep threshold (default: 0.0001)
  --in-place                Replace the input cuts/stage files with merged files
  --force                   Overwrite an existing out-dir
  -h, --help                Show this help
EOF
}

die() {
  echo "error: $*" >&2
  exit 1
}

while [ "$#" -gt 0 ]; do
  case "$1" in
    --catalog) CATALOG="$2"; shift 2 ;;
    --cuts-to-zero) CUTS_TO_ZERO="$2"; shift 2 ;;
    --stage-metrics) STAGE_METRICS="$2"; shift 2 ;;
    --out-dir) OUT_DIR="$2"; shift 2 ;;
    --workers) WORKERS="$2"; shift 2 ;;
    --start) START_SEED="$2"; shift 2 ;;
    --end) END_SEED="$2"; shift 2 ;;
    --max-cuts) MAX_CUTS="$2"; shift 2 ;;
    --max-nodes) MAX_NODES="$2"; shift 2 ;;
    --creep) CREEP="$2"; shift 2 ;;
    --creep-settle-steps) SETTLE_STEPS="$2"; shift 2 ;;
    --creep-threshold) THRESHOLD="$2"; shift 2 ;;
    --force) FORCE=1; shift ;;
    --in-place) IN_PLACE=1; shift ;;
    -h|--help) usage; exit 0 ;;
    *) die "unknown option: $1" ;;
  esac
done

[ -n "$CATALOG" ] || die "--catalog is required"
[ -n "$CUTS_TO_ZERO" ] || die "--cuts-to-zero is required"
[ -n "$STAGE_METRICS" ] || die "--stage-metrics is required"

# Workers cd into per-worker directories, so every path handed to the binary
# must be absolute.  realpath -m does not require the leaf to exist yet.
CATALOG="$(realpath -m "$CATALOG")"
CUTS_TO_ZERO="$(realpath -m "$CUTS_TO_ZERO")"
STAGE_METRICS="$(realpath -m "$STAGE_METRICS")"
BIN="$(realpath -m "$BIN")"

[ -f "$CATALOG" ] || die "catalog not found: $CATALOG"
[ -f "$CUTS_TO_ZERO" ] || die "cuts-to-zero file not found: $CUTS_TO_ZERO"
[ -f "$STAGE_METRICS" ] || die "stage-metrics file not found: $STAGE_METRICS"
[ -x "$BIN" ] || die "binary not executable: $BIN"
[ "$WORKERS" -ge 1 ] || die "--workers must be at least 1"

if [ -z "$OUT_DIR" ]; then
  OUT_DIR="$ROOT/results/oracle_only_$(date +%Y%m%d_%H%M%S)"
fi
OUT_DIR="$(realpath -m "$OUT_DIR")"

if [ -e "$OUT_DIR" ] && [ "$FORCE" -eq 0 ]; then
  die "out-dir already exists (use --force): $OUT_DIR"
fi
if [ "$FORCE" -eq 1 ]; then
  rm -rf "$OUT_DIR"
fi
mkdir -p "$OUT_DIR/parts" "$OUT_DIR/logs" "$OUT_DIR/work"

echo "[oracle-only] splitting catalog"
"$PYTHON_BIN" - "$CATALOG" "$WORKERS" "$START_SEED" "$END_SEED" "$OUT_DIR/parts" <<'PY'
import csv
import sys
from pathlib import Path

catalog = Path(sys.argv[1])
workers = int(sys.argv[2])
start = sys.argv[3]
end = sys.argv[4]
out_dir = Path(sys.argv[5])

with open(catalog, newline="", encoding="utf-8") as handle:
    reader = csv.DictReader(handle, delimiter="\t")
    fieldnames = list(reader.fieldnames or [])
    if "seed" not in fieldnames:
        raise SystemExit("catalog is missing the seed column")
    rows = list(reader)

if start:
    rows = [row for row in rows if int(row["seed"]) >= int(start)]
if end:
    rows = [row for row in rows if int(row["seed"]) <= int(end)]

seen = set()
unique_rows = []
for row in rows:
    seed = row["seed"]
    if seed in seen:
        continue
    seen.add(seed)
    unique_rows.append(row)

buckets = [[] for _ in range(workers)]
for index, row in enumerate(unique_rows):
    buckets[index % workers].append(row)

for index, bucket in enumerate(buckets):
    path = out_dir / f"catalog_{index + 1}.tsv"
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(bucket)
    print(f"  worker {index + 1}: {len(bucket)} seed(s)")

print(f"  total: {len(unique_rows)} seed(s)")
PY

start_time=$(date +%s)
pids=()
worker_numbers=()
for part in "$OUT_DIR"/parts/catalog_*.tsv; do
  worker_name="$(basename "$part" .tsv)"
  worker_name="${worker_name#catalog_}"
  # Skip empty parts (header only).
  if [ "$(wc -l < "$part")" -le 1 ]; then
    continue
  fi
  worker_dir="$OUT_DIR/work/worker_${worker_name}"
  mkdir -p "$worker_dir"
  log="$OUT_DIR/logs/worker_${worker_name}.log"
  (
    cd "$worker_dir"
    xvfb-run -a env OMP_NUM_THREADS=1 "$BIN" \
      --creep="$CREEP" \
      --creep-settle-steps="$SETTLE_STEPS" \
      --creep-threshold="$THRESHOLD" \
      --catalog-file="$part" \
      --eval-ctz-oracle-only \
      --eval-ctz-max-cuts="$MAX_CUTS" \
      --eval-ctz-oracle-max-nodes="$MAX_NODES" \
      --eval-ctz-oracle-output="$OUT_DIR/parts/oracle_${worker_name}.tsv" \
      --eval-ctz-oracle-stage-metrics-output="$OUT_DIR/parts/oracle_stage_${worker_name}.tsv"
  ) > "$log" 2>&1 &
  pids+=("$!")
  worker_numbers+=("$worker_name")
done

status=0
for index in "${!pids[@]}"; do
  pid="${pids[$index]}"
  worker_name="${worker_numbers[$index]}"
  if wait "$pid"; then
    echo "[oracle-only] worker $worker_name finished"
  else
    echo "[oracle-only] worker $worker_name FAILED (see $OUT_DIR/logs/worker_${worker_name}.log)" >&2
    status=1
  fi
done
[ "$status" -eq 0 ] || exit "$status"

end_time=$(date +%s)
elapsed=$((end_time - start_time))

# Concatenate worker outputs, keeping one header.
first=1
: > "$OUT_DIR/oracle_all.tsv"
for part in "$OUT_DIR"/parts/oracle_[0-9]*.tsv; do
  [ -s "$part" ] || continue
  if [ "$first" -eq 1 ]; then
    cat "$part" >> "$OUT_DIR/oracle_all.tsv"
    first=0
  else
    tail -n +2 "$part" >> "$OUT_DIR/oracle_all.tsv"
  fi
done

first=1
: > "$OUT_DIR/oracle_stage_all.tsv"
for part in "$OUT_DIR"/parts/oracle_stage_*.tsv; do
  [ -s "$part" ] || continue
  if [ "$first" -eq 1 ]; then
    cat "$part" >> "$OUT_DIR/oracle_stage_all.tsv"
    first=0
  else
    tail -n +2 "$part" >> "$OUT_DIR/oracle_stage_all.tsv"
  fi
done

[ -s "$OUT_DIR/oracle_all.tsv" ] || die "no oracle rows were produced"
oracle_rows=$(($(wc -l < "$OUT_DIR/oracle_all.tsv") - 1))
[ "$oracle_rows" -ge 1 ] || die "no oracle rows were produced"
echo "[oracle-only] computed $oracle_rows oracle row(s) in ${elapsed}s"

echo "[oracle-only] merging into cuts-to-zero tables"
"$PYTHON_BIN" "$ROOT/scripts/merge_oracle_results.py" \
  --cuts-to-zero "$CUTS_TO_ZERO" \
  --oracle "$OUT_DIR/oracle_all.tsv" \
  --stage-metrics "$STAGE_METRICS" \
  --oracle-stage-metrics "$OUT_DIR/oracle_stage_all.tsv" \
  --output-cuts "$OUT_DIR/cuts_to_zero.tsv" \
  --output-stage-metrics "$OUT_DIR/cuts_to_zero_stage_metrics.tsv"

if [ "$IN_PLACE" -eq 1 ]; then
  cp "$OUT_DIR/cuts_to_zero.tsv" "$CUTS_TO_ZERO"
  cp "$OUT_DIR/cuts_to_zero_stage_metrics.tsv" "$STAGE_METRICS"
  echo "[oracle-only] updated in place:"
  echo "  $CUTS_TO_ZERO"
  echo "  $STAGE_METRICS"
else
  echo "[oracle-only] merged outputs:"
  echo "  $OUT_DIR/cuts_to_zero.tsv"
  echo "  $OUT_DIR/cuts_to_zero_stage_metrics.tsv"
fi
echo "[oracle-only] oracle TSV: $OUT_DIR/oracle_all.tsv"
echo "[oracle-only] elapsed: ${elapsed}s"
