#!/usr/bin/env bash
# Run cut_eval (+ oracle) for a slice of catalog workers, writing each part
# incrementally so partial progress survives. Args: <worker-start> <worker-end>
# Each invocation covers worker parts [start, end] (inclusive, 1-based).
# Uses the completed catalog parts in results/rerun_20261124/parts/catalog_N.tsv.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BIN="$ROOT/bin/MyCableDemo"
OUT="$ROOT/results/rerun_20261124"
START=${1:?worker start}
END=${2:?worker end}
export KYMOKNOT_ROOT="${KYMOKNOT_ROOT:-$ROOT/KymoKnot-clone}"
export OMP_NUM_THREADS=1

for w in $(seq "$START" "$END"); do
  catp="$OUT/parts/catalog_${w}.tsv"
  [ -s "$catp" ] || { echo "skip w$w: no catalog part"; continue; }
  cutp="$OUT/parts/cut_${w}.tsv"
  cutsp="$OUT/parts/cut_stage_${w}.tsv"
  [ -s "$cutp" ] && [ -s "$cutsp" ] && { echo "w$w cut_eval already present, skip"; continue; }
  echo "[chunk] w$w cut_eval $(date +%H:%M:%S)"
  xvfb-run -a "$BIN" \
    --creep=0.02 --creep-settle-steps=300 --creep-threshold=0.0001 \
    --catalog-file="$catp" --eval-run --eval-oracle \
    --eval-output="$cutp" --eval-stage-metrics-output="$cutsp" \
    --eval-random-trials=10 > "$OUT/logs/chunk_$w.log" 2>&1 \
    || echo "[chunk] w$w FAILED (see $OUT/logs/chunk_$w.log)"
done
echo "[chunk] done $START..$END"
