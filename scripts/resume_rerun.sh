#!/usr/bin/env bash
# Resume the full pipeline for results/rerun_20261124 from the completed
# catalog: split the catalog into per-worker parts, run cut evaluation
# (with the exhaustive single-cut oracle) and cuts-to-zero evaluation in
# parallel, merge the parts, then run the Julia complexity pass and the
# oracle-vs-threading analysis. Mirrors scripts/run_full_pipeline.sh --creep 0.02
# with the canonical settings from results/dataset_final/summary.txt.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BIN="$ROOT/bin/MyCableDemo"
OUT="$ROOT/results/rerun_20261124"
WORKERS=60
CREEP=0.02
SETTLE=300
THRESHOLD=0.0001
EVAL_RANDOM_TRIALS=10
CTZ_MAX_CUTS=10
CTZ_RANDOM_TRIALS=20
JULIA_PROJECT="$ROOT/cabletangleanalysis"

[ -x "$BIN" ] || { echo "no binary"; exit 1; }
[ -f "$OUT/catalog.tsv" ] || { echo "no catalog at $OUT/catalog.tsv"; exit 1; }

mkdir -p "$OUT/parts" "$OUT/logs"
rm -f "$OUT/parts/cut_"*.tsv "$OUT/parts/ctz_"*.tsv "$OUT/parts/cut_stage_"*.tsv "$OUT/parts/ctz_stage_"*.tsv
rm -f "$OUT/logs/resume_"*.log

# Original unpadded catalog parts catalog_${worker}.tsv already exist from the
# full-pipeline catalog phase and are complete (1000 seeds), so no re-split is
# needed. Verify cover the full seed range before proceeding.
TOTAL=$(( $(wc -l < "$OUT/catalog.tsv") - 1 ))
FOUND=0
for w in $(seq 1 "$WORKERS"); do [ -s "$OUT/parts/catalog_${w}.tsv" ] && FOUND=$((FOUND+1)); done
[ "$FOUND" -eq "$WORKERS" ] || { echo "expected $WORKERS catalog parts, found $FOUND"; exit 1; }

export KYMOKNOT_ROOT="${KYMOKNOT_ROOT:-$ROOT/KymoKnot-clone}"
export OMP_NUM_THREADS=1

DISPLAY_NUM=90
while [ -e "/tmp/.X11-unix/X$DISPLAY_NUM" ]; do DISPLAY_NUM=$((DISPLAY_NUM+1)); done
Xvfb ":$DISPLAY_NUM" -screen 0 1280x1024x24 > "$OUT/xvfb.log" 2>&1 &
XVFB_PID=$!
export DISPLAY=":$DISPLAY_NUM"
for _ in $(seq 1 50); do [ -e "/tmp/.X11-unix/X$DISPLAY_NUM" ] && break; sleep 0.1; done

run_worker() {
  set -euo pipefail
  local worker="$1"; local part="$2"
  local wdir="$OUT/work/resume_${worker}"
  mkdir -p "$wdir"; cd "$wdir"
  local catp="$OUT/parts/catalog_${worker}.tsv"
  local cutp="$OUT/parts/cut_${worker}.tsv"
  local ctp="$OUT/parts/ctz_${worker}.tsv"
  local cutsp="$OUT/parts/cut_stage_${worker}.tsv"
  local ctzsp="$OUT/parts/ctz_stage_${worker}.tsv"
  {
    echo "[worker $worker] cut evaluation"
    "$BIN" --creep="$CREEP" --creep-settle-steps="$SETTLE" --creep-threshold="$THRESHOLD" \
      --catalog-file="$catp" --eval-run --eval-oracle \
      --eval-output="$cutp" --eval-stage-metrics-output="$cutsp" \
      --eval-random-trials="$EVAL_RANDOM_TRIALS"
    echo "[worker $worker] cuts-to-zero evaluation"
    "$BIN" --creep="$CREEP" --creep-settle-steps="$SETTLE" --creep-threshold="$THRESHOLD" \
      --catalog-file="$catp" --eval-cuts-to-zero \
      --eval-ctz-output="$ctp" --eval-ctz-stage-metrics-output="$ctzsp" \
      --eval-ctz-max-cuts="$CTZ_MAX_CUTS" --eval-ctz-random-trials="$CTZ_RANDOM_TRIALS"
    echo "[worker $worker] done"
  } > "$OUT/logs/resume_${worker}.log" 2>&1
}

PIDS=()
for w in $(seq 1 "$WORKERS"); do
  run_worker "$w" "$OUT/parts/catalog_${w}.tsv" &
  PIDS+=("$!")
done

FAIL=0
for pid in "${PIDS[@]}"; do
  wait "$pid" || FAIL=1
done
kill "$XVFB_PID" 2>/dev/null || true
[ "$FAIL" -eq 0 ] || { echo "RESUME FAILED - see $OUT/logs/resume_*.log"; exit 1; }

merge_parts() {
  local output="$1"; shift
  local first=1
  : > "$output"
  for part in "$@"; do
    [ -s "$part" ] || { echo "missing $part"; exit 1; }
    if [ "$first" -eq 1 ]; then cat "$part" > "$output"; first=0; else tail -n +2 "$part" >> "$output"; fi
  done
}
merge_parts "$OUT/cut_eval.tsv" $(for w in $(seq 1 $WORKERS); do echo "$OUT/parts/cut_${w}.tsv"; done)
merge_parts "$OUT/cuts_to_zero.tsv" $(for w in $(seq 1 $WORKERS); do echo "$OUT/parts/ctz_${w}.tsv"; done)
merge_parts "$OUT/cut_eval_stage_metrics.tsv" $(for w in $(seq 1 $WORKERS); do echo "$OUT/parts/cut_stage_${w}.tsv"; done)
merge_parts "$OUT/cuts_to_zero_stage_metrics.tsv" $(for w in $(seq 1 $WORKERS); do echo "$OUT/parts/ctz_stage_${w}.tsv"; done)

echo "=== merged row counts ==="
for f in cut_eval.tsv cuts_to_zero.tsv cut_eval_stage_metrics.tsv cuts_to_zero_stage_metrics.tsv; do
  echo "$f: $(($(wc -l < "$OUT/$f") - 1))"
done

# Julia complexity + figures
echo "=== running Julia analysis ==="
export CABLESIM_PCA_ORIENTATION=positive
export CABLESIM_COMPLEXITY_FILTER=identified
export CABLESIM_COMPLEXITY_PER_VARIANT=1
julia --project="$JULIA_PROJECT" "$ROOT/cabletangleanalysis/scripts/run_dataset.jl" "$OUT" single="$OUT"
echo "=== oracle-vs-threading analysis ==="
python3 "$ROOT/scripts/analyze_oracle.py" --dataset-root "$OUT" --out-dir "$OUT/analysis_oracle"
echo "DONE"
