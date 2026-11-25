#!/usr/bin/env bash
# Run cut_eval (+ oracle) and cuts-to-zero for a CableTangleAnalysis catalog,
# then the Julia analysis and oracle-vs-threading analysis. Uses the same
# seed-offset approach as run_full_pipeline.sh (--catalog-start/--catalog-max
# per worker, no file splitting). Resumable: re-running skips workers whose
# cut AND ctz outputs already exist.
#
# Usage: run_cut_phases.sh [workers] [results_dir]
#   workers     parallel workers (default 24)
#   results_dir directory with catalog.tsv (default results/rerun_20261124)
set -euo pipefail
W="${1:-42}"
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BIN="$ROOT/bin/MyCableDemo"
OUT="${2:-$ROOT/results/rerun_20261124}"
CREEP=0.02; SETTLE=300; THRESHOLD=0.0001
EVAL_RANDOM_TRIALS=10; CTZ_MAX_CUTS=10; CTZ_RANDOM_TRIALS=20
JULIA_PROJECT="$ROOT/cabletangleanalysis"
export KYMOKNOT_ROOT="${KYMOKNOT_ROOT:-$ROOT/KymoKnot-clone}"
export OMP_NUM_THREADS=1

[ -x "$BIN" ] || { echo "no binary at $BIN (build it first: cmake --build build)"; exit 1; }
[ -f "$OUT/catalog.tsv" ] || { echo "no catalog at $OUT/catalog.tsv"; exit 1; }

mkdir -p "$OUT/parts" "$OUT/logs"
TOTAL=$(( $(wc -l < "$OUT/catalog.tsv") - 1 ))   # data rows (seeds)
echo "[run] workers=$W total_seeds=$TOTAL out=$OUT $(date +%H:%M:%S)"

# Assign contiguous SEED ranges to workers, exactly as run_full_pipeline.sh:
# seeds run 1..TOTAL; worker w takes `count` seeds starting at `first_seed`.
# The catalog is sorted by ascending seed, so seed == file row index here.
BASE=$((TOTAL / W)); REM=$((TOTAL % W))

run_worker() {
  set -euo pipefail
  local w="$1"
  local count=$BASE; [ "$w" -le "$REM" ] && count=$((count+1))
  [ "$count" -ge 1 ] || { echo "[w$w] no seeds"; return 0; }
  local first_seed=$(( (w-1)*BASE + (w<=REM ? w-1 : REM) + 1 ))
  local cutp="$OUT/parts/run_cut_$(printf %02d $w).tsv"
  local cutsp="$OUT/parts/run_cut_stage_$(printf %02d $w).tsv"
  local ctp="$OUT/parts/run_ctz_$(printf %02d $w).tsv"
  local ctzsp="$OUT/parts/run_ctz_stage_$(printf %02d $w).tsv"
  [ -s "$cutp" ] && [ -s "$ctp" ] && { echo "[w$w] already complete ($count seeds), skip"; return 0; }
  local wdir="$OUT/work/run_${w}"; mkdir -p "$wdir"; cd "$wdir"
  {
    echo "[w$w] first_seed=$first_seed count=$count cut_eval begin $(date +%H:%M:%S)"
    [ -s "$cutp" ] || xvfb-run -a "$BIN" --creep="$CREEP" --creep-settle-steps="$SETTLE" --creep-threshold="$THRESHOLD" \
      --catalog-file="$OUT/catalog.tsv" --catalog-start="$first_seed" --catalog-max="$count" \
      --eval-run --eval-oracle \
      --eval-output="$cutp" --eval-stage-metrics-output="$cutsp" --eval-random-trials="$EVAL_RANDOM_TRIALS" \
      || echo "[w$w] cut_eval FAILED"
    echo "[w$w] ctz begin $(date +%H:%M:%S)"
    [ -s "$ctp" ] || xvfb-run -a "$BIN" --creep="$CREEP" --creep-settle-steps="$SETTLE" --creep-threshold="$THRESHOLD" \
      --catalog-file="$OUT/catalog.tsv" --catalog-start="$first_seed" --catalog-max="$count" \
      --eval-cuts-to-zero \
      --eval-ctz-output="$ctp" --eval-ctz-stage-metrics-output="$ctzsp" \
      --eval-ctz-max-cuts="$CTZ_MAX_CUTS" --eval-ctz-random-trials="$CTZ_RANDOM_TRIALS" \
      || echo "[w$w] ctz FAILED"
    echo "[w$w] done $(date +%H:%M:%S)"
  } > "$OUT/logs/run_${w}.log" 2>&1
}

PIDS=()
for w in $(seq 1 "$W"); do run_worker "$w" & PIDS+=("$!"); done
FAIL=0
for pid in "${PIDS[@]}"; do wait "$pid" || FAIL=1; done
[ "$FAIL" -eq 0 ] || { echo "SOME WORKERS FAILED - see $OUT/logs/run_*.log"; exit 1; }

merge_parts() {
  local output="$1"; shift; local first=1
  : > "$output"
  for part in "$@"; do
    [ -s "$part" ] || { echo "missing $part"; exit 1; }
    if [ "$first" -eq 1 ]; then cat "$part" > "$output"; first=0; else tail -n +2 "$part" >> "$output"; fi
  done
}
merge_parts "$OUT/cut_eval.tsv" $(for w in $(seq 1 $W); do echo "$OUT/parts/run_cut_$(printf %02d $w).tsv"; done)
merge_parts "$OUT/cuts_to_zero.tsv" $(for w in $(seq 1 $W); do echo "$OUT/parts/run_ctz_$(printf %02d $w).tsv"; done)
merge_parts "$OUT/cut_eval_stage_metrics.tsv" $(for w in $(seq 1 $W); do echo "$OUT/parts/run_cut_stage_$(printf %02d $w).tsv"; done)
merge_parts "$OUT/cuts_to_zero_stage_metrics.tsv" $(for w in $(seq 1 $W); do echo "$OUT/parts/run_ctz_stage_$(printf %02d $w).tsv"; done)

echo "=== merged row counts ==="
for f in cut_eval.tsv cuts_to_zero.tsv cut_eval_stage_metrics.tsv cuts_to_zero_stage_metrics.tsv; do
  echo "$f: $(($(wc -l < "$OUT/$f") - 1))"
done

echo "=== Julia analysis ==="
export CABLESIM_PCA_ORIENTATION=positive CABLESIM_COMPLEXITY_FILTER=identified CABLESIM_COMPLEXITY_PER_VARIANT=1
julia --project="$JULIA_PROJECT" "$ROOT/cabletangleanalysis/scripts/run_dataset.jl" "$OUT" single="$OUT"
echo "=== oracle-vs-threading analysis ==="
python3 "$ROOT/scripts/analyze_oracle.py" --dataset-root "$OUT" --out-dir "$OUT/analysis_oracle"
echo "DONE"
