#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BIN="${BIN:-$ROOT/bin/MyCableDemo}"
KYMOKNOT_ROOT="${KYMOKNOT_ROOT:-$ROOT/KymoKnot-clone}"
KYMOKNOT_BIN="${KYMOKNOT_BIN:-$KYMOKNOT_ROOT/KymoKnot_linear.x}"
PYTHON_BIN="${PYTHON_BIN:-python3}"
if [ "$PYTHON_BIN" = "python3" ] && [ -x "$ROOT/.venv/bin/python" ]; then
  PYTHON_BIN="$ROOT/.venv/bin/python"
fi

START_SEED=1
END_SEED=1000
WORKERS=8
OUT_DIR=""
FORCE=0
EVAL_RANDOM_TRIALS=10
EVAL_ORACLE=1
CTZ_MAX_CUTS=10
CTZ_RANDOM_TRIALS=20
RUN_COMPLEXITY=1
COMPLEXITY_FILTER=identified
PCA_ORIENTATION=positive
COMPLEXITY_PER_VARIANT=1
RUN_ORACLE_ANALYSIS=1
SINGLE_CREEP=""
CREEP_CURRENT=0.02
CREEP_HIGH=0.2
CREEP_SETTLE_STEPS=300
CREEP_THRESHOLD=0.0001

usage() {
  cat <<'EOF'
Usage: scripts/run_full_pipeline.sh [options]

Runs catalog generation, cut evaluation, and cuts-to-zero evaluation for a
seed range. By default it runs three rest-state-creep variants:

  zero      creep = 0
  current   creep = 0.02 (the code default)
  high      creep = 0.2

Each variant is written to its own subdirectory. Use --creep to run only one
creep value instead.

Options:
  --start N                 First seed (default: 1)
  --end N                   Last seed, inclusive (default: 1000)
  --workers N               Number of parallel seed workers (default: 8)
  --out DIR                 Output root (default: results/full_dataset_<timestamp>)
  --force                   Overwrite existing generated files
  --oracle                  Enable the exhaustive single-cut oracle (default)
  --no-oracle               Skip the exhaustive single-cut oracle
  --no-complexity           Skip the PCA-derived complexity post-processing
  --no-complexity-per-variant
                            Skip the separate-PCA-per-creep pass (default: on)
  --no-oracle-analysis      Skip the oracle vs threading/complexity analysis
  --complexity-filter MODE  PCA training rows: all, identified, unidentified,
                            or non-unknot (default: identified)
  --pca-orientation MODE    positive: absolute PC1 weights (default);
                            signed: keep signed PC1 projection
  --creep V                 Run only this creep value (single-variant mode)
  --creep-current V         "current" creep value (default: 0.02)
  --creep-high V            "high" creep value (default: 0.2)
  --creep-settle-steps N    Steps after creep activation before evaluation (default: 300)
  --creep-threshold V       Creep activation threshold (default: 0.0001)
  --eval-random-trials N    Cut-eval random edge/hoop trials per seed (default: 10)
  --ctz-max-cuts N          Max cuts per seed in cuts-to-zero (default: 10)
  --ctz-random-trials N     Cuts-to-zero random edge/hoop trials per seed (default: 20)
  -h, --help                Show this help

Notes:
  * Seed 0 is normalized to 1 by MyCableDemo. The script clamps --start 0 to 1.
  * Use --start 1 --end 1000 for 1000 unique seeds.
  * The script sets OMP_NUM_THREADS=1 per worker. This is much faster than
    letting every process spawn one OpenMP thread per CPU core.
  * Xvfb is used because MyCableDemo still creates an OpenGL window when run
    in batch mode.
  * Running all three creep variants takes roughly three times as long as a
    single variant.
  * The script lowers the creep activation threshold to 0.0001. With the code
    default of 0.01, the mismatch under gravity-only settling stays below the
    threshold and all creep variants produce identical samples.
  * The complexity pass requires Python 3 with NumPy and Matplotlib.
    Use --no-complexity or --no-oracle-analysis to skip parts of it.
  * Per-variant PCA and the oracle analysis run in multi-variant mode only.
EOF
}

die() {
  echo "error: $*" >&2
  exit 1
}

while [ "$#" -gt 0 ]; do
  case "$1" in
    --start)
      [ "$#" -ge 2 ] || die "--start requires a value"
      START_SEED="$2"
      shift 2
      ;;
    --end)
      [ "$#" -ge 2 ] || die "--end requires a value"
      END_SEED="$2"
      shift 2
      ;;
    --workers)
      [ "$#" -ge 2 ] || die "--workers requires a value"
      WORKERS="$2"
      shift 2
      ;;
    --out)
      [ "$#" -ge 2 ] || die "--out requires a value"
      OUT_DIR="$2"
      shift 2
      ;;
    --force)
      FORCE=1
      shift
      ;;
    --oracle)
      EVAL_ORACLE=1
      shift
      ;;
    --no-oracle)
      EVAL_ORACLE=0
      shift
      ;;
    --no-complexity)
      RUN_COMPLEXITY=0
      shift
      ;;
    --no-complexity-per-variant)
      COMPLEXITY_PER_VARIANT=0
      shift
      ;;
    --no-oracle-analysis)
      RUN_ORACLE_ANALYSIS=0
      shift
      ;;
    --complexity-filter)
      [ "$#" -ge 2 ] || die "--complexity-filter requires a value"
      COMPLEXITY_FILTER="$2"
      shift 2
      ;;
    --pca-orientation)
      [ "$#" -ge 2 ] || die "--pca-orientation requires a value"
      PCA_ORIENTATION="$2"
      shift 2
      ;;
    --creep)
      [ "$#" -ge 2 ] || die "--creep requires a value"
      SINGLE_CREEP="$2"
      shift 2
      ;;
    --creep-current)
      [ "$#" -ge 2 ] || die "--creep-current requires a value"
      CREEP_CURRENT="$2"
      shift 2
      ;;
    --creep-high)
      [ "$#" -ge 2 ] || die "--creep-high requires a value"
      CREEP_HIGH="$2"
      shift 2
      ;;
    --creep-settle-steps)
      [ "$#" -ge 2 ] || die "--creep-settle-steps requires a value"
      CREEP_SETTLE_STEPS="$2"
      shift 2
      ;;
    --creep-threshold)
      [ "$#" -ge 2 ] || die "--creep-threshold requires a value"
      CREEP_THRESHOLD="$2"
      shift 2
      ;;
    --eval-random-trials)
      [ "$#" -ge 2 ] || die "--eval-random-trials requires a value"
      EVAL_RANDOM_TRIALS="$2"
      shift 2
      ;;
    --ctz-max-cuts)
      [ "$#" -ge 2 ] || die "--ctz-max-cuts requires a value"
      CTZ_MAX_CUTS="$2"
      shift 2
      ;;
    --ctz-random-trials)
      [ "$#" -ge 2 ] || die "--ctz-random-trials requires a value"
      CTZ_RANDOM_TRIALS="$2"
      shift 2
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      die "unknown option: $1 (use --help)"
      ;;
  esac
done

for value_name in START_SEED END_SEED WORKERS CREEP_SETTLE_STEPS EVAL_RANDOM_TRIALS CTZ_MAX_CUTS CTZ_RANDOM_TRIALS; do
  value="${!value_name}"
  [[ "$value" =~ ^[0-9]+$ ]] || die "$value_name must be a non-negative integer (got '$value')"
done

is_creep_number() {
  [[ "$1" =~ ^[0-9]*\.?[0-9]+$ ]]
}

if [ -n "$SINGLE_CREEP" ]; then
  is_creep_number "$SINGLE_CREEP" || die "--creep must be a non-negative number (got '$SINGLE_CREEP')"
else
  is_creep_number "$CREEP_CURRENT" || die "--creep-current must be a non-negative number (got '$CREEP_CURRENT')"
  is_creep_number "$CREEP_HIGH" || die "--creep-high must be a non-negative number (got '$CREEP_HIGH')"
fi
is_creep_number "$CREEP_THRESHOLD" || die "--creep-threshold must be a non-negative number (got '$CREEP_THRESHOLD')"

if [ "$START_SEED" -eq 0 ]; then
  echo "note: seed 0 is normalized to 1 by MyCableDemo; starting at 1."
  START_SEED=1
fi

[ "$END_SEED" -ge "$START_SEED" ] || die "--end must be >= --start"
[ "$WORKERS" -ge 1 ] || die "--workers must be at least 1"
case "$COMPLEXITY_FILTER" in
  all|identified|unidentified|non-unknot) ;;
  *) die "invalid --complexity-filter '$COMPLEXITY_FILTER'" ;;
esac
case "$PCA_ORIENTATION" in
  positive|signed) ;;
  *) die "invalid --pca-orientation '$PCA_ORIENTATION'" ;;
esac
[ "$EVAL_RANDOM_TRIALS" -ge 1 ] || die "--eval-random-trials must be at least 1"
[ "$CTZ_MAX_CUTS" -ge 1 ] || die "--ctz-max-cuts must be at least 1"
[ "$CTZ_RANDOM_TRIALS" -ge 1 ] || die "--ctz-random-trials must be at least 1"

TOTAL=$((END_SEED - START_SEED + 1))
if [ "$WORKERS" -gt "$TOTAL" ]; then
  WORKERS="$TOTAL"
fi

[ -x "$BIN" ] || die "MyCableDemo is not executable: $BIN (build it first)"
[ -x "$KYMOKNOT_BIN" ] || die "KymoKnot is not executable: $KYMOKNOT_BIN (run: make -C KymoKnot-clone)"
command -v Xvfb >/dev/null 2>&1 || die "Xvfb not found; install it (e.g. apt install xvfb)"

if [ -n "$SINGLE_CREEP" ]; then
  MODE="single"
  VARIANTS=("single:$SINGLE_CREEP")
  DEFAULT_OUT="$ROOT/results/catalog_$(date +%Y%m%d_%H%M%S)"
else
  MODE="multi"
  VARIANTS=("zero:0" "current:$CREEP_CURRENT" "high:$CREEP_HIGH")
  DEFAULT_OUT="$ROOT/results/full_dataset_$(date +%Y%m%d_%H%M%S)"
fi

if [ -z "$OUT_DIR" ]; then
  OUT_DIR="$DEFAULT_OUT"
fi
mkdir -p "$OUT_DIR"
OUT_DIR="$(cd "$OUT_DIR" && pwd)"
mkdir -p "$OUT_DIR/logs"

MANIFEST="$OUT_DIR/manifest.tsv"
if [ "$FORCE" -eq 1 ]; then
  rm -f "$MANIFEST"
else
  [ ! -e "$MANIFEST" ] || \
    die "$MANIFEST already exists; use --force or a different --out"
fi

export KYMOKNOT_ROOT
export OMP_NUM_THREADS=1

DISPLAY_NUM=90
while [ -e "/tmp/.X11-unix/X$DISPLAY_NUM" ] || pgrep -f "Xvfb :$DISPLAY_NUM " >/dev/null 2>&1; do
  DISPLAY_NUM=$((DISPLAY_NUM + 1))
  [ "$DISPLAY_NUM" -lt 250 ] || die "could not find a free X display number"
done

XVFB_LOG="$OUT_DIR/xvfb.log"
Xvfb ":$DISPLAY_NUM" -screen 0 1280x1024x24 > "$XVFB_LOG" 2>&1 &
XVFB_PID=$!

cleanup() {
  kill "$XVFB_PID" 2>/dev/null || true
}
trap cleanup EXIT

for _ in $(seq 1 50); do
  [ -e "/tmp/.X11-unix/X$DISPLAY_NUM" ] && break
  sleep 0.1
done
[ -e "/tmp/.X11-unix/X$DISPLAY_NUM" ] || {
  cat "$XVFB_LOG" >&2 || true
  die "Xvfb failed to start on display :$DISPLAY_NUM"
}
export DISPLAY=":$DISPLAY_NUM"

run_worker() {
  set -euo pipefail

  local worker="$1"
  local first_seed="$2"
  local last_seed="$3"
  local variant_out="$4"
  local creep="$5"
  local count=$((last_seed - first_seed + 1))

  local work_dir="$variant_out/work/worker_${worker}"
  local log="$variant_out/logs/worker_${worker}.log"
  local catalog_part="$variant_out/parts/catalog_${worker}.tsv"
  local cut_part="$variant_out/parts/cut_${worker}.tsv"
  local ctz_part="$variant_out/parts/ctz_${worker}.tsv"
  local seed_bank="$variant_out/parts/seed_bank_${worker}.txt"
  local oracle_args=()
  if [ "$EVAL_ORACLE" -eq 1 ]; then
    oracle_args+=(--eval-oracle)
  fi

  mkdir -p "$work_dir"
  cd "$work_dir"

  {
    echo "[worker $worker] creep=$creep seeds ${first_seed}-${last_seed} (${count})"
    echo "[worker $worker] catalog generation"
    "$BIN" \
      --creep="$creep" \
      --creep-settle-steps="$CREEP_SETTLE_STEPS" \
      --creep-threshold="$CREEP_THRESHOLD" \
      --catalog-run \
      --catalog-start="$first_seed" \
      --catalog-max="$count" \
      --catalog-file="$catalog_part" \
      --seed-bank="$seed_bank"

    echo "[worker $worker] cut evaluation"
    "$BIN" \
      --creep="$creep" \
      --creep-settle-steps="$CREEP_SETTLE_STEPS" \
      --creep-threshold="$CREEP_THRESHOLD" \
      --catalog-file="$catalog_part" \
      --eval-run \
      --eval-output="$cut_part" \
      --eval-random-trials="$EVAL_RANDOM_TRIALS" \
      "${oracle_args[@]}"

    echo "[worker $worker] cuts-to-zero evaluation"
    "$BIN" \
      --creep="$creep" \
      --creep-settle-steps="$CREEP_SETTLE_STEPS" \
      --creep-threshold="$CREEP_THRESHOLD" \
      --catalog-file="$catalog_part" \
      --eval-cuts-to-zero \
      --eval-ctz-output="$ctz_part" \
      --eval-ctz-max-cuts="$CTZ_MAX_CUTS" \
      --eval-ctz-random-trials="$CTZ_RANDOM_TRIALS"

    echo "[worker $worker] done"
  } > "$log" 2>&1
}

merge_parts() {
  local output_file="$1"
  shift
  local first=1
  local part
  : > "$output_file"
  for part in "$@"; do
    [ -s "$part" ] || {
      echo "error: missing or empty part: $part" >&2
      return 1
    }
    if [ "$first" -eq 1 ]; then
      cat "$part" >> "$output_file"
      first=0
    else
      tail -n +2 "$part" >> "$output_file"
    fi
  done
}

run_variant() {
  local label="$1"
  local creep="$2"
  local variant_out="$3"

  mkdir -p "$variant_out"
  variant_out="$(cd "$variant_out" && pwd)"

  local start_ts
  start_ts=$(date +%s)

  echo "variant: label=$label creep=$creep output=$variant_out"

  if [ "$FORCE" -eq 1 ]; then
    rm -rf "$variant_out/parts" "$variant_out/logs" "$variant_out/work"
    rm -f "$variant_out/catalog.tsv" "$variant_out/cut_eval.tsv" \
      "$variant_out/cuts_to_zero.tsv" "$variant_out/summary.txt"
  else
    for output_name in catalog.tsv cut_eval.tsv cuts_to_zero.tsv; do
      [ ! -e "$variant_out/$output_name" ] || \
        die "$variant_out/$output_name already exists; use --force or a different --out"
    done
  fi

  mkdir -p "$variant_out/parts" "$variant_out/logs" "$variant_out/work"

  local base=$((TOTAL / WORKERS))
  local remainder=$((TOTAL % WORKERS))
  local next_seed="$START_SEED"
  local pids=()

  for worker in $(seq 1 "$WORKERS"); do
    local count="$base"
    if [ "$worker" -le "$remainder" ]; then
      count=$((count + 1))
    fi
    local first_seed="$next_seed"
    local last_seed=$((next_seed + count - 1))
    next_seed=$((last_seed + 1))

    run_worker "$worker" "$first_seed" "$last_seed" "$variant_out" "$creep" &
    pids+=("$!")
  done

  local failed=0
  for pid in "${pids[@]}"; do
    if ! wait "$pid"; then
      failed=1
    fi
  done

  if [ "$failed" -ne 0 ]; then
    echo "error: at least one worker failed; see $variant_out/logs" >&2
    return 1
  fi

  local catalog_parts=()
  local cut_parts=()
  local ctz_parts=()
  for worker in $(seq 1 "$WORKERS"); do
    catalog_parts+=("$variant_out/parts/catalog_${worker}.tsv")
    cut_parts+=("$variant_out/parts/cut_${worker}.tsv")
    ctz_parts+=("$variant_out/parts/ctz_${worker}.tsv")
  done

  merge_parts "$variant_out/catalog.tsv" "${catalog_parts[@]}"
  merge_parts "$variant_out/cut_eval.tsv" "${cut_parts[@]}"
  merge_parts "$variant_out/cuts_to_zero.tsv" "${ctz_parts[@]}"

  local unique_seeds
  unique_seeds="$(awk -F'\t' 'NR > 1 && $2 ~ /^[0-9]+$/ { print $2 }' "$variant_out/catalog.tsv" | sort -n | uniq | wc -l)"
  local catalog_rows
  catalog_rows=$(($(wc -l < "$variant_out/catalog.tsv") - 1))
  local cut_rows
  cut_rows=$(($(wc -l < "$variant_out/cut_eval.tsv") - 1))
  local ctz_rows
  ctz_rows=$(($(wc -l < "$variant_out/cuts_to_zero.tsv") - 1))

  [ "$unique_seeds" -eq "$TOTAL" ] || \
    die "catalog has $unique_seeds unique seeds, expected $TOTAL"
  [ "$cut_rows" -eq "$TOTAL" ] || \
    die "cut evaluation has $cut_rows rows, expected $TOTAL"
  [ "$ctz_rows" -eq "$TOTAL" ] || \
    die "cuts-to-zero evaluation has $ctz_rows rows, expected $TOTAL"

  local end_ts
  end_ts=$(date +%s)
  local elapsed=$((end_ts - start_ts))

  {
    echo "variant=$label"
    echo "creep=$creep"
    echo "start_seed=$START_SEED"
    echo "end_seed=$END_SEED"
    echo "seeds=$TOTAL"
    echo "workers=$WORKERS"
    echo "elapsed_seconds=$elapsed"
    echo "catalog_rows=$catalog_rows"
    echo "unique_seeds=$unique_seeds"
    echo "cut_rows=$cut_rows"
    echo "ctz_rows=$ctz_rows"
    echo "creep_settle_steps=$CREEP_SETTLE_STEPS"
    echo "creep_threshold=$CREEP_THRESHOLD"
    echo "eval_random_trials=$EVAL_RANDOM_TRIALS"
    echo "eval_oracle=$EVAL_ORACLE"
    echo "ctz_max_cuts=$CTZ_MAX_CUTS"
    echo "ctz_random_trials=$CTZ_RANDOM_TRIALS"
  } | tee "$variant_out/summary.txt"

  echo "catalog:       $variant_out/catalog.tsv"
  echo "cut eval:      $variant_out/cut_eval.tsv"
  echo "cuts to zero:  $variant_out/cuts_to_zero.tsv"
  echo "logs:          $variant_out/logs"
}

START_ALL_TS=$(date +%s)
printf 'variant\tcreep\tsettle_steps\tthreshold\tseeds\telapsed_seconds\tcatalog\tcut_eval\tcuts_to_zero\n' > "$MANIFEST"
variant_catalogs=()

for variant_spec in "${VARIANTS[@]}"; do
  label="${variant_spec%%:*}"
  creep="${variant_spec#*:}"

  if [ "$MODE" = "single" ]; then
    variant_out="$OUT_DIR"
  else
    variant_out="$OUT_DIR/$label"
  fi

  run_variant "$label" "$creep" "$variant_out"
  variant_catalogs+=("$variant_out/catalog.tsv")

  elapsed="$(awk -F= '$1 == "elapsed_seconds" { print $2 }' "$variant_out/summary.txt")"
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$label" "$creep" "$CREEP_SETTLE_STEPS" "$CREEP_THRESHOLD" "$TOTAL" "$elapsed" \
    "$variant_out/catalog.tsv" "$variant_out/cut_eval.tsv" \
    "$variant_out/cuts_to_zero.tsv" >> "$MANIFEST"
done

END_ALL_TS=$(date +%s)
TOTAL_ELAPSED=$((END_ALL_TS - START_ALL_TS))

if [ "$RUN_COMPLEXITY" -eq 1 ] && [ "$TOTAL" -lt 2 ]; then
  echo "warning: skipping complexity pass; at least two seeds are required for PCA"
elif [ "$RUN_COMPLEXITY" -eq 1 ]; then
  { [ -x "$PYTHON_BIN" ] || command -v "$PYTHON_BIN" >/dev/null 2>&1; } || \
    die "python3 not found; use --no-complexity to skip the PCA pass"
  "$PYTHON_BIN" -c 'import numpy' >/dev/null 2>&1 || \
    die "NumPy not found; install requirements.txt or use --no-complexity"

  echo "running pooled PCA complexity pass"
  "$PYTHON_BIN" "$ROOT/scripts/compute_complexity.py" \
    "${variant_catalogs[@]}" \
    --out-dir "$OUT_DIR/complexity" \
    --knot-filter "$COMPLEXITY_FILTER" \
    --pca-orientation "$PCA_ORIENTATION"

  if [ "$MODE" = "multi" ] && [ "$COMPLEXITY_PER_VARIANT" -eq 1 ]; then
    echo "running per-variant PCA complexity pass"
    "$PYTHON_BIN" "$ROOT/scripts/compute_complexity.py" \
      "${variant_catalogs[@]}" \
      --out-dir "$OUT_DIR/complexity/per_variant" \
      --knot-filter "$COMPLEXITY_FILTER" \
      --pca-orientation "$PCA_ORIENTATION" \
      --per-variant

    echo "per-variant:   $OUT_DIR/complexity/per_variant"
  fi
fi

if [ "$MODE" = "multi" ] && [ "$RUN_ORACLE_ANALYSIS" -eq 1 ] && [ "$EVAL_ORACLE" -eq 1 ] && [ "$TOTAL" -ge 2 ]; then
  { [ -x "$PYTHON_BIN" ] || command -v "$PYTHON_BIN" >/dev/null 2>&1; } || \
    die "python3 not found; use --no-oracle-analysis"
  "$PYTHON_BIN" -c 'import matplotlib' >/dev/null 2>&1 || \
    die "matplotlib not found; install requirements.txt or use --no-oracle-analysis"

  echo "running oracle threading/complexity analysis"
  "$PYTHON_BIN" "$ROOT/scripts/analyze_oracle.py" \
    --dataset-root "$OUT_DIR" \
    --complexity-root "$OUT_DIR/complexity/per_variant" \
    --out-dir "$OUT_DIR/analysis_oracle"

  echo "oracle analysis: $OUT_DIR/analysis_oracle"
fi

echo "dataset root:  $OUT_DIR"
echo "manifest:      $MANIFEST"
echo "total elapsed: ${TOTAL_ELAPSED}s"
