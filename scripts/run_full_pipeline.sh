#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BIN="${BIN:-$ROOT/bin/MyCableDemo}"
KYMOKNOT_ROOT="${KYMOKNOT_ROOT:-$ROOT/KymoKnot-clone}"
KYMOKNOT_BIN="${KYMOKNOT_BIN:-$KYMOKNOT_ROOT/KymoKnot_linear.x}"

START_SEED=1
END_SEED=1000
WORKERS=8
OUT_DIR=""
FORCE=0
EVAL_RANDOM_TRIALS=10
CTZ_MAX_CUTS=10
CTZ_RANDOM_TRIALS=20

usage() {
  cat <<'EOF'
Usage: scripts/run_full_pipeline.sh [options]

Runs catalog generation, cut-heuristic evaluation, and cuts-to-zero
evaluation for a seed range, with the range split across parallel workers.

Options:
  --start N                 First seed (default: 1)
  --end N                   Last seed, inclusive (default: 1000)
  --workers N               Number of parallel workers (default: 8)
  --out DIR                 Output directory (default: results/catalog_<timestamp>)
  --force                   Overwrite output files in an existing --out directory
  --eval-random-trials N    Cut-eval random trials per seed (default: 10)
  --ctz-max-cuts N          Max cuts per seed in cuts-to-zero (default: 10)
  --ctz-random-trials N     Cuts-to-zero random trials per seed (default: 20)
  -h, --help                Show this help

Notes:
  * Seed 0 is normalized to 1 by MyCableDemo. The script clamps --start 0 to 1.
  * Use --start 1 --end 1000 for 1000 unique seeds.
  * The script sets OMP_NUM_THREADS=1 per worker. This is much faster than
    letting every process spawn one OpenMP thread per CPU core.
  * Xvfb is used because MyCableDemo still creates an OpenGL window when run
    in batch mode.
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

for value_name in START_SEED END_SEED WORKERS EVAL_RANDOM_TRIALS CTZ_MAX_CUTS CTZ_RANDOM_TRIALS; do
  value="${!value_name}"
  [[ "$value" =~ ^[0-9]+$ ]] || die "$value_name must be a non-negative integer (got '$value')"
done

if [ "$START_SEED" -eq 0 ]; then
  echo "note: seed 0 is normalized to 1 by MyCableDemo; starting at 1."
  START_SEED=1
fi

[ "$END_SEED" -ge "$START_SEED" ] || die "--end must be >= --start"
[ "$WORKERS" -ge 1 ] || die "--workers must be at least 1"
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

if [ -z "$OUT_DIR" ]; then
  OUT_DIR="$ROOT/results/catalog_$(date +%Y%m%d_%H%M%S)"
fi
mkdir -p "$OUT_DIR"
OUT_DIR="$(cd "$OUT_DIR" && pwd)"

if [ "$FORCE" -eq 1 ]; then
  rm -rf "$OUT_DIR/parts" "$OUT_DIR/logs" "$OUT_DIR/work"
  rm -f "$OUT_DIR/catalog.tsv" "$OUT_DIR/cut_eval.tsv" "$OUT_DIR/cuts_to_zero.tsv"
else
  for output_name in catalog.tsv cut_eval.tsv cuts_to_zero.tsv; do
    [ ! -e "$OUT_DIR/$output_name" ] || \
      die "$OUT_DIR/$output_name already exists; use --force or a different --out"
  done
fi

mkdir -p "$OUT_DIR/parts" "$OUT_DIR/logs" "$OUT_DIR/work"

export KYMOKNOT_ROOT
export OMP_NUM_THREADS=1

DISPLAY_NUM=90
while [ -e "/tmp/.X11-unix/X$DISPLAY_NUM" ] || pgrep -f "Xvfb :$DISPLAY_NUM " >/dev/null 2>&1; do
  DISPLAY_NUM=$((DISPLAY_NUM + 1))
  [ "$DISPLAY_NUM" -lt 250 ] || die "could not find a free X display number"
done

Xvfb ":$DISPLAY_NUM" -screen 0 1280x1024x24 > "$OUT_DIR/logs/xvfb.log" 2>&1 &
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
  cat "$OUT_DIR/logs/xvfb.log" >&2 || true
  die "Xvfb failed to start on display :$DISPLAY_NUM"
}
export DISPLAY=":$DISPLAY_NUM"

run_worker() {
  set -euo pipefail

  local worker="$1"
  local first_seed="$2"
  local last_seed="$3"
  local count=$((last_seed - first_seed + 1))

  local work_dir="$OUT_DIR/work/worker_${worker}"
  local log="$OUT_DIR/logs/worker_${worker}.log"
  local catalog_part="$OUT_DIR/parts/catalog_${worker}.tsv"
  local cut_part="$OUT_DIR/parts/cut_${worker}.tsv"
  local ctz_part="$OUT_DIR/parts/ctz_${worker}.tsv"
  local seed_bank="$OUT_DIR/parts/seed_bank_${worker}.txt"

  mkdir -p "$work_dir"
  cd "$work_dir"

  {
    echo "[worker $worker] seeds ${first_seed}-${last_seed} (${count})"
    echo "[worker $worker] catalog generation"
    "$BIN" \
      --catalog-run \
      --catalog-start="$first_seed" \
      --catalog-max="$count" \
      --catalog-file="$catalog_part" \
      --seed-bank="$seed_bank"

    echo "[worker $worker] cut-heuristic evaluation"
    "$BIN" \
      --catalog-file="$catalog_part" \
      --eval-run \
      --eval-output="$cut_part" \
      --eval-random-trials="$EVAL_RANDOM_TRIALS"

    echo "[worker $worker] cuts-to-zero evaluation"
    "$BIN" \
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

START_TS=$(date +%s)
echo "run: seeds ${START_SEED}-${END_SEED} (${TOTAL}), workers=${WORKERS}"
echo "output: $OUT_DIR"

base=$((TOTAL / WORKERS))
remainder=$((TOTAL % WORKERS))
next_seed="$START_SEED"
pids=()

for worker in $(seq 1 "$WORKERS"); do
  count="$base"
  if [ "$worker" -le "$remainder" ]; then
    count=$((count + 1))
  fi
  first_seed="$next_seed"
  last_seed=$((next_seed + count - 1))
  next_seed=$((last_seed + 1))

  run_worker "$worker" "$first_seed" "$last_seed" &
  pids+=("$!")
done

failed=0
for pid in "${pids[@]}"; do
  if ! wait "$pid"; then
    failed=1
  fi
done

if [ "$failed" -ne 0 ]; then
  echo "error: at least one worker failed; see $OUT_DIR/logs" >&2
  exit 1
fi

catalog_parts=()
cut_parts=()
ctz_parts=()
for worker in $(seq 1 "$WORKERS"); do
  catalog_parts+=("$OUT_DIR/parts/catalog_${worker}.tsv")
  cut_parts+=("$OUT_DIR/parts/cut_${worker}.tsv")
  ctz_parts+=("$OUT_DIR/parts/ctz_${worker}.tsv")
done

merge_parts "$OUT_DIR/catalog.tsv" "${catalog_parts[@]}"
merge_parts "$OUT_DIR/cut_eval.tsv" "${cut_parts[@]}"
merge_parts "$OUT_DIR/cuts_to_zero.tsv" "${ctz_parts[@]}"

unique_seeds="$(awk -F'\t' 'NR > 1 && $2 ~ /^[0-9]+$/ { print $2 }' "$OUT_DIR/catalog.tsv" | sort -n | uniq | wc -l)"
catalog_rows=$(($(wc -l < "$OUT_DIR/catalog.tsv") - 1))
cut_rows=$(($(wc -l < "$OUT_DIR/cut_eval.tsv") - 1))
ctz_rows=$(($(wc -l < "$OUT_DIR/cuts_to_zero.tsv") - 1))

[ "$unique_seeds" -eq "$TOTAL" ] || \
  die "catalog has $unique_seeds unique seeds, expected $TOTAL"
[ "$cut_rows" -eq "$TOTAL" ] || \
  die "cut evaluation has $cut_rows rows, expected $TOTAL"
[ "$ctz_rows" -eq "$TOTAL" ] || \
  die "cuts-to-zero evaluation has $ctz_rows rows, expected $TOTAL"

END_TS=$(date +%s)
ELAPSED=$((END_TS - START_TS))

{
  echo "start_seed=$START_SEED"
  echo "end_seed=$END_SEED"
  echo "seeds=$TOTAL"
  echo "workers=$WORKERS"
  echo "elapsed_seconds=$ELAPSED"
  echo "catalog_rows=$catalog_rows"
  echo "unique_seeds=$unique_seeds"
  echo "cut_rows=$cut_rows"
  echo "ctz_rows=$ctz_rows"
  echo "eval_random_trials=$EVAL_RANDOM_TRIALS"
  echo "ctz_max_cuts=$CTZ_MAX_CUTS"
  echo "ctz_random_trials=$CTZ_RANDOM_TRIALS"
} | tee "$OUT_DIR/summary.txt"

echo "catalog:       $OUT_DIR/catalog.tsv"
echo "cut eval:      $OUT_DIR/cut_eval.tsv"
echo "cuts to zero:  $OUT_DIR/cuts_to_zero.tsv"
echo "logs:          $OUT_DIR/logs"
