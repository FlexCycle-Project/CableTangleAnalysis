#!/usr/bin/env bash
# Run one cablegen process per seed, with a bounded number of parallel jobs.
#
# Usage:
#   scripts/run_seeds.sh --seeds 1:1000 --jobs 92 --out results/seeds
#   scripts/run_seeds.sh --seeds seeds.txt --jobs 92 --out results/seeds -- --creep=0
#
# The application itself only accepts a single --seed. This driver owns the
# seed iteration and scheduling. Completed seed directories are skipped unless
# --force is given.

set -uo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BIN="${BIN:-$ROOT/bin/cablegen}"

SEEDS=""
JOBS="$(nproc 2>/dev/null || echo 1)"
OUT="$ROOT/results/seeds"
FORCE=0
PASSTHRU=()

usage() {
    cat <<'EOF'
Usage: scripts/run_seeds.sh --seeds <file|start:end> [options] [-- cablegen args...]

Options:
  --seeds SPEC      A file with one seed per line, or a numeric range "start:end".
  --jobs N          Maximum parallel seeds (default: nproc).
  --out DIR         Dataset root (default: <repo>/results/seeds).
  --force           Re-run seeds whose meta.tsv already says "complete".
  -h, --help        Show this help.

Everything after "--" is passed to bin/cablegen for every seed.
The driver sets OMP_NUM_THREADS=1; the intended parallelism unit is one seed
per process.
EOF
}

die() {
    echo "error: $*" >&2
    exit 2
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --seeds)
            [[ $# -ge 2 ]] || die "--seeds needs a value"
            SEEDS="$2"
            shift 2
            ;;
        --seeds=*)
            SEEDS="${1#--seeds=}"
            shift
            ;;
        --jobs)
            [[ $# -ge 2 ]] || die "--jobs needs a value"
            JOBS="$2"
            shift 2
            ;;
        --jobs=*)
            JOBS="${1#--jobs=}"
            shift
            ;;
        --out)
            [[ $# -ge 2 ]] || die "--out needs a value"
            OUT="$2"
            shift 2
            ;;
        --out=*)
            OUT="${1#--out=}"
            shift
            ;;
        --force)
            FORCE=1
            shift
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        --)
            shift
            PASSTHRU=("$@")
            break
            ;;
        *)
            die "unknown argument: $1 (use -- before cablegen arguments)"
            ;;
    esac
done

[[ -n "$SEEDS" ]] || die "--seeds is required"
[[ "$JOBS" =~ ^[0-9]+$ ]] && (( JOBS > 0 )) || die "--jobs must be a positive integer"
[[ -x "$BIN" ]] || die "cablegen not found or not executable: $BIN"

mkdir -p "$OUT"
FAILED="$OUT/failed_seeds.txt"
: > "$FAILED"

seeds=()
if [[ -f "$SEEDS" ]]; then
    while IFS= read -r line; do
        line="${line%%#*}"
        line="${line//[[:space:]]/}"
        [[ -z "$line" ]] && continue
        seeds+=("$line")
    done < "$SEEDS"
else
    if [[ "$SEEDS" =~ ^([0-9]+):([0-9]+)$ ]]; then
        start="${BASH_REMATCH[1]}"
        end="${BASH_REMATCH[2]}"
        if (( start < 1 )); then start=1; fi
        if (( end < start )); then die "empty seed range: $SEEDS"; fi
        for ((s=start; s<=end; ++s)); do
            seeds+=("$s")
        done
    else
        die "--seeds must be a file path or start:end"
    fi
fi

[[ ${#seeds[@]} -gt 0 ]] || die "no seeds selected"

run_one() {
    local seed="$1"
    local dir="$OUT/seed_${seed}"

    # --force replaces the seed directory so no stale optional files (for
    # example an old audit trace) survive from a previous configuration.
    if (( FORCE == 1 )) && [[ -d "$dir" ]]; then
        rm -rf -- "$dir"
    fi

    if (( FORCE == 0 )) && [[ -f "$dir/meta.tsv" ]] && \
       grep -q $'^status\tcomplete$' "$dir/meta.tsv"; then
        echo "[run_seeds] skip seed=$seed (complete)"
        return 0
    fi

    mkdir -p "$dir"
    echo "[run_seeds] start seed=$seed"
    OMP_NUM_THREADS=1 "$BIN" --seed="$seed" --out="$dir" "${PASSTHRU[@]}" \
        > "$dir/run.log" 2>&1
    local rc=$?
    if (( rc != 0 )); then
        echo "[run_seeds] failed seed=$seed rc=$rc"
        printf '%s\n' "$seed" >> "$FAILED"
    else
        echo "[run_seeds] done seed=$seed"
    fi
    return $rc
}

echo "[run_seeds] seeds=${#seeds[@]} jobs=$JOBS out=$OUT"
for seed in "${seeds[@]}"; do
    while (( $(jobs -rp | wc -l) >= JOBS )); do
        wait -n 2>/dev/null || true
    done
    ( run_one "$seed" ) &
done
wait

failed=0
if [[ -s "$FAILED" ]]; then
    failed=$(wc -l < "$FAILED")
fi

{
    printf 'seed\tstatus\n'
    for seed in "${seeds[@]}"; do
        if [[ -f "$OUT/seed_${seed}/meta.tsv" ]] && \
           grep -q $'^status\tcomplete$' "$OUT/seed_${seed}/meta.tsv"; then
            printf '%s\tcomplete\n' "$seed"
        else
            printf '%s\tfailed\n' "$seed"
        fi
    done
} > "$OUT/summary.tsv"

if (( failed > 0 )); then
    echo "[run_seeds] $failed seed(s) failed; see $FAILED" >&2
    exit 1
fi

echo "[run_seeds] all ${#seeds[@]} seed(s) complete; summary: $OUT/summary.tsv"
