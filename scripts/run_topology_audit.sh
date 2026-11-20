#!/usr/bin/env bash
set -euo pipefail

# Runs the MyCableDemo topology audit over a catalog in parallel shards and
# merges the per-shard results into one <prefix>_{seeds,events,summary}.
#
# Each shard is an independent process with its own Xvfb display. The shards
# are disjoint catalog slices, so the merged files are equivalent to a single
# sequential run over the whole catalog.

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BIN="${BIN:-$ROOT/bin/MyCableDemo}"

CATALOG="$ROOT/results/full_dataset_20261106_112203/current/catalog.tsv"
PREFIX=""
WORKERS=32
STEPS=300
CUTS=3
CREEP=0.02
THRESHOLD=0.0001
SETTLE=300
CONTACT_EPS=1e-4
STEPS_PER_FRAME=20
TRACE=0
DT=""
COLLISION_ITERATIONS=""
COLLISION_COMPLIANCE=""
COLLISION_MAX_CORRECTION=""
TANGENTIAL_DAMPING=""
NORMAL_DAMPING=""
HASH_CELL_SIZE=""
START_OFFSET=0
MAX_SEEDS=""

usage() {
  cat <<'EOF'
Usage: scripts/run_topology_audit.sh [options]

  --catalog PATH             Catalog TSV with knotType and seed columns
  --output-prefix PATH       Output prefix (writes <prefix>_seeds.tsv, ...)
  --workers N                Parallel worker processes (default 32)
  --steps N                  Audited substeps per phase (default 300)
  --cuts N                   Heuristic cuts per seed (default 3)
  --creep V                  Rest-state creep rate (default 0.02)
  --threshold V              Creep activation threshold (default 0.0001)
  --settle-steps N           Steps after creep activation before auditing (default 300)
  --contact-eps V            Centre-line contact threshold (default 1e-4)
  --steps-per-frame N        Physics substeps per render/frame tick (default 20)
  --dt V                     Simulation timestep
  --collision-iterations N   Self-collision solver iterations per substep
  --collision-compliance V   Contact compliance
  --collision-max-correction V  Max per-endpoint contact correction
  --tangential-damping V     Contact tangential damping (velocity-proportional)
  --normal-damping V         Contact normal damping
  --hash-cell-size V         Self-collision spatial-hash cell size
  --start N                  Catalog row offset (0-based data rows)
  --max-seeds N              Limit number of catalog seeds audited
  --trace                    Also write per-step trace files
  --help                     Show this message

The merged summary preserves solver/config values from the first shard and uses
max for extrema/penetration fields and sums for counters.
EOF
}

while [ "$#" -gt 0 ]; do
  case "$1" in
    --catalog) CATALOG="$2"; shift 2 ;;
    --output-prefix) PREFIX="$2"; shift 2 ;;
    --workers) WORKERS="$2"; shift 2 ;;
    --steps) STEPS="$2"; shift 2 ;;
    --cuts) CUTS="$2"; shift 2 ;;
    --creep) CREEP="$2"; shift 2 ;;
    --threshold) THRESHOLD="$2"; shift 2 ;;
    --settle-steps) SETTLE="$2"; shift 2 ;;
    --contact-eps) CONTACT_EPS="$2"; shift 2 ;;
    --steps-per-frame) STEPS_PER_FRAME="$2"; shift 2 ;;
    --dt) DT="$2"; shift 2 ;;
    --collision-iterations) COLLISION_ITERATIONS="$2"; shift 2 ;;
    --collision-compliance) COLLISION_COMPLIANCE="$2"; shift 2 ;;
    --collision-max-correction) COLLISION_MAX_CORRECTION="$2"; shift 2 ;;
    --tangential-damping) TANGENTIAL_DAMPING="$2"; shift 2 ;;
    --normal-damping) NORMAL_DAMPING="$2"; shift 2 ;;
    --hash-cell-size) HASH_CELL_SIZE="$2"; shift 2 ;;
    --start) START_OFFSET="$2"; shift 2 ;;
    --max-seeds) MAX_SEEDS="$2"; shift 2 ;;
    --trace) TRACE=1; shift ;;
    --help|-h) usage; exit 0 ;;
    *) echo "unknown option: $1" >&2; usage >&2; exit 2 ;;
  esac
done

if [ ! -f "$CATALOG" ]; then
  echo "catalog not found: $CATALOG" >&2
  exit 1
fi

if [ -z "$PREFIX" ]; then
  PREFIX="$(cd "$(dirname "$CATALOG")" && pwd)/topology_audit"
fi

mkdir -p "$(dirname "$PREFIX")"
WORK="${PREFIX}_work"
rm -rf "$WORK"
mkdir -p "$WORK"

TOTAL_ROWS=$(( $(wc -l < "$CATALOG") - 1 ))
if [ "$TOTAL_ROWS" -le 0 ]; then
  echo "catalog has no data rows: $CATALOG" >&2
  exit 1
fi

SELECTED_CATALOG="$CATALOG"
if [ "$START_OFFSET" != "0" ] || [ -n "$MAX_SEEDS" ]; then
  SELECTED_CATALOG="$WORK/selected.tsv"
  awk -v start="$START_OFFSET" -v max="$MAX_SEEDS" '
    NR == 1 { print; next }
    {
      row = NR - 2
      if (row < start) next
      if (max != "" && (row - start) >= max) exit
      print
    }
  ' "$CATALOG" > "$SELECTED_CATALOG"
  TOTAL_ROWS=$(( $(wc -l < "$SELECTED_CATALOG") - 1 ))
  if [ "$TOTAL_ROWS" -le 0 ]; then
    echo "selected catalog has no data rows" >&2
    exit 1
  fi
fi

if [ "$WORKERS" -gt "$TOTAL_ROWS" ]; then
  WORKERS="$TOTAL_ROWS"
fi

echo "[audit-run] catalog=$CATALOG selectedRows=$TOTAL_ROWS workers=$WORKERS"
echo "[audit-run] stepsPerPhase=$STEPS cuts=$CUTS creep=$CREEP threshold=$THRESHOLD settle=$SETTLE dt=${DT:-default}"

# Split the catalog into contiguous shards, repeating the header in each file.
awk -v workers="$WORKERS" -v work="$WORK" -v total="$TOTAL_ROWS" '
  NR == 1 { header = $0; next }
  {
    row = NR - 2
    part = int(row * workers / total) + 1
    if (part > workers) part = workers
    file = sprintf("%s/part%02d.tsv", work, part)
    if (!(part in seen)) {
      seen[part] = 1
      print header > file
    }
    print > file
  }
' "$SELECTED_CATALOG"

PIDS=()
for ((i = 1; i <= WORKERS; i++)); do
  part=$(printf "%s/part%02d" "$WORK" "$i")
  extra=()
  if [ "$TRACE" -eq 1 ]; then extra+=(--topology-audit-trace); fi
  if [ -n "$DT" ]; then extra+=(--dt="$DT"); fi
  if [ -n "$COLLISION_ITERATIONS" ]; then extra+=(--collision-iterations="$COLLISION_ITERATIONS"); fi
  if [ -n "$COLLISION_COMPLIANCE" ]; then extra+=(--collision-compliance="$COLLISION_COMPLIANCE"); fi
  if [ -n "$COLLISION_MAX_CORRECTION" ]; then extra+=(--collision-max-correction="$COLLISION_MAX_CORRECTION"); fi
  if [ -n "$TANGENTIAL_DAMPING" ]; then extra+=(--tangential-damping="$TANGENTIAL_DAMPING"); fi
  if [ -n "$NORMAL_DAMPING" ]; then extra+=(--normal-damping="$NORMAL_DAMPING"); fi
  if [ -n "$HASH_CELL_SIZE" ]; then extra+=(--hash-cell-size="$HASH_CELL_SIZE"); fi
  if [ -n "$MAX_SEEDS" ] || [ "$START_OFFSET" != "0" ]; then
    # Each shard is already a selected slice; avoid re-limiting inside workers.
    :
  fi

  OMP_NUM_THREADS=1 xvfb-run -a "$BIN" \
    --topology-audit-run \
    --catalog-file="$part.tsv" \
    --topology-audit-steps="$STEPS" \
    --topology-audit-cuts="$CUTS" \
    --topology-audit-contact-eps="$CONTACT_EPS" \
    --creep="$CREEP" \
    --creep-threshold="$THRESHOLD" \
    --creep-settle-steps="$SETTLE" \
    --catalog-steps-per-frame="$STEPS_PER_FRAME" \
    --topology-audit-output="$part" \
    "${extra[@]}" \
    > "$part.log" 2>&1 &

  PIDS+=("$!")
  # Stagger Xvfb startup so auto display selection cannot race.
  sleep 0.4
done

status=0
for pid in "${PIDS[@]}"; do
  if ! wait "$pid"; then
    status=1
  fi
done

if [ "$status" -ne 0 ]; then
  echo "[audit-run] at least one worker failed; see $WORK/part*.log" >&2
  exit 1
fi

# Merge per-shard outputs. The header is taken from the first shard.
merge_tsv() {
  local suffix="$1"
  local target="$2"
  local first=1
  : > "$target"
  for ((i = 1; i <= WORKERS; i++)); do
    local part
    part=$(printf "%s/part%02d_%s" "$WORK" "$i" "$suffix")
    if [ ! -f "$part" ]; then
      echo "[audit-run] missing shard output: $part" >&2
      return 1
    fi
    if [ "$first" -eq 1 ]; then
      head -1 "$part" > "$target"
      first=0
    fi
    tail -n +2 "$part" >> "$target"
  done
  return 0
}

merge_tsv "seeds.tsv" "${PREFIX}_seeds.tsv"
merge_tsv "events.tsv" "${PREFIX}_events.tsv"
if [ "$TRACE" -eq 1 ]; then
  merge_tsv "trace.tsv" "${PREFIX}_trace.tsv"
fi

# Merge shard summaries. Disjoint shards make counters additive; extrema use
# max; solver/config parameters are preserved from the first shard. Values are
# printed in high precision so float parameters and diagnostics survive merge.
python3 - "$WORK" "$PREFIX" <<'PY'
import glob
import os
import sys

work, prefix = sys.argv[1], sys.argv[2]
paths = sorted(glob.glob(os.path.join(work, "part*_summary.txt")))
if not paths:
    raise SystemExit("no shard summaries found")

config_keys = {
    "dt", "creep", "collision_iterations", "collision_compliance",
    "collision_max_correction", "normal_damping", "tangential_damping",
    "hash_cell_size",
}

data = {}
order = []
for path in paths:
    with open(path, "r", encoding="utf-8") as f:
        for raw in f:
            line = raw.strip()
            if not line or "=" not in line:
                continue
            key, val = line.split("=", 1)
            try:
                num = float(val)
            except ValueError:
                continue
            if key not in data:
                data[key] = {"vals": [], "first": num}
                order.append(key)
            data[key]["vals"].append(num)

def is_max_key(key: str) -> bool:
    return key.startswith("max_") or key.endswith("_max") or "maxPenetration" in key or "MaxPenetration" in key or key == "max_applied_correction"

def fmt(num: float) -> str:
    if abs(num - round(num)) < 1e-12 and abs(num) < 2**53:
        return str(int(round(num)))
    return format(num, ".17g")

out_path = f"{prefix}_summary.txt"
with open(out_path, "w", encoding="utf-8") as out:
    for key in order:
        vals = data[key]["vals"]
        first = data[key]["first"]
        if key in config_keys:
            value = first
        elif is_max_key(key):
            value = max(vals) if vals else 0.0
        else:
            value = sum(vals)
        out.write(f"{key}={fmt(value)}\n")
PY

echo "[audit-run] wrote ${PREFIX}_seeds.tsv, ${PREFIX}_events.tsv, ${PREFIX}_summary.txt"
if [ "$TRACE" -eq 1 ]; then
  echo "[audit-run] wrote ${PREFIX}_trace.tsv"
fi
cat "${PREFIX}_summary.txt"
