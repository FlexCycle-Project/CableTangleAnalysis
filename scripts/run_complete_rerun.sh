#!/usr/bin/env bash
# Self-contained orchestrator for the complete pipeline rerun. Runs stages
# sequentially, logging each to its own file under logs/. Exits nonzero on
# stage failure. Designed to be launched fully detached (see launch_cmd).
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
OUT="$ROOT/results/rerun_20261124"
LOGS="$OUT/stage_logs"
mkdir -p "$LOGS"

log() { echo "[$(date +%H:%M:%S)] $*" | tee -a "$LOGS/orchestrator.log"; }

# ---------------------------------------------------------------- stage 1
log "STAGE 1: cut evaluation + cuts-to-zero (60 workers)"
"$ROOT/scripts/resume_rerun.sh" 2>&1 | tee "$LOGS/stage1_resume.log"
log "STAGE 1 done"

# ---------------------------------------------------------------- stage 2
log "STAGE 2: Julia complexity + figures already produced by resume; verifying outputs"
[ -f "$OUT/catalog_with_complexity.tsv" ] && log "  catalog_with_complexity present"
[ -d "$OUT/analysis_julia/single" ] && log "  analysis_julia/single present"
[ -d "$OUT/analysis_oracle" ] && log "  analysis_oracle present"

# ---------------------------------------------------------------- stage 3
log "STAGE 3: exhaustive cuts-to-zero oracle (run_oracle_only.sh, 60 workers)"
"$ROOT/scripts/run_oracle_only.sh" \
  --catalog "$OUT/catalog.tsv" \
  --cuts-to-zero "$OUT/cuts_to_zero.tsv" \
  --stage-metrics "$OUT/cuts_to_zero_stage_metrics.tsv" \
  --out-dir "$OUT/oracle_only_run" \
  --workers 60 \
  --creep 0.02 --creep-settle-steps 300 --creep-threshold 0.0001 \
  --in-place 2>&1 | tee "$LOGS/stage3_oracle.log"
log "STAGE 3 done"

# ---------------------------------------------------------------- stage 4
log "STAGE 4: topology audit (60 workers, 300 steps, 3 cuts)"
"$ROOT/scripts/run_topology_audit.sh" \
  --catalog "$OUT/catalog.tsv" \
  --output-prefix "$OUT/topology_audit" \
  --workers 60 \
  --steps 300 --cuts 3 \
  --creep 0.02 --threshold 0.0001 --settle-steps 300 2>&1 | tee "$LOGS/stage4_audit.log"
log "STAGE 4 done"

# ---------------------------------------------------------------- stage 5
log "STAGE 5: all pipeline stages complete"
echo "DATASET ROOT: $OUT"
echo "DONE"
