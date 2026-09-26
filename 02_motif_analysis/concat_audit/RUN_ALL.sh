#!/usr/bin/env bash
# ============================================================
# Motif2 concatenation audit, run in order.
#   T0 info content; T1 STREME sweep ~30-45 min; T2 core scan ~15-25 min; T3 core or/TOMTOM ~10-20 min; T4 verdict; T5 E-value per run
# Outputs: data/motif2_concat_audit/
# Re-run safe: FIMO/STREME cached, R steps overwrite
# ============================================================
set -euo pipefail
cd "$PROJECT_DIR"
D=02_motif_analysis/concat_audit
L=data/motif2_concat_audit/run_log.txt
mkdir -p data/motif2_concat_audit

{
echo "===== $(date) ====="
echo "### T0"; Rscript $D/00_ic_profile.R
echo "### T1"; bash    $D/01_streme_order_sweep.sh
echo "### T2"; Rscript $D/02_split_scan_register.R
echo "### T3"; Rscript $D/03_core_only_OR_and_tomtom.R
echo "### T4"; Rscript $D/04_verdict.R
echo "### T5"; Rscript $D/06_order_sweep_significance.R
} 2>&1 | tee "$L"

echo
echo "Full log: $L"
