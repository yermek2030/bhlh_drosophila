#!/usr/bin/env bash
# ============================================================
# STREME background-order sweep on masked sequences.
# Tests whether 3' A-tail persists as background order increases.
# Orders run: 0, 2 (STREME 5.5.5 default), 3, 4.
# Runtime ~1-3 min per run, 16 runs total.
# ============================================================
set -euo pipefail

PROJECT_DIR="${PROJECT_DIR:-$(pwd)}"
cd "$PROJECT_DIR"

OUT="data/motif2_concat_audit/streme_order_sweep"
mkdir -p "$OUT"

NAMES=(Trh Sim Dys Sima)
FASTAS=(
  "data/chapter3/Trh_Tgo_joint_peaks_filt_repeatmasked.fa"
  "data/chapter4/Sim_Tgo_joint_peaks_filt_repeatmasked.fa"
  "data/chapter4/Dys_Tgo_joint_peaks_filt_repeatmasked.fa"
  "data/chapter4/Sima_Tgo_joint_peaks_filt_repeatmasked.fa"
)

for i in "${!NAMES[@]}"; do
  N="${NAMES[$i]}"; FA="${FASTAS[$i]}"
  [[ -f "$FA" ]] || { echo "MISSING: $FA"; exit 1; }
  for ORD in 0 2 3 4; do
    D="$OUT/${N}_order${ORD}"
    echo ">>> STREME $N  --order $ORD"
    conda run -n phd streme \
      --p "$FA" \
      --oc "$D" \
      --dna --minw 6 --maxw 20 --thresh 0.05 \
      --order "$ORD" \
      --verbosity 1
  done
done

echo
echo "=============== MOTIF2 RANK / CONSENSUS BY BACKGROUND ORDER ==============="
{
printf "%-6s %-6s %-5s %-26s %-12s\n" DIMER ORDER RANK CONSENSUS E_VALUE
for N in "${NAMES[@]}"; do
  for ORD in 0 2 3 4; do
    T="$OUT/${N}_order${ORD}/streme.txt"
    [[ -f "$T" ]] || continue
    # Motif2 = the STREME motif whose consensus contains CACTGA (either strand)
    HIT=$(grep "^MOTIF" "$T" | grep -n -E "CACTGA|TCAGTG" | head -1 || true)
    if [[ -z "$HIT" ]]; then
      printf "%-6s %-6s %-5s %-26s %-12s\n" "$N" "$ORD" "-" "NOT_RECOVERED" "-"
      continue
    fi
    RANK="${HIT%%:*}"; REST="${HIT#*:}"
    ID=$(echo "$REST" | awk '{print $2}')
    CONS="${ID#*-}"
    EV=$(awk -v id="$ID" '$1=="MOTIF" && $2==id {f=1} f && /E= /{ \
           match($0,/E= *[0-9.eE+-]+/); print substr($0,RSTART+3,RLENGTH-3); exit}' "$T")
    printf "%-6s %-6s %-5s %-26s %-12s\n" "$N" "$ORD" "$RANK" "$CONS" "${EV:-NA}"
  done
done
} | tee "$OUT/T1_order_sweep_summary.txt"

echo
echo "READ IT LIKE THIS:"
echo "  tail SURVIVES at order 3/4  -> not a composition artefact; keep 15 bp."
echo "  consensus SHORTENS to ~9 bp -> the tail was composition; report the core."
echo "  motif LOST at order 3/4     -> whole motif is composition; serious problem."
