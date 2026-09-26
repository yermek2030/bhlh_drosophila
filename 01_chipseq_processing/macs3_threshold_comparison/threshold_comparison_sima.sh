#!/usr/bin/env bash
# ============================================================
# SIMA MACS3 threshold comparison
# Single-end, 2 replicates + 1 matched input control
# MACS3 -f BAM --nomodel --extsize (from SPP); thresholds 1e-3, 1e-2, 0.05, 1e-1
# Output: one summary row per threshold to stdout
# ============================================================

set -euo pipefail
mkdir -p 05_peaks_thresh/SIMA
mkdir -p 10_idr_thresh/SIMA

# --- SPP-derived: SIMA_rep1 (230 bp) + SIMA_rep2 (220 bp) -> mean = 225 bp ---
# Matches TRH Chapter 3 approach (TRH used --extsize 210 from 205/215 bp SPP means)
EXTSIZE=225

echo "THRESH | MACS_rep1 | MACS_rep2 | IDR_pass | IDR_mu | IDR_pct | FRIP_rep1 | FRIP_rep2"

# ============================================================
# Reusable helper: run MACS3 + IDR + FRiP for one threshold
# Args: $1 = threshold tag (e.g., "p1e3"), $2 = -p value (e.g., "1e-3")
# ============================================================
run_threshold() {
  local TAG=$1
  local PVAL=$2

  conda activate phd

  # --- MACS3 per-replicate (SE: -f BAM --nomodel --extsize) ---
  for REP in 1 2
  do
    macs3 callpeak \
      -t 07_blacklist/SIMA_rep${REP}.clean.bam \
      -c 07_blacklist/SIMA_control.clean.bam \
      -f BAM -g dm \
      -n SIMA_rep${REP}_${TAG} \
      --outdir 05_peaks_thresh/SIMA \
      --keep-dup all \
      --nomodel --extsize ${EXTSIZE} \
      -p ${PVAL} \
      --call-summits \
      2> 05_peaks_thresh/SIMA/SIMA_rep${REP}_${TAG}.log
  done

  local M1=$(wc -l < 05_peaks_thresh/SIMA/SIMA_rep1_${TAG}_peaks.narrowPeak)
  local M2=$(wc -l < 05_peaks_thresh/SIMA/SIMA_rep2_${TAG}_peaks.narrowPeak)

  # --- Sort by signal score (col 8) for IDR ---
  for REP in 1 2
  do
    sort -k8,8nr 05_peaks_thresh/SIMA/SIMA_rep${REP}_${TAG}_peaks.narrowPeak \
      > 10_idr_thresh/SIMA/SIMA_rep${REP}_${TAG}.sorted.narrowPeak
  done

  conda activate idr

  # --- IDR (single pair, n=2) ---
  idr \
    --samples \
      10_idr_thresh/SIMA/SIMA_rep1_${TAG}.sorted.narrowPeak \
      10_idr_thresh/SIMA/SIMA_rep2_${TAG}.sorted.narrowPeak \
    --input-file-type narrowPeak \
    --rank p.value --idr-threshold 0.05 \
    --output-file 10_idr_thresh/SIMA/SIMA_1v2_${TAG}.txt \
    --plot \
    --log-output-file 10_idr_thresh/SIMA/SIMA_1v2_${TAG}.log

  conda deactivate

  # --- Column-count-aware IDR peak extraction ---
  local NCOLS=$(awk '{print NF; exit}' 10_idr_thresh/SIMA/SIMA_1v2_${TAG}.txt)
  local IDR_PASS
  if [ "$NCOLS" -ge 12 ]; then
    IDR_PASS=$(awk '$12 >= 1.301' 10_idr_thresh/SIMA/SIMA_1v2_${TAG}.txt | wc -l)
    awk '$12 >= 1.301 {print $1"\t"$2"\t"$3}' \
      10_idr_thresh/SIMA/SIMA_1v2_${TAG}.txt \
      | sort -k1,1 -k2,2n | bedtools merge -i stdin \
      > /tmp/SIMA_idr_${TAG}.bed
  else
    IDR_PASS=$(awk '$5 >= 540' 10_idr_thresh/SIMA/SIMA_1v2_${TAG}.txt | wc -l)
    awk '$5 >= 540 {print $1"\t"$2"\t"$3}' \
      10_idr_thresh/SIMA/SIMA_1v2_${TAG}.txt \
      | sort -k1,1 -k2,2n | bedtools merge -i stdin \
      > /tmp/SIMA_idr_${TAG}.bed
  fi

  # --- Parse IDR log for MU and PCT ---
  local MU=$(grep "Final parameter values" 10_idr_thresh/SIMA/SIMA_1v2_${TAG}.log \
    | grep -oP '\[\K[0-9.]+' | head -1)
  local PCT=$(grep "Number of peaks passing" 10_idr_thresh/SIMA/SIMA_1v2_${TAG}.log \
    | grep -oP '[0-9]+\.[0-9]+(?=%)' | head -1)

  # --- FRiP per replicate ---
  local FRIP_R1 FRIP_R2
  for REP in 1 2
  do
    bedtools bamtobed -i 07_blacklist/SIMA_rep${REP}.clean.bam \
      > /tmp/SIMA_r${REP}_${TAG}.bed
    local TOTAL=$(wc -l < /tmp/SIMA_r${REP}_${TAG}.bed)
    local IN=$(bedtools intersect -u -a /tmp/SIMA_r${REP}_${TAG}.bed \
      -b /tmp/SIMA_idr_${TAG}.bed | wc -l)
    local FRIP=$(awk -v a="$IN" -v b="$TOTAL" 'BEGIN{printf "%.4f", a/b}')
    rm /tmp/SIMA_r${REP}_${TAG}.bed
    eval "FRIP_R${REP}=${FRIP}"
  done
  rm /tmp/SIMA_idr_${TAG}.bed

  echo "${TAG} | $M1 | $M2 | $IDR_PASS | $MU | ${PCT}% | $FRIP_R1 | $FRIP_R2"
}

# ============================================================
# Run all four thresholds
# ============================================================
run_threshold "p1e3"  "1e-3"
run_threshold "p1e2"  "1e-2"
run_threshold "p5e2"  "0.05"
run_threshold "p1e1"  "1e-1"

echo ""
echo "=== SIMA threshold comparison complete ==="
echo "All MACS3 peaks: 05_peaks_thresh/SIMA/"
echo "All IDR outputs: 10_idr_thresh/SIMA/"
echo ""
echo "Guidance: select the p-value that maximises IDR_pass count while keeping"
echo "IDR_mu >= 1.3. For n=2 replicates Sima may benefit from pseudoreplicate IDR"
echo "validation (see pseudoreplicate_idr_sima.sh) to guard against threshold-dependent"
echo "rank instability — the same concern that applied to TRH in Chapter 3."

# SIMA threshold comparison: MACS3 peak counts, IDR pass counts, IDR_mu and FRIP for p1e3/p1e2/p5e2/p1e1.
# MACS3 peaks: 05_peaks_thresh/SIMA/; IDR outputs: 10_idr_thresh/SIMA/
