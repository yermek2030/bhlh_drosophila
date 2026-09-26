#!/usr/bin/env bash
# ============================================================
# DYS MACS3 threshold comparison
# Paired-end, 3 replicates + 1 input control; MACS3 -f BAMPE
# Thresholds: 1e-3, 1e-2, 0.05, 1e-1; prints one summary row per threshold
# n=3 replicates: pairwise IDR (1v2,1v3,2v3), consensus = peaks in >=2 of 3 pairs
# ============================================================

set -euo pipefail
mkdir -p 05_peaks_thresh/DYS
mkdir -p 10_idr_thresh/DYS

echo "THRESH | MACS_r1 | MACS_r2 | MACS_r3 | IDR_12 | IDR_13 | IDR_23 | CONSENSUS | MU_12 | MU_13 | MU_23 | PCT_12 | PCT_13 | PCT_23 | FRIP_r1 | FRIP_r2 | FRIP_r3"

# ============================================================
# Reusable helper: run MACS3 + pairwise IDR + consensus + FRiP for one threshold
# Args: $1 = threshold tag (e.g., "p1e3"), $2 = -p value (e.g., "1e-3")
# ============================================================
run_threshold() {
  local TAG=$1
  local PVAL=$2

  conda activate phd

  # --- MACS3 per-replicate ---
  for REP in 1 2 3
  do
    macs3 callpeak \
      -t 07_blacklist/DYS_rep${REP}.clean.bam \
      -c 07_blacklist/DYS_control.clean.bam \
      -f BAMPE -g dm \
      -n DYS_rep${REP}_${TAG} \
      --outdir 05_peaks_thresh/DYS \
      --keep-dup all \
      -p ${PVAL} \
      --call-summits \
      2> 05_peaks_thresh/DYS/DYS_rep${REP}_${TAG}.log
  done

  local M1=$(wc -l < 05_peaks_thresh/DYS/DYS_rep1_${TAG}_peaks.narrowPeak)
  local M2=$(wc -l < 05_peaks_thresh/DYS/DYS_rep2_${TAG}_peaks.narrowPeak)
  local M3=$(wc -l < 05_peaks_thresh/DYS/DYS_rep3_${TAG}_peaks.narrowPeak)

  # --- Sort by signal score (col 8) for IDR ---
  for REP in 1 2 3
  do
    sort -k8,8nr 05_peaks_thresh/DYS/DYS_rep${REP}_${TAG}_peaks.narrowPeak \
      > 10_idr_thresh/DYS/DYS_rep${REP}_${TAG}.sorted.narrowPeak
  done

  conda activate idr

  # --- Pairwise IDR ---
  idr \
    --samples \
      10_idr_thresh/DYS/DYS_rep1_${TAG}.sorted.narrowPeak \
      10_idr_thresh/DYS/DYS_rep2_${TAG}.sorted.narrowPeak \
    --input-file-type narrowPeak \
    --rank p.value --idr-threshold 0.05 \
    --output-file 10_idr_thresh/DYS/DYS_1v2_${TAG}.txt \
    --plot \
    --log-output-file 10_idr_thresh/DYS/DYS_1v2_${TAG}.log

  idr \
    --samples \
      10_idr_thresh/DYS/DYS_rep1_${TAG}.sorted.narrowPeak \
      10_idr_thresh/DYS/DYS_rep3_${TAG}.sorted.narrowPeak \
    --input-file-type narrowPeak \
    --rank p.value --idr-threshold 0.05 \
    --output-file 10_idr_thresh/DYS/DYS_1v3_${TAG}.txt \
    --plot \
    --log-output-file 10_idr_thresh/DYS/DYS_1v3_${TAG}.log

  idr \
    --samples \
      10_idr_thresh/DYS/DYS_rep2_${TAG}.sorted.narrowPeak \
      10_idr_thresh/DYS/DYS_rep3_${TAG}.sorted.narrowPeak \
    --input-file-type narrowPeak \
    --rank p.value --idr-threshold 0.05 \
    --output-file 10_idr_thresh/DYS/DYS_2v3_${TAG}.txt \
    --plot \
    --log-output-file 10_idr_thresh/DYS/DYS_2v3_${TAG}.log

  conda deactivate

  # --- Column-count-aware IDR peak extraction ---
  local NCOLS=$(awk '{print NF; exit}' 10_idr_thresh/DYS/DYS_1v2_${TAG}.txt)
  local FILTER
  if [ "$NCOLS" -ge 12 ]; then
    FILTER='$12 >= 1.301'
  else
    FILTER='$5 >= 540'
  fi

  for PAIR in 1v2 1v3 2v3
  do
    awk "${FILTER} {print \$1\"\t\"\$2\"\t\"\$3}" \
      10_idr_thresh/DYS/DYS_${PAIR}_${TAG}.txt \
      | sort -k1,1 -k2,2n | bedtools merge -i stdin \
      > 10_idr_thresh/DYS/DYS_${PAIR}_${TAG}.merged.bed
  done

  local IDR_12=$(wc -l < 10_idr_thresh/DYS/DYS_1v2_${TAG}.merged.bed)
  local IDR_13=$(wc -l < 10_idr_thresh/DYS/DYS_1v3_${TAG}.merged.bed)
  local IDR_23=$(wc -l < 10_idr_thresh/DYS/DYS_2v3_${TAG}.merged.bed)

  # --- Consensus (>= 2/3 pairs) ---
  cat \
    10_idr_thresh/DYS/DYS_1v2_${TAG}.merged.bed \
    10_idr_thresh/DYS/DYS_1v3_${TAG}.merged.bed \
    10_idr_thresh/DYS/DYS_2v3_${TAG}.merged.bed \
    | sort -k1,1 -k2,2n \
    | bedtools merge -i stdin -c 1 -o count \
    | awk '$4 >= 2 {print $1"\t"$2"\t"$3}' \
    > /tmp/DYS_consensus_${TAG}.bed

  local CONSENSUS=$(wc -l < /tmp/DYS_consensus_${TAG}.bed)

  # --- Parse IDR log for MU and PCT ---
  local MU_12=$(grep "Final parameter values" 10_idr_thresh/DYS/DYS_1v2_${TAG}.log \
    | grep -oP '\[\K[0-9.]+' | head -1)
  local MU_13=$(grep "Final parameter values" 10_idr_thresh/DYS/DYS_1v3_${TAG}.log \
    | grep -oP '\[\K[0-9.]+' | head -1)
  local MU_23=$(grep "Final parameter values" 10_idr_thresh/DYS/DYS_2v3_${TAG}.log \
    | grep -oP '\[\K[0-9.]+' | head -1)

  local PCT_12=$(grep "Number of peaks passing" 10_idr_thresh/DYS/DYS_1v2_${TAG}.log \
    | grep -oP '[0-9]+\.[0-9]+(?=%)' | head -1)
  local PCT_13=$(grep "Number of peaks passing" 10_idr_thresh/DYS/DYS_1v3_${TAG}.log \
    | grep -oP '[0-9]+\.[0-9]+(?=%)' | head -1)
  local PCT_23=$(grep "Number of peaks passing" 10_idr_thresh/DYS/DYS_2v3_${TAG}.log \
    | grep -oP '[0-9]+\.[0-9]+(?=%)' | head -1)

  # --- FRiP per replicate against consensus peak set ---
  # PE: bedtools bamtobed outputs one line per read (not per pair);
  # ratio still meaningful because numerator and denominator use same counting
  local FRIP_R1 FRIP_R2 FRIP_R3
  for REP in 1 2 3
  do
    bedtools bamtobed -i 07_blacklist/DYS_rep${REP}.clean.bam \
      > /tmp/DYS_r${REP}_${TAG}.bed
    local TOTAL=$(wc -l < /tmp/DYS_r${REP}_${TAG}.bed)
    local IN=$(bedtools intersect -u -a /tmp/DYS_r${REP}_${TAG}.bed \
      -b /tmp/DYS_consensus_${TAG}.bed | wc -l)
    local FRIP=$(awk -v a="$IN" -v b="$TOTAL" 'BEGIN{printf "%.4f", a/b}')
    rm /tmp/DYS_r${REP}_${TAG}.bed
    eval "FRIP_R${REP}=${FRIP}"
  done
  rm /tmp/DYS_consensus_${TAG}.bed

  echo "${TAG} | $M1 | $M2 | $M3 | $IDR_12 | $IDR_13 | $IDR_23 | $CONSENSUS | $MU_12 | $MU_13 | $MU_23 | ${PCT_12}% | ${PCT_13}% | ${PCT_23}% | $FRIP_R1 | $FRIP_R2 | $FRIP_R3"
}

# ============================================================
# Run all four thresholds
# ============================================================
run_threshold "p1e3"  "1e-3"
run_threshold "p1e2"  "1e-2"
run_threshold "p5e2"  "0.05"
run_threshold "p1e1"  "1e-1"

echo ""
echo "=== DYS threshold comparison complete ==="
echo "All MACS3 peaks: 05_peaks_thresh/DYS/"
echo "All IDR outputs: 10_idr_thresh/DYS/"
echo ""
echo "Guidance: select the p-value that maximises CONSENSUS peak count while"
echo "keeping IDR_mu >= 1.3 and per-pair PCT >= 30%. If all four thresholds"
echo "produce similar CONSENSUS counts, prefer the more stringent (smaller -p)."

# Results in results/reported_values.tsv: MACS/IDR/consensus counts and FRIP per threshold (p1e3, p1e2, p5e2, p1e1).
# Peaks in 05_peaks_thresh/DYS/, IDR outputs in 10_idr_thresh/DYS/.
