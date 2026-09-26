#!/usr/bin/env bash
# ============================================================
# TGO MACS3 threshold comparison, paired-end, 3 reps + input control
# Pairwise IDR rep1v2, rep1v3, rep2v3; consensus = peaks in >=2/3 pairs
# Thresholds 1e-3, 1e-2, 1e-1; one summary row per threshold to stdout
# Run: bash threshold_comparison_tgo.sh 2> threshold_comparison_tgo.err
# ============================================================

set -euo pipefail
mkdir -p 05_peaks_thresh/TGO
mkdir -p 10_idr_thresh/TGO

echo "THRESH | MACS_r1 | MACS_r2 | MACS_r3 | IDR_1v2 | IDR_1v3 | IDR_2v3 | CONSENSUS | MU_1v2 | MU_1v3 | MU_2v3 | PCT_1v2 | PCT_1v3 | PCT_2v3 | FRIP_r1 | FRIP_r2 | FRIP_r3"
conda activate phd
# ============================================================
# THRESHOLD 1e-3
# ============================================================

# --- MACS3 ---
macs3 callpeak \
  -t 07_blacklist/TGO_rep1.clean.bam \
  -c 07_blacklist/TGO_control.clean.bam \
  -f BAMPE -g dm \
  -n TGO_rep1_p1e3 \
  --outdir 05_peaks_thresh/TGO \
  --keep-dup all \
  -p 1e-3 --call-summits \
  2> 05_peaks_thresh/TGO/TGO_rep1_p1e3.log

macs3 callpeak \
  -t 07_blacklist/TGO_rep2.clean.bam \
  -c 07_blacklist/TGO_control.clean.bam \
  -f BAMPE -g dm \
  -n TGO_rep2_p1e3 \
  --outdir 05_peaks_thresh/TGO \
  --keep-dup all \
  -p 1e-3 --call-summits \
  2> 05_peaks_thresh/TGO/TGO_rep2_p1e3.log

macs3 callpeak \
  -t 07_blacklist/TGO_rep3.clean.bam \
  -c 07_blacklist/TGO_control.clean.bam \
  -f BAMPE -g dm \
  -n TGO_rep3_p1e3 \
  --outdir 05_peaks_thresh/TGO \
  --keep-dup all \
  -p 1e-3 --call-summits \
  2> 05_peaks_thresh/TGO/TGO_rep3_p1e3.log

M1_1e3=$(wc -l < 05_peaks_thresh/TGO/TGO_rep1_p1e3_peaks.narrowPeak)
M2_1e3=$(wc -l < 05_peaks_thresh/TGO/TGO_rep2_p1e3_peaks.narrowPeak)
M3_1e3=$(wc -l < 05_peaks_thresh/TGO/TGO_rep3_p1e3_peaks.narrowPeak)

# --- Sort ---
sort -k8,8nr 05_peaks_thresh/TGO/TGO_rep1_p1e3_peaks.narrowPeak \
  > 10_idr_thresh/TGO/TGO_rep1_p1e3.sorted.narrowPeak
sort -k8,8nr 05_peaks_thresh/TGO/TGO_rep2_p1e3_peaks.narrowPeak \
  > 10_idr_thresh/TGO/TGO_rep2_p1e3.sorted.narrowPeak
sort -k8,8nr 05_peaks_thresh/TGO/TGO_rep3_p1e3_peaks.narrowPeak \
  > 10_idr_thresh/TGO/TGO_rep3_p1e3.sorted.narrowPeak
conda activate idr
# --- IDR rep1 vs rep2 ---
idr \
  --samples \
    10_idr_thresh/TGO/TGO_rep1_p1e3.sorted.narrowPeak \
    10_idr_thresh/TGO/TGO_rep2_p1e3.sorted.narrowPeak \
  --input-file-type narrowPeak \
  --rank p.value --idr-threshold 0.05 \
  --output-file 10_idr_thresh/TGO/TGO_1v2_p1e3.txt \
  --plot \
  --log-output-file 10_idr_thresh/TGO/TGO_1v2_p1e3.log

# --- IDR rep1 vs rep3 ---
idr \
  --samples \
    10_idr_thresh/TGO/TGO_rep1_p1e3.sorted.narrowPeak \
    10_idr_thresh/TGO/TGO_rep3_p1e3.sorted.narrowPeak \
  --input-file-type narrowPeak \
  --rank p.value --idr-threshold 0.05 \
  --output-file 10_idr_thresh/TGO/TGO_1v3_p1e3.txt \
  --plot \
  --log-output-file 10_idr_thresh/TGO/TGO_1v3_p1e3.log

# --- IDR rep2 vs rep3 ---
idr \
  --samples \
    10_idr_thresh/TGO/TGO_rep2_p1e3.sorted.narrowPeak \
    10_idr_thresh/TGO/TGO_rep3_p1e3.sorted.narrowPeak \
  --input-file-type narrowPeak \
  --rank p.value --idr-threshold 0.05 \
  --output-file 10_idr_thresh/TGO/TGO_2v3_p1e3.txt \
  --plot \
  --log-output-file 10_idr_thresh/TGO/TGO_2v3_p1e3.log
conda deactivate
# --- Extract metrics ---
NCOLS_1e3=$(awk '{print NF; exit}' 10_idr_thresh/TGO/TGO_1v2_p1e3.txt)
if [ "$NCOLS_1e3" -ge 12 ]; then
  IDR_12_1e3=$(awk '$12 >= 1.301' 10_idr_thresh/TGO/TGO_1v2_p1e3.txt | wc -l)
  IDR_13_1e3=$(awk '$12 >= 1.301' 10_idr_thresh/TGO/TGO_1v3_p1e3.txt | wc -l)
  IDR_23_1e3=$(awk '$12 >= 1.301' 10_idr_thresh/TGO/TGO_2v3_p1e3.txt | wc -l)
else
  IDR_12_1e3=$(awk '$5 >= 540' 10_idr_thresh/TGO/TGO_1v2_p1e3.txt | wc -l)
  IDR_13_1e3=$(awk '$5 >= 540' 10_idr_thresh/TGO/TGO_1v3_p1e3.txt | wc -l)
  IDR_23_1e3=$(awk '$5 >= 540' 10_idr_thresh/TGO/TGO_2v3_p1e3.txt | wc -l)
fi

MU_12_1e3=$(grep "Final parameter values" 10_idr_thresh/TGO/TGO_1v2_p1e3.log \
  | grep -oP '\[\K[0-9.]+' | head -1)
MU_13_1e3=$(grep "Final parameter values" 10_idr_thresh/TGO/TGO_1v3_p1e3.log \
  | grep -oP '\[\K[0-9.]+' | head -1)
MU_23_1e3=$(grep "Final parameter values" 10_idr_thresh/TGO/TGO_2v3_p1e3.log \
  | grep -oP '\[\K[0-9.]+' | head -1)

PCT_12_1e3=$(grep "Number of peaks passing" 10_idr_thresh/TGO/TGO_1v2_p1e3.log \
  | grep -oP '[0-9]+\.[0-9]+(?=%)' | head -1)
PCT_13_1e3=$(grep "Number of peaks passing" 10_idr_thresh/TGO/TGO_1v3_p1e3.log \
  | grep -oP '[0-9]+\.[0-9]+(?=%)' | head -1)
PCT_23_1e3=$(grep "Number of peaks passing" 10_idr_thresh/TGO/TGO_2v3_p1e3.log \
  | grep -oP '[0-9]+\.[0-9]+(?=%)' | head -1)

# --- Consensus >=2/3 ---
if [ "$NCOLS_1e3" -ge 12 ]; then
  awk '$12 >= 1.301 {print $1"\t"$2"\t"$3}' \
    10_idr_thresh/TGO/TGO_1v2_p1e3.txt \
    | sort -k1,1 -k2,2n | bedtools merge -i stdin \
    > /tmp/TGO_1v2_p1e3.merged.bed
  awk '$12 >= 1.301 {print $1"\t"$2"\t"$3}' \
    10_idr_thresh/TGO/TGO_1v3_p1e3.txt \
    | sort -k1,1 -k2,2n | bedtools merge -i stdin \
    > /tmp/TGO_1v3_p1e3.merged.bed
  awk '$12 >= 1.301 {print $1"\t"$2"\t"$3}' \
    10_idr_thresh/TGO/TGO_2v3_p1e3.txt \
    | sort -k1,1 -k2,2n | bedtools merge -i stdin \
    > /tmp/TGO_2v3_p1e3.merged.bed
else
  awk '$5 >= 540 {print $1"\t"$2"\t"$3}' \
    10_idr_thresh/TGO/TGO_1v2_p1e3.txt \
    | sort -k1,1 -k2,2n | bedtools merge -i stdin \
    > /tmp/TGO_1v2_p1e3.merged.bed
  awk '$5 >= 540 {print $1"\t"$2"\t"$3}' \
    10_idr_thresh/TGO/TGO_1v3_p1e3.txt \
    | sort -k1,1 -k2,2n | bedtools merge -i stdin \
    > /tmp/TGO_1v3_p1e3.merged.bed
  awk '$5 >= 540 {print $1"\t"$2"\t"$3}' \
    10_idr_thresh/TGO/TGO_2v3_p1e3.txt \
    | sort -k1,1 -k2,2n | bedtools merge -i stdin \
    > /tmp/TGO_2v3_p1e3.merged.bed
fi

CONSENSUS_1e3=$(cat \
  /tmp/TGO_1v2_p1e3.merged.bed \
  /tmp/TGO_1v3_p1e3.merged.bed \
  /tmp/TGO_2v3_p1e3.merged.bed \
  | sort -k1,1 -k2,2n \
  | bedtools merge -i stdin -c 1 -o count \
  | awk '$4 >= 2' | wc -l)

# Keep consensus BED for FRiP
cat \
  /tmp/TGO_1v2_p1e3.merged.bed \
  /tmp/TGO_1v3_p1e3.merged.bed \
  /tmp/TGO_2v3_p1e3.merged.bed \
  | sort -k1,1 -k2,2n \
  | bedtools merge -i stdin -c 1 -o count \
  | awk '$4 >= 2 {print $1"\t"$2"\t"$3}' \
  > /tmp/TGO_consensus_p1e3.bed
rm /tmp/TGO_1v2_p1e3.merged.bed /tmp/TGO_1v3_p1e3.merged.bed /tmp/TGO_2v3_p1e3.merged.bed

# --- FRiP ---
bedtools bamtobed -i 07_blacklist/TGO_rep1.clean.bam > /tmp/TGO_r1_p1e3.bed
TOTAL_R1=$(wc -l < /tmp/TGO_r1_p1e3.bed)
IN_R1=$(bedtools intersect -u -a /tmp/TGO_r1_p1e3.bed \
  -b /tmp/TGO_consensus_p1e3.bed | wc -l)
FRIP_R1_1e3=$(awk -v a="$IN_R1" -v b="$TOTAL_R1" 'BEGIN{printf "%.4f", a/b}')
rm /tmp/TGO_r1_p1e3.bed

bedtools bamtobed -i 07_blacklist/TGO_rep2.clean.bam > /tmp/TGO_r2_p1e3.bed
TOTAL_R2=$(wc -l < /tmp/TGO_r2_p1e3.bed)
IN_R2=$(bedtools intersect -u -a /tmp/TGO_r2_p1e3.bed \
  -b /tmp/TGO_consensus_p1e3.bed | wc -l)
FRIP_R2_1e3=$(awk -v a="$IN_R2" -v b="$TOTAL_R2" 'BEGIN{printf "%.4f", a/b}')
rm /tmp/TGO_r2_p1e3.bed

bedtools bamtobed -i 07_blacklist/TGO_rep3.clean.bam > /tmp/TGO_r3_p1e3.bed
TOTAL_R3=$(wc -l < /tmp/TGO_r3_p1e3.bed)
IN_R3=$(bedtools intersect -u -a /tmp/TGO_r3_p1e3.bed \
  -b /tmp/TGO_consensus_p1e3.bed | wc -l)
FRIP_R3_1e3=$(awk -v a="$IN_R3" -v b="$TOTAL_R3" 'BEGIN{printf "%.4f", a/b}')
rm /tmp/TGO_r3_p1e3.bed /tmp/TGO_consensus_p1e3.bed

echo "1e-3 | $M1_1e3 | $M2_1e3 | $M3_1e3 | $IDR_12_1e3 | $IDR_13_1e3 | $IDR_23_1e3 | $CONSENSUS_1e3 | $MU_12_1e3 | $MU_13_1e3 | $MU_23_1e3 | ${PCT_12_1e3}% | ${PCT_13_1e3}% | ${PCT_23_1e3}% | $FRIP_R1_1e3 | $FRIP_R2_1e3 | $FRIP_R3_1e3"

# ============================================================
# THRESHOLD 1e-2
# ============================================================

macs3 callpeak \
  -t 07_blacklist/TGO_rep1.clean.bam \
  -c 07_blacklist/TGO_control.clean.bam \
  -f BAMPE -g dm \
  -n TGO_rep1_p1e2 \
  --outdir 05_peaks_thresh/TGO \
  --keep-dup all \
  -p 1e-2 --call-summits \
  2> 05_peaks_thresh/TGO/TGO_rep1_p1e2.log

macs3 callpeak \
  -t 07_blacklist/TGO_rep2.clean.bam \
  -c 07_blacklist/TGO_control.clean.bam \
  -f BAMPE -g dm \
  -n TGO_rep2_p1e2 \
  --outdir 05_peaks_thresh/TGO \
  --keep-dup all \
  -p 1e-2 --call-summits \
  2> 05_peaks_thresh/TGO/TGO_rep2_p1e2.log

macs3 callpeak \
  -t 07_blacklist/TGO_rep3.clean.bam \
  -c 07_blacklist/TGO_control.clean.bam \
  -f BAMPE -g dm \
  -n TGO_rep3_p1e2 \
  --outdir 05_peaks_thresh/TGO \
  --keep-dup all \
  -p 1e-2 --call-summits \
  2> 05_peaks_thresh/TGO/TGO_rep3_p1e2.log

M1_1e2=$(wc -l < 05_peaks_thresh/TGO/TGO_rep1_p1e2_peaks.narrowPeak)
M2_1e2=$(wc -l < 05_peaks_thresh/TGO/TGO_rep2_p1e2_peaks.narrowPeak)
M3_1e2=$(wc -l < 05_peaks_thresh/TGO/TGO_rep3_p1e2_peaks.narrowPeak)

sort -k8,8nr 05_peaks_thresh/TGO/TGO_rep1_p1e2_peaks.narrowPeak \
  > 10_idr_thresh/TGO/TGO_rep1_p1e2.sorted.narrowPeak
sort -k8,8nr 05_peaks_thresh/TGO/TGO_rep2_p1e2_peaks.narrowPeak \
  > 10_idr_thresh/TGO/TGO_rep2_p1e2.sorted.narrowPeak
sort -k8,8nr 05_peaks_thresh/TGO/TGO_rep3_p1e2_peaks.narrowPeak \
  > 10_idr_thresh/TGO/TGO_rep3_p1e2.sorted.narrowPeak
conda activate idr
idr \
  --samples \
    10_idr_thresh/TGO/TGO_rep1_p1e2.sorted.narrowPeak \
    10_idr_thresh/TGO/TGO_rep2_p1e2.sorted.narrowPeak \
  --input-file-type narrowPeak \
  --rank p.value --idr-threshold 0.05 \
  --output-file 10_idr_thresh/TGO/TGO_1v2_p1e2.txt \
  --plot \
  --log-output-file 10_idr_thresh/TGO/TGO_1v2_p1e2.log

idr \
  --samples \
    10_idr_thresh/TGO/TGO_rep1_p1e2.sorted.narrowPeak \
    10_idr_thresh/TGO/TGO_rep3_p1e2.sorted.narrowPeak \
  --input-file-type narrowPeak \
  --rank p.value --idr-threshold 0.05 \
  --output-file 10_idr_thresh/TGO/TGO_1v3_p1e2.txt \
  --plot \
  --log-output-file 10_idr_thresh/TGO/TGO_1v3_p1e2.log

idr \
  --samples \
    10_idr_thresh/TGO/TGO_rep2_p1e2.sorted.narrowPeak \
    10_idr_thresh/TGO/TGO_rep3_p1e2.sorted.narrowPeak \
  --input-file-type narrowPeak \
  --rank p.value --idr-threshold 0.05 \
  --output-file 10_idr_thresh/TGO/TGO_2v3_p1e2.txt \
  --plot \
  --log-output-file 10_idr_thresh/TGO/TGO_2v3_p1e2.log
conda deactivate 
NCOLS_1e2=$(awk '{print NF; exit}' 10_idr_thresh/TGO/TGO_1v2_p1e2.txt)
if [ "$NCOLS_1e2" -ge 12 ]; then
  IDR_12_1e2=$(awk '$12 >= 1.301' 10_idr_thresh/TGO/TGO_1v2_p1e2.txt | wc -l)
  IDR_13_1e2=$(awk '$12 >= 1.301' 10_idr_thresh/TGO/TGO_1v3_p1e2.txt | wc -l)
  IDR_23_1e2=$(awk '$12 >= 1.301' 10_idr_thresh/TGO/TGO_2v3_p1e2.txt | wc -l)
else
  IDR_12_1e2=$(awk '$5 >= 540' 10_idr_thresh/TGO/TGO_1v2_p1e2.txt | wc -l)
  IDR_13_1e2=$(awk '$5 >= 540' 10_idr_thresh/TGO/TGO_1v3_p1e2.txt | wc -l)
  IDR_23_1e2=$(awk '$5 >= 540' 10_idr_thresh/TGO/TGO_2v3_p1e2.txt | wc -l)
fi

MU_12_1e2=$(grep "Final parameter values" 10_idr_thresh/TGO/TGO_1v2_p1e2.log \
  | grep -oP '\[\K[0-9.]+' | head -1)
MU_13_1e2=$(grep "Final parameter values" 10_idr_thresh/TGO/TGO_1v3_p1e2.log \
  | grep -oP '\[\K[0-9.]+' | head -1)
MU_23_1e2=$(grep "Final parameter values" 10_idr_thresh/TGO/TGO_2v3_p1e2.log \
  | grep -oP '\[\K[0-9.]+' | head -1)

PCT_12_1e2=$(grep "Number of peaks passing" 10_idr_thresh/TGO/TGO_1v2_p1e2.log \
  | grep -oP '[0-9]+\.[0-9]+(?=%)' | head -1)
PCT_13_1e2=$(grep "Number of peaks passing" 10_idr_thresh/TGO/TGO_1v3_p1e2.log \
  | grep -oP '[0-9]+\.[0-9]+(?=%)' | head -1)
PCT_23_1e2=$(grep "Number of peaks passing" 10_idr_thresh/TGO/TGO_2v3_p1e2.log \
  | grep -oP '[0-9]+\.[0-9]+(?=%)' | head -1)

if [ "$NCOLS_1e2" -ge 12 ]; then
  awk '$12 >= 1.301 {print $1"\t"$2"\t"$3}' \
    10_idr_thresh/TGO/TGO_1v2_p1e2.txt \
    | sort -k1,1 -k2,2n | bedtools merge -i stdin \
    > /tmp/TGO_1v2_p1e2.merged.bed
  awk '$12 >= 1.301 {print $1"\t"$2"\t"$3}' \
    10_idr_thresh/TGO/TGO_1v3_p1e2.txt \
    | sort -k1,1 -k2,2n | bedtools merge -i stdin \
    > /tmp/TGO_1v3_p1e2.merged.bed
  awk '$12 >= 1.301 {print $1"\t"$2"\t"$3}' \
    10_idr_thresh/TGO/TGO_2v3_p1e2.txt \
    | sort -k1,1 -k2,2n | bedtools merge -i stdin \
    > /tmp/TGO_2v3_p1e2.merged.bed
else
  awk '$5 >= 540 {print $1"\t"$2"\t"$3}' \
    10_idr_thresh/TGO/TGO_1v2_p1e2.txt \
    | sort -k1,1 -k2,2n | bedtools merge -i stdin \
    > /tmp/TGO_1v2_p1e2.merged.bed
  awk '$5 >= 540 {print $1"\t"$2"\t"$3}' \
    10_idr_thresh/TGO/TGO_1v3_p1e2.txt \
    | sort -k1,1 -k2,2n | bedtools merge -i stdin \
    > /tmp/TGO_1v3_p1e2.merged.bed
  awk '$5 >= 540 {print $1"\t"$2"\t"$3}' \
    10_idr_thresh/TGO/TGO_2v3_p1e2.txt \
    | sort -k1,1 -k2,2n | bedtools merge -i stdin \
    > /tmp/TGO_2v3_p1e2.merged.bed
fi

CONSENSUS_1e2=$(cat \
  /tmp/TGO_1v2_p1e2.merged.bed \
  /tmp/TGO_1v3_p1e2.merged.bed \
  /tmp/TGO_2v3_p1e2.merged.bed \
  | sort -k1,1 -k2,2n \
  | bedtools merge -i stdin -c 1 -o count \
  | awk '$4 >= 2' | wc -l)

cat \
  /tmp/TGO_1v2_p1e2.merged.bed \
  /tmp/TGO_1v3_p1e2.merged.bed \
  /tmp/TGO_2v3_p1e2.merged.bed \
  | sort -k1,1 -k2,2n \
  | bedtools merge -i stdin -c 1 -o count \
  | awk '$4 >= 2 {print $1"\t"$2"\t"$3}' \
  > /tmp/TGO_consensus_p1e2.bed
rm /tmp/TGO_1v2_p1e2.merged.bed /tmp/TGO_1v3_p1e2.merged.bed /tmp/TGO_2v3_p1e2.merged.bed

bedtools bamtobed -i 07_blacklist/TGO_rep1.clean.bam > /tmp/TGO_r1_p1e2.bed
TOTAL_R1=$(wc -l < /tmp/TGO_r1_p1e2.bed)
IN_R1=$(bedtools intersect -u -a /tmp/TGO_r1_p1e2.bed \
  -b /tmp/TGO_consensus_p1e2.bed | wc -l)
FRIP_R1_1e2=$(awk -v a="$IN_R1" -v b="$TOTAL_R1" 'BEGIN{printf "%.4f", a/b}')
rm /tmp/TGO_r1_p1e2.bed

bedtools bamtobed -i 07_blacklist/TGO_rep2.clean.bam > /tmp/TGO_r2_p1e2.bed
TOTAL_R2=$(wc -l < /tmp/TGO_r2_p1e2.bed)
IN_R2=$(bedtools intersect -u -a /tmp/TGO_r2_p1e2.bed \
  -b /tmp/TGO_consensus_p1e2.bed | wc -l)
FRIP_R2_1e2=$(awk -v a="$IN_R2" -v b="$TOTAL_R2" 'BEGIN{printf "%.4f", a/b}')
rm /tmp/TGO_r2_p1e2.bed

bedtools bamtobed -i 07_blacklist/TGO_rep3.clean.bam > /tmp/TGO_r3_p1e2.bed
TOTAL_R3=$(wc -l < /tmp/TGO_r3_p1e2.bed)
IN_R3=$(bedtools intersect -u -a /tmp/TGO_r3_p1e2.bed \
  -b /tmp/TGO_consensus_p1e2.bed | wc -l)
FRIP_R3_1e2=$(awk -v a="$IN_R3" -v b="$TOTAL_R3" 'BEGIN{printf "%.4f", a/b}')
rm /tmp/TGO_r3_p1e2.bed /tmp/TGO_consensus_p1e2.bed

echo "1e-2 | $M1_1e2 | $M2_1e2 | $M3_1e2 | $IDR_12_1e2 | $IDR_13_1e2 | $IDR_23_1e2 | $CONSENSUS_1e2 | $MU_12_1e2 | $MU_13_1e2 | $MU_23_1e2 | ${PCT_12_1e2}% | ${PCT_13_1e2}% | ${PCT_23_1e2}% | $FRIP_R1_1e2 | $FRIP_R2_1e2 | $FRIP_R3_1e2"

# ============================================================
# THRESHOLD 0.05
# ============================================================
macs3 callpeak \
  -t 07_blacklist/TGO_rep1.clean.bam \
  -c 07_blacklist/TGO_control.clean.bam \
  -f BAMPE -g dm \
  -n TGO_rep1_p5e2 \
  --outdir 05_peaks_thresh/TGO \
  --keep-dup all \
  -p 0.05 --call-summits \
  2> 05_peaks_thresh/TGO/TGO_rep1_p5e2.log

macs3 callpeak \
  -t 07_blacklist/TGO_rep2.clean.bam \
  -c 07_blacklist/TGO_control.clean.bam \
  -f BAMPE -g dm \
  -n TGO_rep2_p5e2 \
  --outdir 05_peaks_thresh/TGO \
  --keep-dup all \
  -p 0.05 --call-summits \
  2> 05_peaks_thresh/TGO/TGO_rep2_p5e2.log

macs3 callpeak \
  -t 07_blacklist/TGO_rep3.clean.bam \
  -c 07_blacklist/TGO_control.clean.bam \
  -f BAMPE -g dm \
  -n TGO_rep3_p5e2 \
  --outdir 05_peaks_thresh/TGO \
  --keep-dup all \
  -p 0.05 --call-summits \
  2> 05_peaks_thresh/TGO/TGO_rep3_p5e2.log

TGO_M1=$(wc -l < 05_peaks_thresh/TGO/TGO_rep1_p5e2_peaks.narrowPeak)
TGO_M2=$(wc -l < 05_peaks_thresh/TGO/TGO_rep2_p5e2_peaks.narrowPeak)
TGO_M3=$(wc -l < 05_peaks_thresh/TGO/TGO_rep3_p5e2_peaks.narrowPeak)

# --- Sort ---
sort -k8,8nr 05_peaks_thresh/TGO/TGO_rep1_p5e2_peaks.narrowPeak \
  > 10_idr_thresh/TGO/TGO_rep1_p5e2.sorted.narrowPeak
sort -k8,8nr 05_peaks_thresh/TGO/TGO_rep2_p5e2_peaks.narrowPeak \
  > 10_idr_thresh/TGO/TGO_rep2_p5e2.sorted.narrowPeak
sort -k8,8nr 05_peaks_thresh/TGO/TGO_rep3_p5e2_peaks.narrowPeak \
  > 10_idr_thresh/TGO/TGO_rep3_p5e2.sorted.narrowPeak
conda activate idr
# --- IDR ---
idr \
  --samples \
    10_idr_thresh/TGO/TGO_rep1_p5e2.sorted.narrowPeak \
    10_idr_thresh/TGO/TGO_rep2_p5e2.sorted.narrowPeak \
  --input-file-type narrowPeak \
  --rank p.value --idr-threshold 0.05 \
  --output-file 10_idr_thresh/TGO/TGO_1v2_p5e2.txt \
  --plot \
  --log-output-file 10_idr_thresh/TGO/TGO_1v2_p5e2.log

idr \
  --samples \
    10_idr_thresh/TGO/TGO_rep1_p5e2.sorted.narrowPeak \
    10_idr_thresh/TGO/TGO_rep3_p5e2.sorted.narrowPeak \
  --input-file-type narrowPeak \
  --rank p.value --idr-threshold 0.05 \
  --output-file 10_idr_thresh/TGO/TGO_1v3_p5e2.txt \
  --plot \
  --log-output-file 10_idr_thresh/TGO/TGO_1v3_p5e2.log

idr \
  --samples \
    10_idr_thresh/TGO/TGO_rep2_p5e2.sorted.narrowPeak \
    10_idr_thresh/TGO/TGO_rep3_p5e2.sorted.narrowPeak \
  --input-file-type narrowPeak \
  --rank p.value --idr-threshold 0.05 \
  --output-file 10_idr_thresh/TGO/TGO_2v3_p5e2.txt \
  --plot \
  --log-output-file 10_idr_thresh/TGO/TGO_2v3_p5e2.log
conda deactivate
# --- Metrics ---
NCOLS=$(awk '{print NF; exit}' 10_idr_thresh/TGO/TGO_1v2_p5e2.txt)
if [ "$NCOLS" -ge 12 ]; then
  TGO_12=$(awk '$12 >= 1.301' 10_idr_thresh/TGO/TGO_1v2_p5e2.txt | wc -l)
  TGO_13=$(awk '$12 >= 1.301' 10_idr_thresh/TGO/TGO_1v3_p5e2.txt | wc -l)
  TGO_23=$(awk '$12 >= 1.301' 10_idr_thresh/TGO/TGO_2v3_p5e2.txt | wc -l)
  awk '$12 >= 1.301 {print $1"\t"$2"\t"$3}' \
    10_idr_thresh/TGO/TGO_1v2_p5e2.txt \
    | sort -k1,1 -k2,2n | bedtools merge -i stdin \
    > /tmp/TGO_1v2_p5e2.merged.bed
  awk '$12 >= 1.301 {print $1"\t"$2"\t"$3}' \
    10_idr_thresh/TGO/TGO_1v3_p5e2.txt \
    | sort -k1,1 -k2,2n | bedtools merge -i stdin \
    > /tmp/TGO_1v3_p5e2.merged.bed
  awk '$12 >= 1.301 {print $1"\t"$2"\t"$3}' \
    10_idr_thresh/TGO/TGO_2v3_p5e2.txt \
    | sort -k1,1 -k2,2n | bedtools merge -i stdin \
    > /tmp/TGO_2v3_p5e2.merged.bed
else
  TGO_12=$(awk '$5 >= 540' 10_idr_thresh/TGO/TGO_1v2_p5e2.txt | wc -l)
  TGO_13=$(awk '$5 >= 540' 10_idr_thresh/TGO/TGO_1v3_p5e2.txt | wc -l)
  TGO_23=$(awk '$5 >= 540' 10_idr_thresh/TGO/TGO_2v3_p5e2.txt | wc -l)
  awk '$5 >= 540 {print $1"\t"$2"\t"$3}' \
    10_idr_thresh/TGO/TGO_1v2_p5e2.txt \
    | sort -k1,1 -k2,2n | bedtools merge -i stdin \
    > /tmp/TGO_1v2_p5e2.merged.bed
  awk '$5 >= 540 {print $1"\t"$2"\t"$3}' \
    10_idr_thresh/TGO/TGO_1v3_p5e2.txt \
    | sort -k1,1 -k2,2n | bedtools merge -i stdin \
    > /tmp/TGO_1v3_p5e2.merged.bed
  awk '$5 >= 540 {print $1"\t"$2"\t"$3}' \
    10_idr_thresh/TGO/TGO_2v3_p5e2.txt \
    | sort -k1,1 -k2,2n | bedtools merge -i stdin \
    > /tmp/TGO_2v3_p5e2.merged.bed
fi

TGO_MU_12=$(grep "Final parameter values" 10_idr_thresh/TGO/TGO_1v2_p5e2.log \
  | grep -oP '\[\K[0-9.]+' | head -1)
TGO_MU_13=$(grep "Final parameter values" 10_idr_thresh/TGO/TGO_1v3_p5e2.log \
  | grep -oP '\[\K[0-9.]+' | head -1)
TGO_MU_23=$(grep "Final parameter values" 10_idr_thresh/TGO/TGO_2v3_p5e2.log \
  | grep -oP '\[\K[0-9.]+' | head -1)

TGO_PCT_12=$(grep "Number of peaks passing" 10_idr_thresh/TGO/TGO_1v2_p5e2.log \
  | grep -oP '[0-9]+\.[0-9]+(?=%)' | head -1)
TGO_PCT_13=$(grep "Number of peaks passing" 10_idr_thresh/TGO/TGO_1v3_p5e2.log \
  | grep -oP '[0-9]+\.[0-9]+(?=%)' | head -1)
TGO_PCT_23=$(grep "Number of peaks passing" 10_idr_thresh/TGO/TGO_2v3_p5e2.log \
  | grep -oP '[0-9]+\.[0-9]+(?=%)' | head -1)

TGO_CONS=$(cat \
  /tmp/TGO_1v2_p5e2.merged.bed \
  /tmp/TGO_1v3_p5e2.merged.bed \
  /tmp/TGO_2v3_p5e2.merged.bed \
  | sort -k1,1 -k2,2n \
  | bedtools merge -i stdin -c 1 -o count \
  | awk '$4 >= 2' | wc -l)

cat \
  /tmp/TGO_1v2_p5e2.merged.bed \
  /tmp/TGO_1v3_p5e2.merged.bed \
  /tmp/TGO_2v3_p5e2.merged.bed \
  | sort -k1,1 -k2,2n \
  | bedtools merge -i stdin -c 1 -o count \
  | awk '$4 >= 2 {print $1"\t"$2"\t"$3}' \
  > /tmp/TGO_consensus_p5e2.bed
rm /tmp/TGO_1v2_p5e2.merged.bed \
   /tmp/TGO_1v3_p5e2.merged.bed \
   /tmp/TGO_2v3_p5e2.merged.bed

# --- FRiP ---
bedtools bamtobed -i 07_blacklist/TGO_rep1.clean.bam > /tmp/TGO_r1.bed
TOTAL=$(wc -l < /tmp/TGO_r1.bed)
IN=$(bedtools intersect -u -a /tmp/TGO_r1.bed \
  -b /tmp/TGO_consensus_p5e2.bed | wc -l)
TGO_FRIP1=$(awk -v a="$IN" -v b="$TOTAL" 'BEGIN{printf "%.4f", a/b}')
rm /tmp/TGO_r1.bed

bedtools bamtobed -i 07_blacklist/TGO_rep2.clean.bam > /tmp/TGO_r2.bed
TOTAL=$(wc -l < /tmp/TGO_r2.bed)
IN=$(bedtools intersect -u -a /tmp/TGO_r2.bed \
  -b /tmp/TGO_consensus_p5e2.bed | wc -l)
TGO_FRIP2=$(awk -v a="$IN" -v b="$TOTAL" 'BEGIN{printf "%.4f", a/b}')
rm /tmp/TGO_r2.bed

bedtools bamtobed -i 07_blacklist/TGO_rep3.clean.bam > /tmp/TGO_r3.bed
TOTAL=$(wc -l < /tmp/TGO_r3.bed)
IN=$(bedtools intersect -u -a /tmp/TGO_r3.bed \
  -b /tmp/TGO_consensus_p5e2.bed | wc -l)
TGO_FRIP3=$(awk -v a="$IN" -v b="$TOTAL" 'BEGIN{printf "%.4f", a/b}')
rm /tmp/TGO_r3.bed /tmp/TGO_consensus_p5e2.bed

echo "THRESH | MACS_r1 | MACS_r2 | MACS_r3 | IDR_1v2 | IDR_1v3 | IDR_2v3 | CONSENSUS | MU_1v2 | MU_1v3 | MU_2v3 | PCT_1v2 | PCT_1v3 | PCT_2v3 | FRIP_r1 | FRIP_r2 | FRIP_r3"
echo "0.05   | $TGO_M1 | $TGO_M2 | $TGO_M3 | $TGO_12 | $TGO_13 | $TGO_23 | $TGO_CONS | $TGO_MU_12 | $TGO_MU_13 | $TGO_MU_23 | ${TGO_PCT_12}% | ${TGO_PCT_13}% | ${TGO_PCT_23}% | $TGO_FRIP1 | $TGO_FRIP2 | $TGO_FRIP3"
echo ""
echo "--- TGO context (from previous runs) ---"
echo "1e-3 | 4317  | 6606  | 4017  | 1531 | 1096 | 1114 | 1052 | 1.44 | 1.14 | 1.30 | 46.5% | 40.7% | 37.4% | 0.0257 | 0.0281 | 0.0229"
echo "1e-2 | 8128  | 11903 | 7894  | 1700 | 1781 | 1407 | 1403 | 1.48 | 1.59 | 1.35 | 36.3% | 44.6% | 31.1% | 0.0296 | 0.0327 | 0.0269"
echo "0.05 | (above)"
echo "1e-1 | 35434 | 36750 | 34288 | 1739 | 2114 | 1423 | 1265 | 1.73 | 1.95 | 1.77 | 13.7% | 15.6% | 10.8% | 0.0235 | 0.0260 | 0.0218"

echo ""
echo "=== p=0.05 test complete ==="

# ============================================================
# THRESHOLD 1e-1
# ============================================================

macs3 callpeak \
  -t 07_blacklist/TGO_rep1.clean.bam \
  -c 07_blacklist/TGO_control.clean.bam \
  -f BAMPE -g dm \
  -n TGO_rep1_p1e1 \
  --outdir 05_peaks_thresh/TGO \
  --keep-dup all \
  -p 1e-1 --call-summits \
  2> 05_peaks_thresh/TGO/TGO_rep1_p1e1.log

macs3 callpeak \
  -t 07_blacklist/TGO_rep2.clean.bam \
  -c 07_blacklist/TGO_control.clean.bam \
  -f BAMPE -g dm \
  -n TGO_rep2_p1e1 \
  --outdir 05_peaks_thresh/TGO \
  --keep-dup all \
  -p 1e-1 --call-summits \
  2> 05_peaks_thresh/TGO/TGO_rep2_p1e1.log

macs3 callpeak \
  -t 07_blacklist/TGO_rep3.clean.bam \
  -c 07_blacklist/TGO_control.clean.bam \
  -f BAMPE -g dm \
  -n TGO_rep3_p1e1 \
  --outdir 05_peaks_thresh/TGO \
  --keep-dup all \
  -p 1e-1 --call-summits \
  2> 05_peaks_thresh/TGO/TGO_rep3_p1e1.log

M1_1e1=$(wc -l < 05_peaks_thresh/TGO/TGO_rep1_p1e1_peaks.narrowPeak)
M2_1e1=$(wc -l < 05_peaks_thresh/TGO/TGO_rep2_p1e1_peaks.narrowPeak)
M3_1e1=$(wc -l < 05_peaks_thresh/TGO/TGO_rep3_p1e1_peaks.narrowPeak)

sort -k8,8nr 05_peaks_thresh/TGO/TGO_rep1_p1e1_peaks.narrowPeak \
  > 10_idr_thresh/TGO/TGO_rep1_p1e1.sorted.narrowPeak
sort -k8,8nr 05_peaks_thresh/TGO/TGO_rep2_p1e1_peaks.narrowPeak \
  > 10_idr_thresh/TGO/TGO_rep2_p1e1.sorted.narrowPeak
sort -k8,8nr 05_peaks_thresh/TGO/TGO_rep3_p1e1_peaks.narrowPeak \
  > 10_idr_thresh/TGO/TGO_rep3_p1e1.sorted.narrowPeak
conda activate idr
idr \
  --samples \
    10_idr_thresh/TGO/TGO_rep1_p1e1.sorted.narrowPeak \
    10_idr_thresh/TGO/TGO_rep2_p1e1.sorted.narrowPeak \
  --input-file-type narrowPeak \
  --rank p.value --idr-threshold 0.05 \
  --output-file 10_idr_thresh/TGO/TGO_1v2_p1e1.txt \
  --plot \
  --log-output-file 10_idr_thresh/TGO/TGO_1v2_p1e1.log

idr \
  --samples \
    10_idr_thresh/TGO/TGO_rep1_p1e1.sorted.narrowPeak \
    10_idr_thresh/TGO/TGO_rep3_p1e1.sorted.narrowPeak \
  --input-file-type narrowPeak \
  --rank p.value --idr-threshold 0.05 \
  --output-file 10_idr_thresh/TGO/TGO_1v3_p1e1.txt \
  --plot \
  --log-output-file 10_idr_thresh/TGO/TGO_1v3_p1e1.log

idr \
  --samples \
    10_idr_thresh/TGO/TGO_rep2_p1e1.sorted.narrowPeak \
    10_idr_thresh/TGO/TGO_rep3_p1e1.sorted.narrowPeak \
  --input-file-type narrowPeak \
  --rank p.value --idr-threshold 0.05 \
  --output-file 10_idr_thresh/TGO/TGO_2v3_p1e1.txt \
  --plot \
  --log-output-file 10_idr_thresh/TGO/TGO_2v3_p1e1.log
conda deactivate
NCOLS_1e1=$(awk '{print NF; exit}' 10_idr_thresh/TGO/TGO_1v2_p1e1.txt)
if [ "$NCOLS_1e1" -ge 12 ]; then
  IDR_12_1e1=$(awk '$12 >= 1.301' 10_idr_thresh/TGO/TGO_1v2_p1e1.txt | wc -l)
  IDR_13_1e1=$(awk '$12 >= 1.301' 10_idr_thresh/TGO/TGO_1v3_p1e1.txt | wc -l)
  IDR_23_1e1=$(awk '$12 >= 1.301' 10_idr_thresh/TGO/TGO_2v3_p1e1.txt | wc -l)
else
  IDR_12_1e1=$(awk '$5 >= 540' 10_idr_thresh/TGO/TGO_1v2_p1e1.txt | wc -l)
  IDR_13_1e1=$(awk '$5 >= 540' 10_idr_thresh/TGO/TGO_1v3_p1e1.txt | wc -l)
  IDR_23_1e1=$(awk '$5 >= 540' 10_idr_thresh/TGO/TGO_2v3_p1e1.txt | wc -l)
fi

MU_12_1e1=$(grep "Final parameter values" 10_idr_thresh/TGO/TGO_1v2_p1e1.log \
  | grep -oP '\[\K[0-9.]+' | head -1)
MU_13_1e1=$(grep "Final parameter values" 10_idr_thresh/TGO/TGO_1v3_p1e1.log \
  | grep -oP '\[\K[0-9.]+' | head -1)
MU_23_1e1=$(grep "Final parameter values" 10_idr_thresh/TGO/TGO_2v3_p1e1.log \
  | grep -oP '\[\K[0-9.]+' | head -1)

PCT_12_1e1=$(grep "Number of peaks passing" 10_idr_thresh/TGO/TGO_1v2_p1e1.log \
  | grep -oP '[0-9]+\.[0-9]+(?=%)' | head -1)
PCT_13_1e1=$(grep "Number of peaks passing" 10_idr_thresh/TGO/TGO_1v3_p1e1.log \
  | grep -oP '[0-9]+\.[0-9]+(?=%)' | head -1)
PCT_23_1e1=$(grep "Number of peaks passing" 10_idr_thresh/TGO/TGO_2v3_p1e1.log \
  | grep -oP '[0-9]+\.[0-9]+(?=%)' | head -1)

if [ "$NCOLS_1e1" -ge 12 ]; then
  awk '$12 >= 1.301 {print $1"\t"$2"\t"$3}' \
    10_idr_thresh/TGO/TGO_1v2_p1e1.txt \
    | sort -k1,1 -k2,2n | bedtools merge -i stdin \
    > /tmp/TGO_1v2_p1e1.merged.bed
  awk '$12 >= 1.301 {print $1"\t"$2"\t"$3}' \
    10_idr_thresh/TGO/TGO_1v3_p1e1.txt \
    | sort -k1,1 -k2,2n | bedtools merge -i stdin \
    > /tmp/TGO_1v3_p1e1.merged.bed
  awk '$12 >= 1.301 {print $1"\t"$2"\t"$3}' \
    10_idr_thresh/TGO/TGO_2v3_p1e1.txt \
    | sort -k1,1 -k2,2n | bedtools merge -i stdin \
    > /tmp/TGO_2v3_p1e1.merged.bed
else
  awk '$5 >= 540 {print $1"\t"$2"\t"$3}' \
    10_idr_thresh/TGO/TGO_1v2_p1e1.txt \
    | sort -k1,1 -k2,2n | bedtools merge -i stdin \
    > /tmp/TGO_1v2_p1e1.merged.bed
  awk '$5 >= 540 {print $1"\t"$2"\t"$3}' \
    10_idr_thresh/TGO/TGO_1v3_p1e1.txt \
    | sort -k1,1 -k2,2n | bedtools merge -i stdin \
    > /tmp/TGO_1v3_p1e1.merged.bed
  awk '$5 >= 540 {print $1"\t"$2"\t"$3}' \
    10_idr_thresh/TGO/TGO_2v3_p1e1.txt \
    | sort -k1,1 -k2,2n | bedtools merge -i stdin \
    > /tmp/TGO_2v3_p1e1.merged.bed
fi

CONSENSUS_1e1=$(cat \
  /tmp/TGO_1v2_p1e1.merged.bed \
  /tmp/TGO_1v3_p1e1.merged.bed \
  /tmp/TGO_2v3_p1e1.merged.bed \
  | sort -k1,1 -k2,2n \
  | bedtools merge -i stdin -c 1 -o count \
  | awk '$4 >= 2' | wc -l)

cat \
  /tmp/TGO_1v2_p1e1.merged.bed \
  /tmp/TGO_1v3_p1e1.merged.bed \
  /tmp/TGO_2v3_p1e1.merged.bed \
  | sort -k1,1 -k2,2n \
  | bedtools merge -i stdin -c 1 -o count \
  | awk '$4 >= 2 {print $1"\t"$2"\t"$3}' \
  > /tmp/TGO_consensus_p1e1.bed
rm /tmp/TGO_1v2_p1e1.merged.bed /tmp/TGO_1v3_p1e1.merged.bed /tmp/TGO_2v3_p1e1.merged.bed

bedtools bamtobed -i 07_blacklist/TGO_rep1.clean.bam > /tmp/TGO_r1_p1e1.bed
TOTAL_R1=$(wc -l < /tmp/TGO_r1_p1e1.bed)
IN_R1=$(bedtools intersect -u -a /tmp/TGO_r1_p1e1.bed \
  -b /tmp/TGO_consensus_p1e1.bed | wc -l)
FRIP_R1_1e1=$(awk -v a="$IN_R1" -v b="$TOTAL_R1" 'BEGIN{printf "%.4f", a/b}')
rm /tmp/TGO_r1_p1e1.bed

bedtools bamtobed -i 07_blacklist/TGO_rep3.clean.bam > /tmp/TGO_r3_p1e1.bed
TOTAL_R3=$(wc -l < /tmp/TGO_r3_p1e1.bed)
IN_R3=$(bedtools intersect -u -a /tmp/TGO_r3_p1e1.bed \
  -b /tmp/TGO_consensus_p1e1.bed | wc -l)
FRIP_R3_1e1=$(awk -v a="$IN_R3" -v b="$TOTAL_R3" 'BEGIN{printf "%.4f", a/b}')
rm /tmp/TGO_r3_p1e1.bed

bedtools bamtobed -i 07_blacklist/TGO_rep2.clean.bam > /tmp/TGO_r2_p1e1.bed
TOTAL_R2=$(wc -l < /tmp/TGO_r2_p1e1.bed)
IN_R2=$(bedtools intersect -u -a /tmp/TGO_r2_p1e1.bed \
  -b /tmp/TGO_consensus_p1e1.bed | wc -l)
FRIP_R2_1e1=$(awk -v a="$IN_R2" -v b="$TOTAL_R2" 'BEGIN{printf "%.4f", a/b}')
rm /tmp/TGO_r2_p1e1.bed /tmp/TGO_consensus_p1e1.bed

echo "1e-1 | $M1_1e1 | $M2_1e1 | $M3_1e1 | $IDR_12_1e1 | $IDR_13_1e1 | $IDR_23_1e1 | $CONSENSUS_1e1 | $MU_12_1e1 | $MU_13_1e1 | $MU_23_1e1 | ${PCT_12_1e1}% | ${PCT_13_1e1}% | ${PCT_23_1e1}% | $FRIP_R1_1e1 | $FRIP_R2_1e1 | $FRIP_R3_1e1"

echo ""
echo "=== TGO threshold comparison complete ==="
echo "All MACS3 peaks: 05_peaks_thresh/TGO/"
echo "All IDR outputs: 10_idr_thresh/TGO/"
