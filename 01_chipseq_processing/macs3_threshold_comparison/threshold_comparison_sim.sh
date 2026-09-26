#!/usr/bin/env bash
# ============================================================
# SIM MACS3 threshold comparison
# Paired-end, 2 replicates + matched input control, -f BAMPE
# Thresholds: 1e-3, 1e-2, 1e-1; one summary row per threshold to stdout
# Run: bash threshold_comparison_sim.sh 2> threshold_comparison_sim.err
# ============================================================

set -euo pipefail
mkdir -p 05_peaks_thresh/SIM
mkdir -p 10_idr_thresh/SIM

echo "THRESH | MACS_rep1 | MACS_rep2 | IDR_pass | IDR_mu | IDR_pct | FRIP_rep1 | FRIP_rep2"

# ============================================================
# THRESHOLD 1e-3
# ============================================================
conda activate phd
macs3 callpeak \
  -t 07_blacklist/SIM_rep1.clean.bam \
  -c 07_blacklist/SIM_control.clean.bam \
  -f BAMPE -g dm \
  -n SIM_rep1_p1e3 \
  --outdir 05_peaks_thresh/SIM \
  --keep-dup all \
  -p 1e-3 \
  --call-summits \
  2> 05_peaks_thresh/SIM/SIM_rep1_p1e3.log

macs3 callpeak \
  -t 07_blacklist/SIM_rep2.clean.bam \
  -c 07_blacklist/SIM_control.clean.bam \
  -f BAMPE -g dm \
  -n SIM_rep2_p1e3 \
  --outdir 05_peaks_thresh/SIM \
  --keep-dup all \
  -p 1e-3 \
  --call-summits \
  2> 05_peaks_thresh/SIM/SIM_rep2_p1e3.log

M1_1e3=$(wc -l < 05_peaks_thresh/SIM/SIM_rep1_p1e3_peaks.narrowPeak)
M2_1e3=$(wc -l < 05_peaks_thresh/SIM/SIM_rep2_p1e3_peaks.narrowPeak)

sort -k8,8nr 05_peaks_thresh/SIM/SIM_rep1_p1e3_peaks.narrowPeak \
  > 10_idr_thresh/SIM/SIM_rep1_p1e3.sorted.narrowPeak

sort -k8,8nr 05_peaks_thresh/SIM/SIM_rep2_p1e3_peaks.narrowPeak \
  > 10_idr_thresh/SIM/SIM_rep2_p1e3.sorted.narrowPeak
conda activate idr
idr \
  --samples \
    10_idr_thresh/SIM/SIM_rep1_p1e3.sorted.narrowPeak \
    10_idr_thresh/SIM/SIM_rep2_p1e3.sorted.narrowPeak \
  --input-file-type narrowPeak \
  --rank p.value \
  --idr-threshold 0.05 \
  --output-file 10_idr_thresh/SIM/SIM_1v2_p1e3.txt \
  --plot \
  --log-output-file 10_idr_thresh/SIM/SIM_1v2_p1e3.log
conda deactivate
NCOLS_1e3=$(awk '{print NF; exit}' 10_idr_thresh/SIM/SIM_1v2_p1e3.txt)
if [ "$NCOLS_1e3" -ge 12 ]; then
  IDR_1e3=$(awk '$12 >= 1.301' 10_idr_thresh/SIM/SIM_1v2_p1e3.txt | wc -l)
else
  IDR_1e3=$(awk '$5 >= 540'   10_idr_thresh/SIM/SIM_1v2_p1e3.txt | wc -l)
fi

MU_1e3=$(grep "Final parameter values" 10_idr_thresh/SIM/SIM_1v2_p1e3.log \
  | grep -oP '\[\K[0-9.]+' | head -1)
PCT_1e3=$(grep "Number of peaks passing" 10_idr_thresh/SIM/SIM_1v2_p1e3.log \
  | grep -oP '[0-9]+\.[0-9]+(?=%)' | head -1)

if [ "$NCOLS_1e3" -ge 12 ]; then
  awk '$12 >= 1.301 {print $1"\t"$2"\t"$3}' \
    10_idr_thresh/SIM/SIM_1v2_p1e3.txt \
    | sort -k1,1 -k2,2n | bedtools merge -i stdin \
    > /tmp/SIM_idr_p1e3.bed
else
  awk '$5 >= 540 {print $1"\t"$2"\t"$3}' \
    10_idr_thresh/SIM/SIM_1v2_p1e3.txt \
    | sort -k1,1 -k2,2n | bedtools merge -i stdin \
    > /tmp/SIM_idr_p1e3.bed
fi

# FRiP for PE: bamtobed outputs one line per read (not per pair)
# Total read count and intersect count are still comparable as a ratio
bedtools bamtobed -i 07_blacklist/SIM_rep1.clean.bam > /tmp/SIM_r1_p1e3.bed
TOTAL_R1=$(wc -l < /tmp/SIM_r1_p1e3.bed)
IN_R1=$(bedtools intersect -u -a /tmp/SIM_r1_p1e3.bed \
  -b /tmp/SIM_idr_p1e3.bed | wc -l)
FRIP_R1_1e3=$(awk -v a="$IN_R1" -v b="$TOTAL_R1" 'BEGIN{printf "%.4f", a/b}')
rm /tmp/SIM_r1_p1e3.bed

bedtools bamtobed -i 07_blacklist/SIM_rep2.clean.bam > /tmp/SIM_r2_p1e3.bed
TOTAL_R2=$(wc -l < /tmp/SIM_r2_p1e3.bed)
IN_R2=$(bedtools intersect -u -a /tmp/SIM_r2_p1e3.bed \
  -b /tmp/SIM_idr_p1e3.bed | wc -l)
FRIP_R2_1e3=$(awk -v a="$IN_R2" -v b="$TOTAL_R2" 'BEGIN{printf "%.4f", a/b}')
rm /tmp/SIM_r2_p1e3.bed /tmp/SIM_idr_p1e3.bed

echo "1e-3 | $M1_1e3 | $M2_1e3 | $IDR_1e3 | $MU_1e3 | ${PCT_1e3}% | $FRIP_R1_1e3 | $FRIP_R2_1e3"

# ============================================================
# THRESHOLD 1e-2
# ============================================================

macs3 callpeak \
  -t 07_blacklist/SIM_rep1.clean.bam \
  -c 07_blacklist/SIM_control.clean.bam \
  -f BAMPE -g dm \
  -n SIM_rep1_p1e2 \
  --outdir 05_peaks_thresh/SIM \
  --keep-dup all \
  -p 1e-2 \
  --call-summits \
  2> 05_peaks_thresh/SIM/SIM_rep1_p1e2.log

macs3 callpeak \
  -t 07_blacklist/SIM_rep2.clean.bam \
  -c 07_blacklist/SIM_control.clean.bam \
  -f BAMPE -g dm \
  -n SIM_rep2_p1e2 \
  --outdir 05_peaks_thresh/SIM \
  --keep-dup all \
  -p 1e-2 \
  --call-summits \
  2> 05_peaks_thresh/SIM/SIM_rep2_p1e2.log

M1_1e2=$(wc -l < 05_peaks_thresh/SIM/SIM_rep1_p1e2_peaks.narrowPeak)
M2_1e2=$(wc -l < 05_peaks_thresh/SIM/SIM_rep2_p1e2_peaks.narrowPeak)

sort -k8,8nr 05_peaks_thresh/SIM/SIM_rep1_p1e2_peaks.narrowPeak \
  > 10_idr_thresh/SIM/SIM_rep1_p1e2.sorted.narrowPeak

sort -k8,8nr 05_peaks_thresh/SIM/SIM_rep2_p1e2_peaks.narrowPeak \
  > 10_idr_thresh/SIM/SIM_rep2_p1e2.sorted.narrowPeak
conda activate idr
idr \
  --samples \
    10_idr_thresh/SIM/SIM_rep1_p1e2.sorted.narrowPeak \
    10_idr_thresh/SIM/SIM_rep2_p1e2.sorted.narrowPeak \
  --input-file-type narrowPeak \
  --rank p.value \
  --idr-threshold 0.05 \
  --output-file 10_idr_thresh/SIM/SIM_1v2_p1e2.txt \
  --plot \
  --log-output-file 10_idr_thresh/SIM/SIM_1v2_p1e2.log
conda deactivate
NCOLS_1e2=$(awk '{print NF; exit}' 10_idr_thresh/SIM/SIM_1v2_p1e2.txt)
if [ "$NCOLS_1e2" -ge 12 ]; then
  IDR_1e2=$(awk '$12 >= 1.301' 10_idr_thresh/SIM/SIM_1v2_p1e2.txt | wc -l)
else
  IDR_1e2=$(awk '$5 >= 540'   10_idr_thresh/SIM/SIM_1v2_p1e2.txt | wc -l)
fi

MU_1e2=$(grep "Final parameter values" 10_idr_thresh/SIM/SIM_1v2_p1e2.log \
  | grep -oP '\[\K[0-9.]+' | head -1)
PCT_1e2=$(grep "Number of peaks passing" 10_idr_thresh/SIM/SIM_1v2_p1e2.log \
  | grep -oP '[0-9]+\.[0-9]+(?=%)' | head -1)

if [ "$NCOLS_1e2" -ge 12 ]; then
  awk '$12 >= 1.301 {print $1"\t"$2"\t"$3}' \
    10_idr_thresh/SIM/SIM_1v2_p1e2.txt \
    | sort -k1,1 -k2,2n | bedtools merge -i stdin \
    > /tmp/SIM_idr_p1e2.bed
else
  awk '$5 >= 540 {print $1"\t"$2"\t"$3}' \
    10_idr_thresh/SIM/SIM_1v2_p1e2.txt \
    | sort -k1,1 -k2,2n | bedtools merge -i stdin \
    > /tmp/SIM_idr_p1e2.bed
fi

bedtools bamtobed -i 07_blacklist/SIM_rep1.clean.bam > /tmp/SIM_r1_p1e2.bed
TOTAL_R1=$(wc -l < /tmp/SIM_r1_p1e2.bed)
IN_R1=$(bedtools intersect -u -a /tmp/SIM_r1_p1e2.bed \
  -b /tmp/SIM_idr_p1e2.bed | wc -l)
FRIP_R1_1e2=$(awk -v a="$IN_R1" -v b="$TOTAL_R1" 'BEGIN{printf "%.4f", a/b}')
rm /tmp/SIM_r1_p1e2.bed

bedtools bamtobed -i 07_blacklist/SIM_rep2.clean.bam > /tmp/SIM_r2_p1e2.bed
TOTAL_R2=$(wc -l < /tmp/SIM_r2_p1e2.bed)
IN_R2=$(bedtools intersect -u -a /tmp/SIM_r2_p1e2.bed \
  -b /tmp/SIM_idr_p1e2.bed | wc -l)
FRIP_R2_1e2=$(awk -v a="$IN_R2" -v b="$TOTAL_R2" 'BEGIN{printf "%.4f", a/b}')
rm /tmp/SIM_r2_p1e2.bed /tmp/SIM_idr_p1e2.bed

echo "1e-2 | $M1_1e2 | $M2_1e2 | $IDR_1e2 | $MU_1e2 | ${PCT_1e2}% | $FRIP_R1_1e2 | $FRIP_R2_1e2"


# ============================================================
# THRESHOLD 0.05
# ============================================================

conda activate phd
# --- MACS3 ---
macs3 callpeak \
  -t 07_blacklist/SIM_rep1.clean.bam \
  -c 07_blacklist/SIM_control.clean.bam \
  -f BAMPE -g dm \
  -n SIM_rep1_p5e2 \
  --outdir 05_peaks_thresh/SIM \
  --keep-dup all \
  -p 0.05 --call-summits \
  2> 05_peaks_thresh/SIM/SIM_rep1_p5e2.log

macs3 callpeak \
  -t 07_blacklist/SIM_rep2.clean.bam \
  -c 07_blacklist/SIM_control.clean.bam \
  -f BAMPE -g dm \
  -n SIM_rep2_p5e2 \
  --outdir 05_peaks_thresh/SIM \
  --keep-dup all \
  -p 0.05 --call-summits \
  2> 05_peaks_thresh/SIM/SIM_rep2_p5e2.log

SIM_M1=$(wc -l < 05_peaks_thresh/SIM/SIM_rep1_p5e2_peaks.narrowPeak)
SIM_M2=$(wc -l < 05_peaks_thresh/SIM/SIM_rep2_p5e2_peaks.narrowPeak)

# --- Sort ---
sort -k8,8nr 05_peaks_thresh/SIM/SIM_rep1_p5e2_peaks.narrowPeak \
  > 10_idr_thresh/SIM/SIM_rep1_p5e2.sorted.narrowPeak
sort -k8,8nr 05_peaks_thresh/SIM/SIM_rep2_p5e2_peaks.narrowPeak \
  > 10_idr_thresh/SIM/SIM_rep2_p5e2.sorted.narrowPeak
conda activate idr
# --- IDR ---
idr \
  --samples \
    10_idr_thresh/SIM/SIM_rep1_p5e2.sorted.narrowPeak \
    10_idr_thresh/SIM/SIM_rep2_p5e2.sorted.narrowPeak \
  --input-file-type narrowPeak \
  --rank p.value --idr-threshold 0.05 \
  --output-file 10_idr_thresh/SIM/SIM_1v2_p5e2.txt \
  --plot \
  --log-output-file 10_idr_thresh/SIM/SIM_1v2_p5e2.log
conda deactivate 
# --- Metrics ---
NCOLS=$(awk '{print NF; exit}' 10_idr_thresh/SIM/SIM_1v2_p5e2.txt)
if [ "$NCOLS" -ge 12 ]; then
  SIM_IDR=$(awk '$12 >= 1.301' 10_idr_thresh/SIM/SIM_1v2_p5e2.txt | wc -l)
  awk '$12 >= 1.301 {print $1"\t"$2"\t"$3}' \
    10_idr_thresh/SIM/SIM_1v2_p5e2.txt \
    | sort -k1,1 -k2,2n | bedtools merge -i stdin \
    > /tmp/SIM_idr_p5e2.bed
else
  SIM_IDR=$(awk '$5 >= 540' 10_idr_thresh/SIM/SIM_1v2_p5e2.txt | wc -l)
  awk '$5 >= 540 {print $1"\t"$2"\t"$3}' \
    10_idr_thresh/SIM/SIM_1v2_p5e2.txt \
    | sort -k1,1 -k2,2n | bedtools merge -i stdin \
    > /tmp/SIM_idr_p5e2.bed
fi

SIM_MU=$(grep "Final parameter values" 10_idr_thresh/SIM/SIM_1v2_p5e2.log \
  | grep -oP '\[\K[0-9.]+' | head -1)
SIM_PCT=$(grep "Number of peaks passing" 10_idr_thresh/SIM/SIM_1v2_p5e2.log \
  | grep -oP '[0-9]+\.[0-9]+(?=%)' | head -1)

# --- FRiP ---
bedtools bamtobed -i 07_blacklist/SIM_rep1.clean.bam > /tmp/SIM_r1_p5e2.bed
TOTAL=$(wc -l < /tmp/SIM_r1_p5e2.bed)
IN=$(bedtools intersect -u -a /tmp/SIM_r1_p5e2.bed \
  -b /tmp/SIM_idr_p5e2.bed | wc -l)
SIM_FRIP1=$(awk -v a="$IN" -v b="$TOTAL" 'BEGIN{printf "%.4f", a/b}')
rm /tmp/SIM_r1_p5e2.bed

bedtools bamtobed -i 07_blacklist/SIM_rep2.clean.bam > /tmp/SIM_r2_p5e2.bed
TOTAL=$(wc -l < /tmp/SIM_r2_p5e2.bed)
IN=$(bedtools intersect -u -a /tmp/SIM_r2_p5e2.bed \
  -b /tmp/SIM_idr_p5e2.bed | wc -l)
SIM_FRIP2=$(awk -v a="$IN" -v b="$TOTAL" 'BEGIN{printf "%.4f", a/b}')
rm /tmp/SIM_r2_p5e2.bed /tmp/SIM_idr_p5e2.bed

echo "THRESH | MACS_r1 | MACS_r2 | IDR_pass | IDR_mu | IDR_pct | FRIP_r1 | FRIP_r2"
echo "0.05   | $SIM_M1 | $SIM_M2 | $SIM_IDR | $SIM_MU | ${SIM_PCT}% | $SIM_FRIP1 | $SIM_FRIP2"
echo ""
echo "--- SIM context (from previous runs) ---"
echo "1e-3   | 1922    | 2923    | 493      | 1.54   | 43.5%   | 0.0084    | 0.0098"
echo "1e-2   | 6412    | 8080    | 587      | 1.67   | 22.7%   | 0.0095    | 0.0109"
echo "0.05   | (above)"
echo "1e-1   | 37109   | 37443   | 556      | 1.95   | 4.0%    | 0.0084    | 0.0092"



# ============================================================
# THRESHOLD 1e-1
# ============================================================

macs3 callpeak \
  -t 07_blacklist/SIM_rep1.clean.bam \
  -c 07_blacklist/SIM_control.clean.bam \
  -f BAMPE -g dm \
  -n SIM_rep1_p1e1 \
  --outdir 05_peaks_thresh/SIM \
  --keep-dup all \
  -p 1e-1 \
  --call-summits \
  2> 05_peaks_thresh/SIM/SIM_rep1_p1e1.log

macs3 callpeak \
  -t 07_blacklist/SIM_rep2.clean.bam \
  -c 07_blacklist/SIM_control.clean.bam \
  -f BAMPE -g dm \
  -n SIM_rep2_p1e1 \
  --outdir 05_peaks_thresh/SIM \
  --keep-dup all \
  -p 1e-1 \
  --call-summits \
  2> 05_peaks_thresh/SIM/SIM_rep2_p1e1.log

M1_1e1=$(wc -l < 05_peaks_thresh/SIM/SIM_rep1_p1e1_peaks.narrowPeak)
M2_1e1=$(wc -l < 05_peaks_thresh/SIM/SIM_rep2_p1e1_peaks.narrowPeak)

sort -k8,8nr 05_peaks_thresh/SIM/SIM_rep1_p1e1_peaks.narrowPeak \
  > 10_idr_thresh/SIM/SIM_rep1_p1e1.sorted.narrowPeak

sort -k8,8nr 05_peaks_thresh/SIM/SIM_rep2_p1e1_peaks.narrowPeak \
  > 10_idr_thresh/SIM/SIM_rep2_p1e1.sorted.narrowPeak
conda activate idr
idr \
  --samples \
    10_idr_thresh/SIM/SIM_rep1_p1e1.sorted.narrowPeak \
    10_idr_thresh/SIM/SIM_rep2_p1e1.sorted.narrowPeak \
  --input-file-type narrowPeak \
  --rank p.value \
  --idr-threshold 0.05 \
  --output-file 10_idr_thresh/SIM/SIM_1v2_p1e1.txt \
  --plot \
  --log-output-file 10_idr_thresh/SIM/SIM_1v2_p1e1.log
conda deactivate
NCOLS_1e1=$(awk '{print NF; exit}' 10_idr_thresh/SIM/SIM_1v2_p1e1.txt)
if [ "$NCOLS_1e1" -ge 12 ]; then
  IDR_1e1=$(awk '$12 >= 1.301' 10_idr_thresh/SIM/SIM_1v2_p1e1.txt | wc -l)
else
  IDR_1e1=$(awk '$5 >= 540'   10_idr_thresh/SIM/SIM_1v2_p1e1.txt | wc -l)
fi

MU_1e1=$(grep "Final parameter values" 10_idr_thresh/SIM/SIM_1v2_p1e1.log \
  | grep -oP '\[\K[0-9.]+' | head -1)
PCT_1e1=$(grep "Number of peaks passing" 10_idr_thresh/SIM/SIM_1v2_p1e1.log \
  | grep -oP '[0-9]+\.[0-9]+(?=%)' | head -1)

if [ "$NCOLS_1e1" -ge 12 ]; then
  awk '$12 >= 1.301 {print $1"\t"$2"\t"$3}' \
    10_idr_thresh/SIM/SIM_1v2_p1e1.txt \
    | sort -k1,1 -k2,2n | bedtools merge -i stdin \
    > /tmp/SIM_idr_p1e1.bed
else
  awk '$5 >= 540 {print $1"\t"$2"\t"$3}' \
    10_idr_thresh/SIM/SIM_1v2_p1e1.txt \
    | sort -k1,1 -k2,2n | bedtools merge -i stdin \
    > /tmp/SIM_idr_p1e1.bed
fi

bedtools bamtobed -i 07_blacklist/SIM_rep1.clean.bam > /tmp/SIM_r1_p1e1.bed
TOTAL_R1=$(wc -l < /tmp/SIM_r1_p1e1.bed)
IN_R1=$(bedtools intersect -u -a /tmp/SIM_r1_p1e1.bed \
  -b /tmp/SIM_idr_p1e1.bed | wc -l)
FRIP_R1_1e1=$(awk -v a="$IN_R1" -v b="$TOTAL_R1" 'BEGIN{printf "%.4f", a/b}')
rm /tmp/SIM_r1_p1e1.bed

bedtools bamtobed -i 07_blacklist/SIM_rep2.clean.bam > /tmp/SIM_r2_p1e1.bed
TOTAL_R2=$(wc -l < /tmp/SIM_r2_p1e1.bed)
IN_R2=$(bedtools intersect -u -a /tmp/SIM_r2_p1e1.bed \
  -b /tmp/SIM_idr_p1e1.bed | wc -l)
FRIP_R2_1e1=$(awk -v a="$IN_R2" -v b="$TOTAL_R2" 'BEGIN{printf "%.4f", a/b}')
rm /tmp/SIM_r2_p1e1.bed /tmp/SIM_idr_p1e1.bed

echo "1e-1 | $M1_1e1 | $M2_1e1 | $IDR_1e1 | $MU_1e1 | ${PCT_1e1}% | $FRIP_R1_1e1 | $FRIP_R2_1e1"

echo ""
echo "=== SIM threshold comparison complete ==="
echo "All MACS3 peaks: 05_peaks_thresh/SIM/"
echo "All IDR outputs: 10_idr_thresh/SIM/"
