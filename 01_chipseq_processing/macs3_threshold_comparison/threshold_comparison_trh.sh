#!/usr/bin/env bash
# ============================================================
# TRH MACS3 THRESHOLD COMPARISON
# Single-end library; 2 biological replicates + 1 matched input control
# Thresholds tested: 1e-3, 1e-2, 1e-1
# Output: one summary row per threshold printed to stdout
# Run: bash threshold_comparison_trh.sh 2> threshold_comparison_trh.err
# ============================================================

set -euo pipefail
mkdir -p 05_peaks_thresh/TRH
mkdir -p 10_idr_thresh/TRH

echo "THRESH | MACS_rep1 | MACS_rep2 | IDR_pass | IDR_mu | IDR_pct | FRIP_rep1 | FRIP_rep2"

# ============================================================
# THRESHOLD 1e-3
# ============================================================
conda activate phd
# --- MACS3 ---
macs3 callpeak \
  -t 07_blacklist/TRH_rep1.clean.bam \
  -c 07_blacklist/TRH_control.clean.bam \
  -f BAM -g dm \
  -n TRH_rep1_p1e3 \
  --outdir 05_peaks_thresh/TRH \
  --keep-dup all \
  --nomodel --extsize 210 \
  -p 1e-3 \
  --call-summits \
  2> 05_peaks_thresh/TRH/TRH_rep1_p1e3.log

macs3 callpeak \
  -t 07_blacklist/TRH_rep2.clean.bam \
  -c 07_blacklist/TRH_control.clean.bam \
  -f BAM -g dm \
  -n TRH_rep2_p1e3 \
  --outdir 05_peaks_thresh/TRH \
  --keep-dup all \
  --nomodel --extsize 210 \
  -p 1e-3 \
  --call-summits \
  2> 05_peaks_thresh/TRH/TRH_rep2_p1e3.log

M1_1e3=$(wc -l < 05_peaks_thresh/TRH/TRH_rep1_p1e3_peaks.narrowPeak)
M2_1e3=$(wc -l < 05_peaks_thresh/TRH/TRH_rep2_p1e3_peaks.narrowPeak)

# --- Sort by signal score (col 8) ---
sort -k8,8nr 05_peaks_thresh/TRH/TRH_rep1_p1e3_peaks.narrowPeak \
  > 10_idr_thresh/TRH/TRH_rep1_p1e3.sorted.narrowPeak

sort -k8,8nr 05_peaks_thresh/TRH/TRH_rep2_p1e3_peaks.narrowPeak \
  > 10_idr_thresh/TRH/TRH_rep2_p1e3.sorted.narrowPeak
conda activate idr
# --- IDR ---
idr \
  --samples \
    10_idr_thresh/TRH/TRH_rep1_p1e3.sorted.narrowPeak \
    10_idr_thresh/TRH/TRH_rep2_p1e3.sorted.narrowPeak \
  --input-file-type narrowPeak \
  --rank p.value \
  --idr-threshold 0.05 \
  --output-file 10_idr_thresh/TRH/TRH_1v2_p1e3.txt \
  --plot \
  --log-output-file 10_idr_thresh/TRH/TRH_1v2_p1e3.log
conda deactivate
# --- Detect column count for filter ---
NCOLS_1e3=$(awk '{print NF; exit}' 10_idr_thresh/TRH/TRH_1v2_p1e3.txt)
if [ "$NCOLS_1e3" -ge 12 ]; then
  IDR_1e3=$(awk '$12 >= 1.301' 10_idr_thresh/TRH/TRH_1v2_p1e3.txt | wc -l)
else
  IDR_1e3=$(awk '$5 >= 540'   10_idr_thresh/TRH/TRH_1v2_p1e3.txt | wc -l)
fi

MU_1e3=$(grep "Final parameter values" 10_idr_thresh/TRH/TRH_1v2_p1e3.log \
  | grep -oP '\[\K[0-9.]+' | head -1)
PCT_1e3=$(grep "Number of peaks passing" 10_idr_thresh/TRH/TRH_1v2_p1e3.log \
  | grep -oP '[0-9]+\.[0-9]+(?=%)' | head -1)

# --- Build IDR peak BED for FRiP ---
if [ "$NCOLS_1e3" -ge 12 ]; then
  awk '$12 >= 1.301 {print $1"\t"$2"\t"$3}' \
    10_idr_thresh/TRH/TRH_1v2_p1e3.txt \
    | sort -k1,1 -k2,2n | bedtools merge -i stdin \
    > /tmp/TRH_idr_p1e3.bed
else
  awk '$5 >= 540 {print $1"\t"$2"\t"$3}' \
    10_idr_thresh/TRH/TRH_1v2_p1e3.txt \
    | sort -k1,1 -k2,2n | bedtools merge -i stdin \
    > /tmp/TRH_idr_p1e3.bed
fi

# --- FRiP rep1 ---
bedtools bamtobed -i 07_blacklist/TRH_rep1.clean.bam > /tmp/TRH_r1_p1e3.bed
TOTAL_R1=$(wc -l < /tmp/TRH_r1_p1e3.bed)
IN_R1=$(bedtools intersect -u -a /tmp/TRH_r1_p1e3.bed \
  -b /tmp/TRH_idr_p1e3.bed | wc -l)
FRIP_R1_1e3=$(awk -v a="$IN_R1" -v b="$TOTAL_R1" 'BEGIN{printf "%.4f", a/b}')
rm /tmp/TRH_r1_p1e3.bed

# --- FRiP rep2 ---
bedtools bamtobed -i 07_blacklist/TRH_rep2.clean.bam > /tmp/TRH_r2_p1e3.bed
TOTAL_R2=$(wc -l < /tmp/TRH_r2_p1e3.bed)
IN_R2=$(bedtools intersect -u -a /tmp/TRH_r2_p1e3.bed \
  -b /tmp/TRH_idr_p1e3.bed | wc -l)
FRIP_R2_1e3=$(awk -v a="$IN_R2" -v b="$TOTAL_R2" 'BEGIN{printf "%.4f", a/b}')
rm /tmp/TRH_r2_p1e3.bed /tmp/TRH_idr_p1e3.bed

echo "1e-3 | $M1_1e3 | $M2_1e3 | $IDR_1e3 | $MU_1e3 | ${PCT_1e3}% | $FRIP_R1_1e3 | $FRIP_R2_1e3"

# ============================================================
# THRESHOLD 1e-2
# ============================================================

macs3 callpeak \
  -t 07_blacklist/TRH_rep1.clean.bam \
  -c 07_blacklist/TRH_control.clean.bam \
  -f BAM -g dm \
  -n TRH_rep1_p1e2 \
  --outdir 05_peaks_thresh/TRH \
  --keep-dup all \
  --nomodel --extsize 210 \
  -p 1e-2 \
  --call-summits \
  2> 05_peaks_thresh/TRH/TRH_rep1_p1e2.log

macs3 callpeak \
  -t 07_blacklist/TRH_rep2.clean.bam \
  -c 07_blacklist/TRH_control.clean.bam \
  -f BAM -g dm \
  -n TRH_rep2_p1e2 \
  --outdir 05_peaks_thresh/TRH \
  --keep-dup all \
  --nomodel --extsize 210 \
  -p 1e-2 \
  --call-summits \
  2> 05_peaks_thresh/TRH/TRH_rep2_p1e2.log

M1_1e2=$(wc -l < 05_peaks_thresh/TRH/TRH_rep1_p1e2_peaks.narrowPeak)
M2_1e2=$(wc -l < 05_peaks_thresh/TRH/TRH_rep2_p1e2_peaks.narrowPeak)

sort -k8,8nr 05_peaks_thresh/TRH/TRH_rep1_p1e2_peaks.narrowPeak \
  > 10_idr_thresh/TRH/TRH_rep1_p1e2.sorted.narrowPeak

sort -k8,8nr 05_peaks_thresh/TRH/TRH_rep2_p1e2_peaks.narrowPeak \
  > 10_idr_thresh/TRH/TRH_rep2_p1e2.sorted.narrowPeak
conda activate idr
idr \
  --samples \
    10_idr_thresh/TRH/TRH_rep1_p1e2.sorted.narrowPeak \
    10_idr_thresh/TRH/TRH_rep2_p1e2.sorted.narrowPeak \
  --input-file-type narrowPeak \
  --rank p.value \
  --idr-threshold 0.05 \
  --output-file 10_idr_thresh/TRH/TRH_1v2_p1e2.txt \
  --plot \
  --log-output-file 10_idr_thresh/TRH/TRH_1v2_p1e2.log
conda deactivate 
NCOLS_1e2=$(awk '{print NF; exit}' 10_idr_thresh/TRH/TRH_1v2_p1e2.txt)
if [ "$NCOLS_1e2" -ge 12 ]; then
  IDR_1e2=$(awk '$12 >= 1.301' 10_idr_thresh/TRH/TRH_1v2_p1e2.txt | wc -l)
else
  IDR_1e2=$(awk '$5 >= 540'   10_idr_thresh/TRH/TRH_1v2_p1e2.txt | wc -l)
fi

MU_1e2=$(grep "Final parameter values" 10_idr_thresh/TRH/TRH_1v2_p1e2.log \
  | grep -oP '\[\K[0-9.]+' | head -1)
PCT_1e2=$(grep "Number of peaks passing" 10_idr_thresh/TRH/TRH_1v2_p1e2.log \
  | grep -oP '[0-9]+\.[0-9]+(?=%)' | head -1)

if [ "$NCOLS_1e2" -ge 12 ]; then
  awk '$12 >= 1.301 {print $1"\t"$2"\t"$3}' \
    10_idr_thresh/TRH/TRH_1v2_p1e2.txt \
    | sort -k1,1 -k2,2n | bedtools merge -i stdin \
    > /tmp/TRH_idr_p1e2.bed
else
  awk '$5 >= 540 {print $1"\t"$2"\t"$3}' \
    10_idr_thresh/TRH/TRH_1v2_p1e2.txt \
    | sort -k1,1 -k2,2n | bedtools merge -i stdin \
    > /tmp/TRH_idr_p1e2.bed
fi

bedtools bamtobed -i 07_blacklist/TRH_rep1.clean.bam > /tmp/TRH_r1_p1e2.bed
TOTAL_R1=$(wc -l < /tmp/TRH_r1_p1e2.bed)
IN_R1=$(bedtools intersect -u -a /tmp/TRH_r1_p1e2.bed \
  -b /tmp/TRH_idr_p1e2.bed | wc -l)
FRIP_R1_1e2=$(awk -v a="$IN_R1" -v b="$TOTAL_R1" 'BEGIN{printf "%.4f", a/b}')
rm /tmp/TRH_r1_p1e2.bed

bedtools bamtobed -i 07_blacklist/TRH_rep2.clean.bam > /tmp/TRH_r2_p1e2.bed
TOTAL_R2=$(wc -l < /tmp/TRH_r2_p1e2.bed)
IN_R2=$(bedtools intersect -u -a /tmp/TRH_r2_p1e2.bed \
  -b /tmp/TRH_idr_p1e2.bed | wc -l)
FRIP_R2_1e2=$(awk -v a="$IN_R2" -v b="$TOTAL_R2" 'BEGIN{printf "%.4f", a/b}')
rm /tmp/TRH_r2_p1e2.bed /tmp/TRH_idr_p1e2.bed

echo "1e-2 | $M1_1e2 | $M2_1e2 | $IDR_1e2 | $MU_1e2 | ${PCT_1e2}% | $FRIP_R1_1e2 | $FRIP_R2_1e2"

# ============================================================
# TRH p=0.05 SINGLE THRESHOLD TEST
# Single-end; 2 replicates; extsize 210
# ============================================================

conda activate phd
macs3 callpeak \
  -t 07_blacklist/TRH_rep1.clean.bam \
  -c 07_blacklist/TRH_control.clean.bam \
  -f BAM -g dm \
  -n TRH_rep1_p5e2 \
  --outdir 05_peaks_thresh/TRH \
  --keep-dup all \
  --nomodel --extsize 210 \
  -p 0.05 --call-summits \
  2> 05_peaks_thresh/TRH/TRH_rep1_p5e2.log

macs3 callpeak \
  -t 07_blacklist/TRH_rep2.clean.bam \
  -c 07_blacklist/TRH_control.clean.bam \
  -f BAM -g dm \
  -n TRH_rep2_p5e2 \
  --outdir 05_peaks_thresh/TRH \
  --keep-dup all \
  --nomodel --extsize 210 \
  -p 0.05 --call-summits \
  2> 05_peaks_thresh/TRH/TRH_rep2_p5e2.log

TRH_M1=$(wc -l < 05_peaks_thresh/TRH/TRH_rep1_p5e2_peaks.narrowPeak)
TRH_M2=$(wc -l < 05_peaks_thresh/TRH/TRH_rep2_p5e2_peaks.narrowPeak)

sort -k8,8nr 05_peaks_thresh/TRH/TRH_rep1_p5e2_peaks.narrowPeak \
  > 10_idr_thresh/TRH/TRH_rep1_p5e2.sorted.narrowPeak
sort -k8,8nr 05_peaks_thresh/TRH/TRH_rep2_p5e2_peaks.narrowPeak \
  > 10_idr_thresh/TRH/TRH_rep2_p5e2.sorted.narrowPeak
conda activate idr
idr \
  --samples \
    10_idr_thresh/TRH/TRH_rep1_p5e2.sorted.narrowPeak \
    10_idr_thresh/TRH/TRH_rep2_p5e2.sorted.narrowPeak \
  --input-file-type narrowPeak \
  --rank p.value --idr-threshold 0.05 \
  --output-file 10_idr_thresh/TRH/TRH_1v2_p5e2.txt \
  --plot \
  --log-output-file 10_idr_thresh/TRH/TRH_1v2_p5e2.log
conda deactivate
NCOLS=$(awk '{print NF; exit}' 10_idr_thresh/TRH/TRH_1v2_p5e2.txt)
if [ "$NCOLS" -ge 12 ]; then
  TRH_IDR=$(awk '$12 >= 1.301' 10_idr_thresh/TRH/TRH_1v2_p5e2.txt | wc -l)
  awk '$12 >= 1.301 {print $1"\t"$2"\t"$3}' \
    10_idr_thresh/TRH/TRH_1v2_p5e2.txt \
    | sort -k1,1 -k2,2n | bedtools merge -i stdin \
    > /tmp/TRH_idr_p5e2.bed
else
  TRH_IDR=$(awk '$5 >= 540' 10_idr_thresh/TRH/TRH_1v2_p5e2.txt | wc -l)
  awk '$5 >= 540 {print $1"\t"$2"\t"$3}' \
    10_idr_thresh/TRH/TRH_1v2_p5e2.txt \
    | sort -k1,1 -k2,2n | bedtools merge -i stdin \
    > /tmp/TRH_idr_p5e2.bed
fi

TRH_MU=$(grep "Final parameter values" 10_idr_thresh/TRH/TRH_1v2_p5e2.log \
  | grep -oP '\[\K[0-9.]+' | head -1)
TRH_PCT=$(grep "Number of peaks passing" 10_idr_thresh/TRH/TRH_1v2_p5e2.log \
  | grep -oP '[0-9]+\.[0-9]+(?=%)' | head -1)

bedtools bamtobed -i 07_blacklist/TRH_rep1.clean.bam > /tmp/TRH_r1_p5e2.bed
TOTAL=$(wc -l < /tmp/TRH_r1_p5e2.bed)
IN=$(bedtools intersect -u -a /tmp/TRH_r1_p5e2.bed \
  -b /tmp/TRH_idr_p5e2.bed | wc -l)
TRH_FRIP1=$(awk -v a="$IN" -v b="$TOTAL" 'BEGIN{printf "%.4f", a/b}')
rm /tmp/TRH_r1_p5e2.bed

bedtools bamtobed -i 07_blacklist/TRH_rep2.clean.bam > /tmp/TRH_r2_p5e2.bed
TOTAL=$(wc -l < /tmp/TRH_r2_p5e2.bed)
IN=$(bedtools intersect -u -a /tmp/TRH_r2_p5e2.bed \
  -b /tmp/TRH_idr_p5e2.bed | wc -l)
TRH_FRIP2=$(awk -v a="$IN" -v b="$TOTAL" 'BEGIN{printf "%.4f", a/b}')
rm /tmp/TRH_r2_p5e2.bed /tmp/TRH_idr_p5e2.bed

echo "THRESH | MACS_r1 | MACS_r2 | IDR_pass | IDR_mu | IDR_pct | FRIP_r1 | FRIP_r2"
echo "0.05   | $TRH_M1 | $TRH_M2 | $TRH_IDR | $TRH_MU | ${TRH_PCT}% | $TRH_FRIP1 | $TRH_FRIP2"
echo ""
echo "--- TRH context (from previous runs) ---"
echo "1e-3 | 3709  | 1899  | 707  | 1.63 | 43.1% | 0.0115 | 0.0097"
echo "1e-2 | 6218  | 3560  | 1033 | 1.64 | 39.3% | 0.0160 | 0.0135"
echo "0.05 | (above)"
echo "1e-1 | 17864 | 17460 | 1511 | 1.68 | 23.1% | 0.0222 | 0.0186"
echo ""
echo "=== TRH p=0.05 test complete ==="


# ============================================================
# THRESHOLD 1e-1
# ============================================================

macs3 callpeak \
  -t 07_blacklist/TRH_rep1.clean.bam \
  -c 07_blacklist/TRH_control.clean.bam \
  -f BAM -g dm \
  -n TRH_rep1_p1e1 \
  --outdir 05_peaks_thresh/TRH \
  --keep-dup all \
  --nomodel --extsize 210 \
  -p 1e-1 \
  --call-summits \
  2> 05_peaks_thresh/TRH/TRH_rep1_p1e1.log

macs3 callpeak \
  -t 07_blacklist/TRH_rep2.clean.bam \
  -c 07_blacklist/TRH_control.clean.bam \
  -f BAM -g dm \
  -n TRH_rep2_p1e1 \
  --outdir 05_peaks_thresh/TRH \
  --keep-dup all \
  --nomodel --extsize 210 \
  -p 1e-1 \
  --call-summits \
  2> 05_peaks_thresh/TRH/TRH_rep2_p1e1.log

M1_1e1=$(wc -l < 05_peaks_thresh/TRH/TRH_rep1_p1e1_peaks.narrowPeak)
M2_1e1=$(wc -l < 05_peaks_thresh/TRH/TRH_rep2_p1e1_peaks.narrowPeak)

sort -k8,8nr 05_peaks_thresh/TRH/TRH_rep1_p1e1_peaks.narrowPeak \
  > 10_idr_thresh/TRH/TRH_rep1_p1e1.sorted.narrowPeak

sort -k8,8nr 05_peaks_thresh/TRH/TRH_rep2_p1e1_peaks.narrowPeak \
  > 10_idr_thresh/TRH/TRH_rep2_p1e1.sorted.narrowPeak
conda activate idr
idr \
  --samples \
    10_idr_thresh/TRH/TRH_rep1_p1e1.sorted.narrowPeak \
    10_idr_thresh/TRH/TRH_rep2_p1e1.sorted.narrowPeak \
  --input-file-type narrowPeak \
  --rank p.value \
  --idr-threshold 0.05 \
  --output-file 10_idr_thresh/TRH/TRH_1v2_p1e1.txt \
  --plot \
  --log-output-file 10_idr_thresh/TRH/TRH_1v2_p1e1.log
conda deactivate 
NCOLS_1e1=$(awk '{print NF; exit}' 10_idr_thresh/TRH/TRH_1v2_p1e1.txt)
if [ "$NCOLS_1e1" -ge 12 ]; then
  IDR_1e1=$(awk '$12 >= 1.301' 10_idr_thresh/TRH/TRH_1v2_p1e1.txt | wc -l)
else
  IDR_1e1=$(awk '$5 >= 540'   10_idr_thresh/TRH/TRH_1v2_p1e1.txt | wc -l)
fi

MU_1e1=$(grep "Final parameter values" 10_idr_thresh/TRH/TRH_1v2_p1e1.log \
  | grep -oP '\[\K[0-9.]+' | head -1)
PCT_1e1=$(grep "Number of peaks passing" 10_idr_thresh/TRH/TRH_1v2_p1e1.log \
  | grep -oP '[0-9]+\.[0-9]+(?=%)' | head -1)

if [ "$NCOLS_1e1" -ge 12 ]; then
  awk '$12 >= 1.301 {print $1"\t"$2"\t"$3}' \
    10_idr_thresh/TRH/TRH_1v2_p1e1.txt \
    | sort -k1,1 -k2,2n | bedtools merge -i stdin \
    > /tmp/TRH_idr_p1e1.bed
else
  awk '$5 >= 540 {print $1"\t"$2"\t"$3}' \
    10_idr_thresh/TRH/TRH_1v2_p1e1.txt \
    | sort -k1,1 -k2,2n | bedtools merge -i stdin \
    > /tmp/TRH_idr_p1e1.bed
fi

bedtools bamtobed -i 07_blacklist/TRH_rep1.clean.bam > /tmp/TRH_r1_p1e1.bed
TOTAL_R1=$(wc -l < /tmp/TRH_r1_p1e1.bed)
IN_R1=$(bedtools intersect -u -a /tmp/TRH_r1_p1e1.bed \
  -b /tmp/TRH_idr_p1e1.bed | wc -l)
FRIP_R1_1e1=$(awk -v a="$IN_R1" -v b="$TOTAL_R1" 'BEGIN{printf "%.4f", a/b}')
rm /tmp/TRH_r1_p1e1.bed

bedtools bamtobed -i 07_blacklist/TRH_rep2.clean.bam > /tmp/TRH_r2_p1e1.bed
TOTAL_R2=$(wc -l < /tmp/TRH_r2_p1e1.bed)
IN_R2=$(bedtools intersect -u -a /tmp/TRH_r2_p1e1.bed \
  -b /tmp/TRH_idr_p1e1.bed | wc -l)
FRIP_R2_1e1=$(awk -v a="$IN_R2" -v b="$TOTAL_R2" 'BEGIN{printf "%.4f", a/b}')
rm /tmp/TRH_r2_p1e1.bed /tmp/TRH_idr_p1e1.bed

echo "1e-1 | $M1_1e1 | $M2_1e1 | $IDR_1e1 | $MU_1e1 | ${PCT_1e1}% | $FRIP_R1_1e1 | $FRIP_R2_1e1"

echo ""
echo "=== TRH threshold comparison complete ==="
echo "All MACS3 peaks: 05_peaks_thresh/TRH/"
echo "All IDR outputs: 10_idr_thresh/TRH/"
