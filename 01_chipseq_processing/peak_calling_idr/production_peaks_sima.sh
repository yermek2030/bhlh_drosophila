#!/usr/bin/env bash
# ============================================================
# SIMA PRODUCTION PEAK CALLING, IDR, CONSENSUS, SIGNAL NARROWPEAK
# Single-end; 2 biological replicates; p=0.01 (empirically selected)
# IDR strategy: single biological replicate IDR (rep1 vs rep2); n=2
# Output: 11_final/SIMA_IDR_consensus.bed
#         11_final/SIMA_IDR0.05.final.narrowPeak (for R/GRanges import)
# ============================================================

set -euo pipefail

mkdir -p 09_peaks/SIMA
mkdir -p 10_idr/SIMA
mkdir -p 11_final

# --- SPP-derived: SIMA_rep1 (230 bp) + SIMA_rep2 (220 bp) -> mean = 225 bp ---
# Matches TRH Chapter 3 approach (TRH used --extsize 210 from 205/215 bp SPP means)
EXTSIZE=225

# ============================================================
# STAGE 1, MACS3 PEAK CALLING (p=0.01, single-end)
# ============================================================

conda activate phd

for REP in 1 2
do
  macs3 callpeak \
    -t 07_blacklist/SIMA_rep${REP}.clean.bam \
    -c 07_blacklist/SIMA_control.clean.bam \
    -f BAM -g dm \
    -n SIMA_rep${REP} \
    --outdir 09_peaks/SIMA \
    --keep-dup all \
    --nomodel --extsize ${EXTSIZE} \
    -p 0.01 \
    --call-summits \
    2> 09_peaks/SIMA/SIMA_rep${REP}_macs3.log
done

echo "=== MACS3 peak counts ==="
for REP in 1 2
do
  echo -n "SIMA_rep${REP}: "
  wc -l < 09_peaks/SIMA/SIMA_rep${REP}_peaks.narrowPeak
done

# ============================================================
# STAGE 2, SORT BY SIGNAL SCORE FOR IDR
# ============================================================

for REP in 1 2
do
  sort -k8,8nr 09_peaks/SIMA/SIMA_rep${REP}_peaks.narrowPeak \
    > 10_idr/SIMA/SIMA_rep${REP}.sorted.narrowPeak
done

# ============================================================
# STAGE 3, IDR (rep1 vs rep2)
# ============================================================

conda activate idr

idr \
  --samples \
    10_idr/SIMA/SIMA_rep1.sorted.narrowPeak \
    10_idr/SIMA/SIMA_rep2.sorted.narrowPeak \
  --input-file-type narrowPeak \
  --rank p.value \
  --idr-threshold 0.05 \
  --output-file 10_idr/SIMA/SIMA_1v2.txt \
  --plot \
  --log-output-file 10_idr/SIMA/SIMA_1v2.log

conda deactivate

echo "=== IDR model summary ==="
grep "Final parameter values" 10_idr/SIMA/SIMA_1v2.log
grep "Number of peaks passing" 10_idr/SIMA/SIMA_1v2.log

# ============================================================
# STAGE 4, CONSENSUS PEAK SET (globalIDR <= 0.05)
# Column 12 >= 1.301 = globalIDR <= 0.05 in -log10 scale
# Falls back to col 5 >= 540 if col 12 absent
# ============================================================

NCOLS=$(awk '{print NF; exit}' 10_idr/SIMA/SIMA_1v2.txt)

if [ "$NCOLS" -ge 12 ]; then
  echo "Using column 12 (globalIDR -log10 scale)"
  awk '$12 >= 1.301 {print $1"\t"$2"\t"$3}' \
    10_idr/SIMA/SIMA_1v2.txt \
    | sort -k1,1 -k2,2n \
    | bedtools merge -i stdin \
    > 11_final/SIMA_IDR_consensus.bed
else
  echo "Column 12 absent — using column 5 (scaled score >= 540)"
  awk '$5 >= 540 {print $1"\t"$2"\t"$3}' \
    10_idr/SIMA/SIMA_1v2.txt \
    | sort -k1,1 -k2,2n \
    | bedtools merge -i stdin \
    > 11_final/SIMA_IDR_consensus.bed
fi

echo -n "SIMA consensus peaks (IDR <= 0.05): "
wc -l < 11_final/SIMA_IDR_consensus.bed

# ============================================================
# STAGE 5, MEAN SIGNAL FROM MERGED BAM
# ============================================================

bedtools coverage \
  -a 11_final/SIMA_IDR_consensus.bed \
  -b 07_blacklist/SIMA_merged.clean.bam \
  -mean \
  > 11_final/SIMA_IDR_with_signal.txt

echo "Head of signal file:"
head -3 11_final/SIMA_IDR_with_signal.txt

# ============================================================
# STAGE 6, FORMAT AS NARROWPEAK
# cols: chr start end name score strand signal pval qval summit
# summit = midpoint offset from peak start
# ============================================================

awk 'BEGIN{OFS="\t"}{
  print $1, $2, $3, ".", 1000, ".", $4, -1, -1, int(($3-$2)/2)
}' 11_final/SIMA_IDR_with_signal.txt \
  > 11_final/SIMA_IDR0.05.clean.narrowPeak

# ============================================================
# STAGE 7, CANONICAL CHROMOSOME FILTER
# Retain chr2L chr2R chr3L chr3R chrX chr4 only
# ============================================================

grep -E '^chr[234X][LR]?[[:space:]]' \
  11_final/SIMA_IDR0.05.clean.narrowPeak \
  > 11_final/SIMA_IDR0.05.final.narrowPeak

echo "=== Final peak counts ==="
echo -n "Before chr filter: "
wc -l < 11_final/SIMA_IDR0.05.clean.narrowPeak
echo -n "After chr filter:  "
wc -l < 11_final/SIMA_IDR0.05.final.narrowPeak

# ============================================================
# STAGE 8, FRIP ON FINAL PEAK SET
# ============================================================

for REP in 1 2
do
  echo -n "SIMA_rep${REP} FRiP: "
  bedtools bamtobed \
    -i 07_blacklist/SIMA_rep${REP}.clean.bam \
    > /tmp/SIMA_frip_rep${REP}.bed
  total=$(wc -l < /tmp/SIMA_frip_rep${REP}.bed)
  inpeaks=$(bedtools intersect -u \
    -a /tmp/SIMA_frip_rep${REP}.bed \
    -b 11_final/SIMA_IDR0.05.final.narrowPeak | wc -l)
  awk -v a="$inpeaks" -v b="$total" \
    'BEGIN { printf "%.4f\n", a/b }'
  rm /tmp/SIMA_frip_rep${REP}.bed
done

# ============================================================
# STAGE 9, QC CHECKS
# ============================================================

echo ""
echo "=== Signal range check ==="
awk '{print $7}' 11_final/SIMA_IDR0.05.final.narrowPeak \
  | sort -n \
  | awk 'NR==1{min=$1} END{print "min: "min"  max: "$1"  peaks: "NR}'

echo ""
echo "=== Final narrowPeak structure (first 3 lines) ==="
head -3 11_final/SIMA_IDR0.05.final.narrowPeak

echo ""
echo "=== PIPELINE COMPLETE ==="
echo "Consensus BED:      11_final/SIMA_IDR_consensus.bed"
echo "Final narrowPeak:   11_final/SIMA_IDR0.05.final.narrowPeak"
wc -l 11_final/SIMA_IDR0.05.final.narrowPeak

# ============================================================
# NEXT STEP: run pseudoreplicate_idr_sima.sh for pseudoreplicate IDR validation
# This is particularly important for n=2 single-end datasets (same rationale
# as TRH Chapter 3) to check against threshold-dependent rank instability.
# ============================================================