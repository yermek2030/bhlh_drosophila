#!/usr/bin/env bash
# ============================================================
# TRH PRODUCTION PEAK CALLING, IDR, CONSENSUS, SIGNAL NARROWPEAK
# Single-end; 2 biological replicates; p=0.01 (empirically selected)
# Output: 11_final/TRH_IDR_consensus.bed
#         11_final/TRH_IDR0.05.final.narrowPeak (for R/GRanges import)
# ============================================================

set -euo pipefail

mkdir -p 09_peaks/TRH
mkdir -p 10_idr/TRH
mkdir -p 11_final

# ============================================================
# STAGE 1, MACS3 PEAK CALLING (p=0.01, single-end)
# extsize 210 = mean of SPP estimates (rep1: 205 bp, rep2: 215 bp)
# ============================================================

macs3 callpeak \
  -t 07_blacklist/TRH_rep1.clean.bam \
  -c 07_blacklist/TRH_control.clean.bam \
  -f BAM -g dm \
  -n TRH_rep1 \
  --outdir 09_peaks/TRH \
  --keep-dup all \
  --nomodel --extsize 210 \
  -p 0.01 \
  --call-summits \
  2> 09_peaks/TRH/TRH_rep1_macs3.log

macs3 callpeak \
  -t 07_blacklist/TRH_rep2.clean.bam \
  -c 07_blacklist/TRH_control.clean.bam \
  -f BAM -g dm \
  -n TRH_rep2 \
  --outdir 09_peaks/TRH \
  --keep-dup all \
  --nomodel --extsize 210 \
  -p 0.01 \
  --call-summits \
  2> 09_peaks/TRH/TRH_rep2_macs3.log

echo "=== MACS3 peak counts ==="
echo -n "TRH_rep1: "
wc -l < 09_peaks/TRH/TRH_rep1_peaks.narrowPeak
echo -n "TRH_rep2: "
wc -l < 09_peaks/TRH/TRH_rep2_peaks.narrowPeak

# ============================================================
# STAGE 2, SORT BY SIGNAL SCORE FOR IDR
# ============================================================

sort -k8,8nr 09_peaks/TRH/TRH_rep1_peaks.narrowPeak \
  > 10_idr/TRH/TRH_rep1.sorted.narrowPeak

sort -k8,8nr 09_peaks/TRH/TRH_rep2_peaks.narrowPeak \
  > 10_idr/TRH/TRH_rep2.sorted.narrowPeak

# ============================================================
# STAGE 3, IDR (rep1 vs rep2)
# ============================================================

conda activate idr

idr \
  --samples \
    10_idr/TRH/TRH_rep1.sorted.narrowPeak \
    10_idr/TRH/TRH_rep2.sorted.narrowPeak \
  --input-file-type narrowPeak \
  --rank p.value \
  --idr-threshold 0.05 \
  --output-file 10_idr/TRH/TRH_1v2.txt \
  --plot \
  --log-output-file 10_idr/TRH/TRH_1v2.log

conda deactivate

echo "=== IDR model summary ==="
grep "Final parameter values" 10_idr/TRH/TRH_1v2.log
grep "Number of peaks passing" 10_idr/TRH/TRH_1v2.log

# ============================================================
# STAGE 4, CONSENSUS PEAK SET (globalIDR <= 0.05)
# Column 12 >= 1.301 = globalIDR <= 0.05 in -log10 scale
# Falls back to col 5 >= 540 if col 12 absent
# ============================================================

NCOLS=$(awk '{print NF; exit}' 10_idr/TRH/TRH_1v2.txt)

if [ "$NCOLS" -ge 12 ]; then
  echo "Using column 12 (globalIDR -log10 scale)"
  awk '$12 >= 1.301 {print $1"\t"$2"\t"$3}' \
    10_idr/TRH/TRH_1v2.txt \
    | sort -k1,1 -k2,2n \
    | bedtools merge -i stdin \
    > 11_final/TRH_IDR_consensus.bed
else
  echo "Column 12 absent — using column 5 (scaled score >= 540)"
  awk '$5 >= 540 {print $1"\t"$2"\t"$3}' \
    10_idr/TRH/TRH_1v2.txt \
    | sort -k1,1 -k2,2n \
    | bedtools merge -i stdin \
    > 11_final/TRH_IDR_consensus.bed
fi

echo -n "TRH consensus peaks (IDR <= 0.05): "
wc -l < 11_final/TRH_IDR_consensus.bed

# ============================================================
# STAGE 5, MEAN SIGNAL FROM MERGED BAM
# ============================================================

bedtools coverage \
  -a 11_final/TRH_IDR_consensus.bed \
  -b 07_blacklist/TRH_merged.clean.bam \
  -mean \
  > 11_final/TRH_IDR_with_signal.txt

echo "Head of signal file:"
head -3 11_final/TRH_IDR_with_signal.txt

# ============================================================
# STAGE 6, FORMAT AS NARROWPEAK
# cols: chr start end name score strand signal pval qval summit
# summit = midpoint offset from peak start
# ============================================================

awk 'BEGIN{OFS="\t"}{
  print $1, $2, $3, ".", 1000, ".", $4, -1, -1, int(($3-$2)/2)
}' 11_final/TRH_IDR_with_signal.txt \
  > 11_final/TRH_IDR0.05.clean.narrowPeak

# ============================================================
# STAGE 7, CANONICAL CHROMOSOME FILTER
# Retain chr2L chr2R chr3L chr3R chrX only
# ============================================================

grep -E '^chr[234X][LR]?[[:space:]]' \
  11_final/TRH_IDR0.05.clean.narrowPeak \
  > 11_final/TRH_IDR0.05.final.narrowPeak

echo "=== Final peak counts ==="
echo -n "Before chr filter: "
wc -l < 11_final/TRH_IDR0.05.clean.narrowPeak
echo -n "After chr filter:  "
wc -l < 11_final/TRH_IDR0.05.final.narrowPeak

# ============================================================
# STAGE 8, FRIP ON FINAL PEAK SET
# ============================================================

for REP in 1 2
do
  echo -n "TRH_rep${REP} FRiP: "
  bedtools bamtobed \
    -i 07_blacklist/TRH_rep${REP}.clean.bam \
    > /tmp/TRH_frip_rep${REP}.bed
  total=$(wc -l < /tmp/TRH_frip_rep${REP}.bed)
  inpeaks=$(bedtools intersect -u \
    -a /tmp/TRH_frip_rep${REP}.bed \
    -b 11_final/TRH_IDR0.05.final.narrowPeak | wc -l)
  awk -v a="$inpeaks" -v b="$total" \
    'BEGIN { printf "%.4f\n", a/b }'
  rm /tmp/TRH_frip_rep${REP}.bed
done

# ============================================================
# STAGE 9, QC CHECKS
# ============================================================

echo ""
echo "=== Signal range check ==="
awk '{print $7}' 11_final/TRH_IDR0.05.final.narrowPeak \
  | sort -n \
  | awk 'NR==1{min=$1} END{print "min: "min"  max: "$1"  peaks: "NR}'

echo ""
echo "=== Final narrowPeak structure (first 3 lines) ==="
head -3 11_final/TRH_IDR0.05.final.narrowPeak

echo ""
echo "=== PIPELINE COMPLETE ==="
echo "Consensus BED:      11_final/TRH_IDR_consensus.bed"
echo "Final narrowPeak:   11_final/TRH_IDR0.05.final.narrowPeak"
wc -l 11_final/TRH_IDR0.05.final.narrowPeak
