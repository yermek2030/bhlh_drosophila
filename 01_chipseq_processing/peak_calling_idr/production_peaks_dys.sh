#!/usr/bin/env bash
# ============================================================
# DYS peak calling: paired-end, 3 biological replicates, p=0.01. MACS3 -f BAMPE, no --nomodel.
# IDR pairwise (1v2, 1v3, 2v3); consensus = peaks in >=2/3 pairs.
# Output: 11_final/DYS_IDR_consensus.bed
#         11_final/DYS_IDR0.05.final.narrowPeak
# ============================================================

set -euo pipefail

mkdir -p 09_peaks/DYS
mkdir -p 10_idr/DYS
mkdir -p 11_final

# ============================================================
# STAGE 1, MACS3 PEAK CALLING (p=0.01, paired-end BAMPE)
# ============================================================

conda activate phd

for REP in 1 2 3
do
  macs3 callpeak \
    -t 07_blacklist/DYS_rep${REP}.clean.bam \
    -c 07_blacklist/DYS_control.clean.bam \
    -f BAMPE -g dm \
    -n DYS_rep${REP} \
    --outdir 09_peaks/DYS \
    --keep-dup all \
    -p 0.01 \
    --call-summits \
    2> 09_peaks/DYS/DYS_rep${REP}_macs3.log
done

echo "=== MACS3 peak counts ==="
for REP in 1 2 3
do
  echo -n "DYS_rep${REP}: "
  wc -l < 09_peaks/DYS/DYS_rep${REP}_peaks.narrowPeak
done

# ============================================================
# STAGE 2, SORT BY SIGNAL SCORE FOR IDR
# ============================================================

for REP in 1 2 3
do
  sort -k8,8nr 09_peaks/DYS/DYS_rep${REP}_peaks.narrowPeak \
    > 10_idr/DYS/DYS_rep${REP}.sorted.narrowPeak
done

# ============================================================
# STAGE 3, PAIRWISE IDR
# ============================================================

conda activate idr

idr \
  --samples \
    10_idr/DYS/DYS_rep1.sorted.narrowPeak \
    10_idr/DYS/DYS_rep2.sorted.narrowPeak \
  --input-file-type narrowPeak \
  --rank p.value \
  --idr-threshold 0.05 \
  --output-file 10_idr/DYS/DYS_1v2.txt \
  --plot \
  --log-output-file 10_idr/DYS/DYS_1v2.log

idr \
  --samples \
    10_idr/DYS/DYS_rep1.sorted.narrowPeak \
    10_idr/DYS/DYS_rep3.sorted.narrowPeak \
  --input-file-type narrowPeak \
  --rank p.value \
  --idr-threshold 0.05 \
  --output-file 10_idr/DYS/DYS_1v3.txt \
  --plot \
  --log-output-file 10_idr/DYS/DYS_1v3.log

idr \
  --samples \
    10_idr/DYS/DYS_rep2.sorted.narrowPeak \
    10_idr/DYS/DYS_rep3.sorted.narrowPeak \
  --input-file-type narrowPeak \
  --rank p.value \
  --idr-threshold 0.05 \
  --output-file 10_idr/DYS/DYS_2v3.txt \
  --plot \
  --log-output-file 10_idr/DYS/DYS_2v3.log

conda deactivate

echo "=== IDR model summaries ==="
for PAIR in 1v2 1v3 2v3
do
  echo "--- rep${PAIR%v*} vs rep${PAIR#*v} ---"
  grep "Final parameter values" 10_idr/DYS/DYS_${PAIR}.log
  grep "Number of peaks passing" 10_idr/DYS/DYS_${PAIR}.log
done

# ============================================================
# STAGE 4, PAIRWISE MERGED BEDS (globalIDR <= 0.05)
# Column 12 >= 1.301 = globalIDR <= 0.05 in -log10 scale
# Falls back to col 5 >= 540 if col 12 absent
# Saved permanently to 10_idr/DYS/, not /tmp
# ============================================================

NCOLS=$(awk '{print NF; exit}' 10_idr/DYS/DYS_1v2.txt)

if [ "$NCOLS" -ge 12 ]; then
  echo "Using column 12 (globalIDR -log10 scale)"
  FILTER='$12 >= 1.301'
else
  echo "Column 12 absent — using column 5 (scaled score >= 540)"
  FILTER='$5 >= 540'
fi

for PAIR in 1v2 1v3 2v3
do
  awk "${FILTER} {print \$1\"\t\"\$2\"\t\"\$3}" \
    10_idr/DYS/DYS_${PAIR}.txt \
    | sort -k1,1 -k2,2n \
    | bedtools merge -i stdin \
    > 10_idr/DYS/DYS_${PAIR}.merged.bed
done

echo "=== IDR-passing peaks per pair ==="
for PAIR in 1v2 1v3 2v3
do
  echo -n "DYS_${PAIR}: "
  wc -l < 10_idr/DYS/DYS_${PAIR}.merged.bed
done

# ============================================================
# STAGE 5, CONSENSUS PEAKS (>=2 of 3 pairwise comparisons)
# Uses cat + merge + count method (proven approach from TGO Chapter 3)
# ============================================================

cat \
  10_idr/DYS/DYS_1v2.merged.bed \
  10_idr/DYS/DYS_1v3.merged.bed \
  10_idr/DYS/DYS_2v3.merged.bed \
  | sort -k1,1 -k2,2n \
  | bedtools merge -i stdin -c 1 -o count \
  | awk '$4 >= 2 {print $1"\t"$2"\t"$3}' \
  > 11_final/DYS_IDR_consensus.bed

echo -n "DYS consensus peaks (>=2/3, IDR <= 0.05): "
wc -l < 11_final/DYS_IDR_consensus.bed

# ============================================================
# STAGE 6, MEAN SIGNAL FROM MERGED BAM
# Uses samtools bedcov for accuracy with paired-end data
# ============================================================

# Sort and index merged BAM if not already done
if [ ! -f 07_blacklist/DYS_merged.clean.bam.bai ]; then
  samtools sort -@ 8 \
    -o 07_blacklist/DYS_merged.clean.sorted.bam \
    07_blacklist/DYS_merged.clean.bam
  samtools index 07_blacklist/DYS_merged.clean.sorted.bam
  MERGED_BAM="07_blacklist/DYS_merged.clean.sorted.bam"
else
  MERGED_BAM="07_blacklist/DYS_merged.clean.bam"
fi

# Sort consensus BED to match BAM chromosome order
sort -k1,1 -k2,2n \
  11_final/DYS_IDR_consensus.bed \
  > 11_final/DYS_IDR_consensus.sorted.bed

samtools bedcov \
  11_final/DYS_IDR_consensus.sorted.bed \
  ${MERGED_BAM} \
  | awk 'BEGIN{OFS="\t"} {print $1, $2, $3, $NF / ($3 - $2)}' \
  > 11_final/DYS_IDR_with_signal.txt

echo "Head of signal file:"
head -3 11_final/DYS_IDR_with_signal.txt

# ============================================================
# STAGE 7, FORMAT AS NARROWPEAK
# cols: chr start end name score strand signal pval qval summit
# summit = midpoint offset from peak start
# ============================================================

awk 'BEGIN{OFS="\t"}{
  print $1, $2, $3, ".", 1000, ".", $4, -1, -1, int(($3-$2)/2)
}' 11_final/DYS_IDR_with_signal.txt \
  > 11_final/DYS_IDR0.05.clean.narrowPeak

# ============================================================
# STAGE 8, CANONICAL CHROMOSOME FILTER
# Retain chr2L chr2R chr3L chr3R chrX chr4 only
# ============================================================

grep -E '^chr[234X][LR]?[[:space:]]' \
  11_final/DYS_IDR0.05.clean.narrowPeak \
  > 11_final/DYS_IDR0.05.final.narrowPeak

echo "=== Final peak counts ==="
echo -n "Before chr filter: "
wc -l < 11_final/DYS_IDR0.05.clean.narrowPeak
echo -n "After chr filter:  "
wc -l < 11_final/DYS_IDR0.05.final.narrowPeak

# ============================================================
# STAGE 9, FRIP ON FINAL PEAK SET
# ============================================================

for REP in 1 2 3
do
  echo -n "DYS_rep${REP} FRiP: "
  bedtools bamtobed \
    -i 07_blacklist/DYS_rep${REP}.clean.bam \
    > /tmp/DYS_frip_rep${REP}.bed
  total=$(wc -l < /tmp/DYS_frip_rep${REP}.bed)
  inpeaks=$(bedtools intersect -u \
    -a /tmp/DYS_frip_rep${REP}.bed \
    -b 11_final/DYS_IDR0.05.final.narrowPeak | wc -l)
  awk -v a="$inpeaks" -v b="$total" \
    'BEGIN { printf "%.4f\n", a/b }'
  rm /tmp/DYS_frip_rep${REP}.bed
done

# ============================================================
# STAGE 10, QC CHECKS
# ============================================================

echo ""
echo "=== Signal range check ==="
awk '{print $7}' 11_final/DYS_IDR0.05.final.narrowPeak \
  | sort -n \
  | awk 'NR==1{min=$1} END{print "min: "min"  max: "$1"  peaks: "NR}'

echo ""
echo "=== Final narrowPeak structure (first 3 lines) ==="
head -3 11_final/DYS_IDR0.05.final.narrowPeak

echo ""
echo "=== PIPELINE COMPLETE ==="
echo "Consensus BED:      11_final/DYS_IDR_consensus.bed"
echo "Final narrowPeak:   11_final/DYS_IDR0.05.final.narrowPeak"
wc -l 11_final/DYS_IDR0.05.final.narrowPeak