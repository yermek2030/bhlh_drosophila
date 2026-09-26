#!/usr/bin/env bash
# ============================================================
# TGO peak calling, IDR, consensus, signal narrowPeak
# Paired-end; 3 biological replicates; p=0.01
# IDR pairwise (1v2, 1v3, 2v3); consensus = peaks in >=2/3 pairs, from merged overlapping intervals
# Output: 11_final/TGO_IDR_consensus.bed, 11_final/TGO_IDR0.05.final.narrowPeak
# ============================================================

set -euo pipefail

mkdir -p 09_peaks/TGO
mkdir -p 10_idr/TGO
mkdir -p 11_final

# ============================================================
# STAGE 1, MACS3 PEAK CALLING (p=0.01, paired-end BAMPE)
# ============================================================

macs3 callpeak \
  -t 07_blacklist/TGO_rep1.clean.bam \
  -c 07_blacklist/TGO_control.clean.bam \
  -f BAMPE -g dm \
  -n TGO_rep1 \
  --outdir 09_peaks/TGO \
  --keep-dup all \
  -p 0.01 \
  --call-summits \
  2> 09_peaks/TGO/TGO_rep1_macs3.log

macs3 callpeak \
  -t 07_blacklist/TGO_rep2.clean.bam \
  -c 07_blacklist/TGO_control.clean.bam \
  -f BAMPE -g dm \
  -n TGO_rep2 \
  --outdir 09_peaks/TGO \
  --keep-dup all \
  -p 0.01 \
  --call-summits \
  2> 09_peaks/TGO/TGO_rep2_macs3.log

macs3 callpeak \
  -t 07_blacklist/TGO_rep3.clean.bam \
  -c 07_blacklist/TGO_control.clean.bam \
  -f BAMPE -g dm \
  -n TGO_rep3 \
  --outdir 09_peaks/TGO \
  --keep-dup all \
  -p 0.01 \
  --call-summits \
  2> 09_peaks/TGO/TGO_rep3_macs3.log

echo "=== MACS3 peak counts ==="
echo -n "TGO_rep1: "
wc -l < 09_peaks/TGO/TGO_rep1_peaks.narrowPeak
echo -n "TGO_rep2: "
wc -l < 09_peaks/TGO/TGO_rep2_peaks.narrowPeak
echo -n "TGO_rep3: "
wc -l < 09_peaks/TGO/TGO_rep3_peaks.narrowPeak

# ============================================================
# STAGE 2, SORT BY SIGNAL SCORE FOR IDR
# ============================================================

sort -k8,8nr 09_peaks/TGO/TGO_rep1_peaks.narrowPeak \
  > 10_idr/TGO/TGO_rep1.sorted.narrowPeak

sort -k8,8nr 09_peaks/TGO/TGO_rep2_peaks.narrowPeak \
  > 10_idr/TGO/TGO_rep2.sorted.narrowPeak

sort -k8,8nr 09_peaks/TGO/TGO_rep3_peaks.narrowPeak \
  > 10_idr/TGO/TGO_rep3.sorted.narrowPeak

# ============================================================
# STAGE 3, PAIRWISE IDR
# ============================================================

conda activate idr

idr \
  --samples \
    10_idr/TGO/TGO_rep1.sorted.narrowPeak \
    10_idr/TGO/TGO_rep2.sorted.narrowPeak \
  --input-file-type narrowPeak \
  --rank p.value \
  --idr-threshold 0.05 \
  --output-file 10_idr/TGO/TGO_1v2.txt \
  --plot \
  --log-output-file 10_idr/TGO/TGO_1v2.log

idr \
  --samples \
    10_idr/TGO/TGO_rep1.sorted.narrowPeak \
    10_idr/TGO/TGO_rep3.sorted.narrowPeak \
  --input-file-type narrowPeak \
  --rank p.value \
  --idr-threshold 0.05 \
  --output-file 10_idr/TGO/TGO_1v3.txt \
  --plot \
  --log-output-file 10_idr/TGO/TGO_1v3.log

idr \
  --samples \
    10_idr/TGO/TGO_rep2.sorted.narrowPeak \
    10_idr/TGO/TGO_rep3.sorted.narrowPeak \
  --input-file-type narrowPeak \
  --rank p.value \
  --idr-threshold 0.05 \
  --output-file 10_idr/TGO/TGO_2v3.txt \
  --plot \
  --log-output-file 10_idr/TGO/TGO_2v3.log

conda deactivate

echo "=== IDR model summaries ==="
echo "--- rep1 vs rep2 ---"
grep "Final parameter values" 10_idr/TGO/TGO_1v2.log
grep "Number of peaks passing" 10_idr/TGO/TGO_1v2.log
echo "--- rep1 vs rep3 ---"
grep "Final parameter values" 10_idr/TGO/TGO_1v3.log
grep "Number of peaks passing" 10_idr/TGO/TGO_1v3.log
echo "--- rep2 vs rep3 ---"
grep "Final parameter values" 10_idr/TGO/TGO_2v3.log
grep "Number of peaks passing" 10_idr/TGO/TGO_2v3.log

# ============================================================
# STAGE 4, PAIRWISE MERGED BEDS (globalIDR <= 0.05)
# Saved permanently to 10_idr/TGO/, not /tmp
# ============================================================

NCOLS=$(awk '{print NF; exit}' 10_idr/TGO/TGO_1v2.txt)

if [ "$NCOLS" -ge 12 ]; then
  echo "Using column 12 (globalIDR -log10 scale)"

  awk '$12 >= 1.301 {print $1"\t"$2"\t"$3}' \
    10_idr/TGO/TGO_1v2.txt \
    | sort -k1,1 -k2,2n \
    | bedtools merge -i stdin \
    > 10_idr/TGO/TGO_1v2.merged.bed

  awk '$12 >= 1.301 {print $1"\t"$2"\t"$3}' \
    10_idr/TGO/TGO_1v3.txt \
    | sort -k1,1 -k2,2n \
    | bedtools merge -i stdin \
    > 10_idr/TGO/TGO_1v3.merged.bed

  awk '$12 >= 1.301 {print $1"\t"$2"\t"$3}' \
    10_idr/TGO/TGO_2v3.txt \
    | sort -k1,1 -k2,2n \
    | bedtools merge -i stdin \
    > 10_idr/TGO/TGO_2v3.merged.bed

else
  echo "Column 12 absent — using column 5 (scaled score >= 540)"

  awk '$5 >= 540 {print $1"\t"$2"\t"$3}' \
    10_idr/TGO/TGO_1v2.txt \
    | sort -k1,1 -k2,2n \
    | bedtools merge -i stdin \
    > 10_idr/TGO/TGO_1v2.merged.bed

  awk '$5 >= 540 {print $1"\t"$2"\t"$3}' \
    10_idr/TGO/TGO_1v3.txt \
    | sort -k1,1 -k2,2n \
    | bedtools merge -i stdin \
    > 10_idr/TGO/TGO_1v3.merged.bed

  awk '$5 >= 540 {print $1"\t"$2"\t"$3}' \
    10_idr/TGO/TGO_2v3.txt \
    | sort -k1,1 -k2,2n \
    | bedtools merge -i stdin \
    > 10_idr/TGO/TGO_2v3.merged.bed
fi

echo "=== IDR-passing peaks per pair ==="
echo -n "TGO_1v2: "
wc -l < 10_idr/TGO/TGO_1v2.merged.bed
echo -n "TGO_1v3: "
wc -l < 10_idr/TGO/TGO_1v3.merged.bed
echo -n "TGO_2v3: "
wc -l < 10_idr/TGO/TGO_2v3.merged.bed

# ============================================================
# STAGE 5, CONSENSUS PEAKS (>=2 of 3 pairwise comparisons)
# Uses a cat+merge+count method
# ============================================================

cat \
  10_idr/TGO/TGO_1v2.merged.bed \
  10_idr/TGO/TGO_1v3.merged.bed \
  10_idr/TGO/TGO_2v3.merged.bed \
  | sort -k1,1 -k2,2n \
  | bedtools merge -i stdin -c 1 -o count \
  | awk '$4 >= 2 {print $1"\t"$2"\t"$3}' \
  > 11_final/TGO_IDR_consensus.bed

echo -n "TGO consensus peaks (>=2/3, IDR <= 0.05): "
wc -l < 11_final/TGO_IDR_consensus.bed

# ============================================================
# STAGE 6, MEAN SIGNAL FROM MERGED BAM
# Uses samtools bedcov for accuracy with paired-end data
# ============================================================

# Sort and index merged BAM if not already done
if [ ! -f 07_blacklist/TGO_merged.clean.bam.bai ]; then
  samtools sort -@ 8 \
    -o 07_blacklist/TGO_merged.clean.sorted.bam \
    07_blacklist/TGO_merged.clean.bam
  samtools index 07_blacklist/TGO_merged.clean.sorted.bam
  MERGED_BAM="07_blacklist/TGO_merged.clean.sorted.bam"
else
  MERGED_BAM="07_blacklist/TGO_merged.clean.bam"
fi

# Sort consensus BED to match BAM chromosome order
sort -k1,1 -k2,2n \
  11_final/TGO_IDR_consensus.bed \
  > 11_final/TGO_IDR_consensus.sorted.bed

samtools bedcov \
  11_final/TGO_IDR_consensus.sorted.bed \
  ${MERGED_BAM} \
  | awk 'BEGIN{OFS="\t"} {print $1, $2, $3, $NF / ($3 - $2)}' \
  > 11_final/TGO_IDR_with_signal.txt

echo "Head of signal file:"
head -3 11_final/TGO_IDR_with_signal.txt

# ============================================================
# STAGE 7, FORMAT AS NARROWPEAK
# ============================================================

awk 'BEGIN{OFS="\t"}{
  print $1, $2, $3, ".", 1000, ".", $4, -1, -1, int(($3-$2)/2)
}' 11_final/TGO_IDR_with_signal.txt \
  > 11_final/TGO_IDR0.05.clean.narrowPeak

# ============================================================
# STAGE 8, CANONICAL CHROMOSOME FILTER
# ============================================================

grep -E '^chr[234X][LR]?[[:space:]]' \
  11_final/TGO_IDR0.05.clean.narrowPeak \
  > 11_final/TGO_IDR0.05.final.narrowPeak

echo "=== Final peak counts ==="
echo -n "Before chr filter: "
wc -l < 11_final/TGO_IDR0.05.clean.narrowPeak
echo -n "After chr filter:  "
wc -l < 11_final/TGO_IDR0.05.final.narrowPeak

# ============================================================
# STAGE 9, FRIP ON FINAL PEAK SET
# ============================================================

for REP in 1 2 3
do
  echo -n "TGO_rep${REP} FRiP: "
  bedtools bamtobed \
    -i 07_blacklist/TGO_rep${REP}.clean.bam \
    > /tmp/TGO_frip_rep${REP}.bed
  total=$(wc -l < /tmp/TGO_frip_rep${REP}.bed)
  inpeaks=$(bedtools intersect -u \
    -a /tmp/TGO_frip_rep${REP}.bed \
    -b 11_final/TGO_IDR0.05.final.narrowPeak | wc -l)
  awk -v a="$inpeaks" -v b="$total" \
    'BEGIN { printf "%.4f\n", a/b }'
  rm /tmp/TGO_frip_rep${REP}.bed
done

# ============================================================
# STAGE 10, QC CHECKS
# ============================================================

echo ""
echo "=== Signal range check ==="
awk '{print $7}' 11_final/TGO_IDR0.05.final.narrowPeak \
  | sort -n \
  | awk 'NR==1{min=$1} END{print "min: "min"  max: "$1"  peaks: "NR}'

echo ""
echo "=== Final narrowPeak structure (first 3 lines) ==="
head -3 11_final/TGO_IDR0.05.final.narrowPeak

echo ""
echo "=== PIPELINE COMPLETE ==="
echo "Consensus BED:      11_final/TGO_IDR_consensus.bed"
echo "Final narrowPeak:   11_final/TGO_IDR0.05.final.narrowPeak"
wc -l 11_final/TGO_IDR0.05.final.narrowPeak
