#!/usr/bin/env bash
# ADAPTATION_OF_CONSENSUS_PEAKS_OF_TGO_TO_R_PIPELINE

# Step 1: sort the consensus BED to match BAM chromosome order
sort -k1,1 -k2,2n \
11_final/TGO_IDR_consensus.bed \
> 09_peaks/TGO/TGO_idr_consensus.sorted.bed

# Step 2: confirm merged BAM is sorted and indexed
# (skip if TGO_merged.clean.bam.bai already exists)
samtools sort -@ 8 -o 07_blacklist/TGO_merged.clean.sorted.bam \
07_blacklist/TGO_merged.clean.bam
samtools index 07_blacklist/TGO_merged.clean.sorted.bam

# Step 3: samtools bedcov (returns sum of base coverage per interval)
# then awk divides by interval length to get mean
samtools bedcov \
09_peaks/TGO/TGO_idr_consensus.sorted.bed \
07_blacklist/TGO_merged.clean.sorted.bam \
| awk 'BEGIN{OFS="\t"} {print $1, $2, $3, $NF / ($3 - $2)}' \
> 11_final/TGO_IDR_with_signal.txt


awk 'BEGIN{OFS="\t"}{
  print $1, $2, $3, ".", 1000, ".", $4, -1, -1, int(($2+$3)/2)
}' 11_final/TGO_IDR_with_signal.txt \
> 11_final/TGO_IDR0.05.clean.narrowPeak  


grep -E '^chr[234X]' 11_final/TGO_IDR0.05.clean.narrowPeak \
> 11_final/TGO_IDR0.05.final.narrowPeak
echo "=== Final peak counts ==="
echo -n "Before chr filter: "
wc -l < 11_final/TGO_IDR0.05.clean.narrowPeak
echo -n "After chr filter:  "
wc -l < 11_final/TGO_IDR0.05.final.narrowPeak

echo ""
echo "=== Signal range check ==="
awk '{print $7}' 11_final/TGO_IDR0.05.final.narrowPeak \
| sort -n \
| awk 'NR==1{min=$1} END{print "min: "min"  max: "$1"  peaks: "NR}'
# Expected: min > 0, max < 10000, peaks = final count

echo ""
echo "=== Final narrowPeak structure (first 3 lines) ==="
head -3 11_final/TGO_IDR0.05.final.narrowPeak

echo ""
echo "=== PIPELINE COMPLETE ==="
echo "Final TGO peak set: 11_final/TGO_IDR0.05.final.narrowPeak"

wc -l 11_final/TGO_IDR0.05.final.narrowPeak
