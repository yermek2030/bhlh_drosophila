#!/usr/bin/env bash
# ADAPTATION_OF_CONSENSUS_PEAKS_OF_TRH_TO_R_PIPELINE

bedtools coverage \
-a 11_final/TRH_IDR_consensus.bed \
-b 07_blacklist/TRH_merged.clean.bam \
-mean \
> 11_final/TRH_IDR_with_signal.txt

echo "Head of signal file:"
head -3 11_final/TRH_IDR_with_signal.txt

awk 'BEGIN{OFS="\t"}{
  print $1, $2, $3, ".", 1000, ".", $4, -1, -1, int(($2+$3)/2)
}' 11_final/TRH_IDR_with_signal.txt \
> 11_final/TRH_IDR0.05.clean.narrowPeak

grep -E '^chr[234X]' 11_final/TRH_IDR0.05.clean.narrowPeak \
> 11_final/TRH_IDR0.05.final.narrowPeak

echo "=== Final peak counts ==="
echo -n "Before chr filter: "
wc -l < 11_final/TRH_IDR0.05.clean.narrowPeak
echo -n "After chr filter:  "
wc -l < 11_final/TRH_IDR0.05.final.narrowPeak

for REP in 1 2
do
echo -n "TRH_rep${REP} FRiP: "
bedtools bamtobed -i 07_blacklist/TRH_rep${REP}.clean.bam > tmp_trh.bed
total=$(wc -l < tmp_trh.bed)
inpeaks=$(bedtools intersect -u \
          -a tmp_trh.bed \
          -b 11_final/TRH_IDR0.05.final.narrowPeak | wc -l)
awk -v a="$inpeaks" -v b="$total" 'BEGIN { printf "%.4f\n", a/b }'
rm tmp_trh.bed
done
# TRH_rep1 FRiP: 0.0219
# TRH_rep2 FRiP: 0.0185

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
echo "Final TRH peak set: 11_final/TRH_IDR0.05.final.narrowPeak"

wc -l 11_final/TRH_IDR0.05.final.narrowPeak
