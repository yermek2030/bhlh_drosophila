#!/usr/bin/env bash
# ADAPTATION_OF_CONSENSUS_PEAKS_OF_SIM_TO_R_PIPELINE---------------------------

bedtools coverage \
-a 11_final/SIM_IDR_consensus.bed \
-b 07_blacklist/SIM_merged.clean.bam \
-mean \
> 11_final/SIM_IDR_with_signal.txt

echo "Head of signal file:"
head -3 11_final/SIM_IDR_with_signal.txt

awk 'BEGIN{OFS="\t"}{
  print $1, $2, $3, ".", 1000, ".", $4, -1, -1, int(($2+$3)/2)
}' 11_final/SIM_IDR_with_signal.txt \
> 11_final/SIM_IDR0.05.clean.narrowPeak

grep -E '^chr[234X]' 11_final/SIM_IDR0.05.clean.narrowPeak \
> 11_final/SIM_IDR0.05.final.narrowPeak

echo "=== Final peak counts ==="
echo -n "Before chr filter: "
wc -l < 11_final/SIM_IDR0.05.clean.narrowPeak
echo -n "After chr filter:  "
wc -l < 11_final/SIM_IDR0.05.final.narrowPeak
# Head of signal file:
# chr2L   223088  223764  54.9526634
# chr2L   267404  268097  37.8759003
# chr2L   512028  512435  36.2972984
# === Final peak counts ===
# Before chr filter: 556
# After chr filter:  556

for REP in 1 2
do
echo -n "SIM_rep${REP} FRiP: "
bedtools bamtobed -i 07_blacklist/SIM_rep${REP}.clean.bam > tmp_SIM.bed
total=$(wc -l < tmp_SIM.bed)
inpeaks=$(bedtools intersect -u \
          -a tmp_SIM.bed \
          -b 11_final/SIM_IDR0.05.final.narrowPeak | wc -l)
awk -v a="$inpeaks" -v b="$total" 'BEGIN { printf "%.4f\n", a/b }'
rm tmp_SIM.bed
done
# SIM_rep1 FRiP: 0.0084
# SIM_rep2 FRiP: 0.0092

echo ""
echo "=== Signal range check ==="
awk '{print $7}' 11_final/SIM_IDR0.05.final.narrowPeak \
| sort -n \
| awk 'NR==1{min=$1} END{print "min: "min"  max: "$1"  peaks: "NR}'

echo ""
echo "=== Final narrowPeak structure (first 3 lines) ==="
head -3 11_final/SIM_IDR0.05.final.narrowPeak

echo ""
echo "=== PIPELINE COMPLETE ==="
echo "Final SIM peak set: 11_final/SIM_IDR0.05.final.narrowPeak"

wc -l 11_final/SIM_IDR0.05.final.narrowPeak
