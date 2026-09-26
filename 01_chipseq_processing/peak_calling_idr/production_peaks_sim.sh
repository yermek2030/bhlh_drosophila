#!/usr/bin/env bash
# ============================================================
# SIM PRODUCTION PEAK CALLING, IDR, CONSENSUS, SIGNAL NARROWPEAK
# Paired-end; 2 biological replicates; p=0.01 (empirically selected)
# MACS3 uses -f BAMPE, no --nomodel, no --extsize needed
# Output: 11_final/SIM_IDR_consensus.bed
#         11_final/SIM_IDR0.05.final.narrowPeak (for R/GRanges import)
# ============================================================

set -euo pipefail

mkdir -p 09_peaks/SIM
mkdir -p 10_idr/SIM
mkdir -p 11_final

# ============================================================
# STAGE 1, MACS3 PEAK CALLING (p=0.01, paired-end BAMPE)
# ============================================================

macs3 callpeak \
  -t 07_blacklist/SIM_rep1.clean.bam \
  -c 07_blacklist/SIM_control.clean.bam \
  -f BAMPE -g dm \
  -n SIM_rep1 \
  --outdir 09_peaks/SIM \
  --keep-dup all \
  -p 0.01 \
  --call-summits \
  2> 09_peaks/SIM/SIM_rep1_macs3.log

macs3 callpeak \
  -t 07_blacklist/SIM_rep2.clean.bam \
  -c 07_blacklist/SIM_control.clean.bam \
  -f BAMPE -g dm \
  -n SIM_rep2 \
  --outdir 09_peaks/SIM \
  --keep-dup all \
  -p 0.01 \
  --call-summits \
  2> 09_peaks/SIM/SIM_rep2_macs3.log

echo "=== MACS3 peak counts ==="
echo -n "SIM_rep1: "
wc -l < 09_peaks/SIM/SIM_rep1_peaks.narrowPeak
echo -n "SIM_rep2: "
wc -l < 09_peaks/SIM/SIM_rep2_peaks.narrowPeak

# ============================================================
# STAGE 2, SORT BY SIGNAL SCORE FOR IDR
# ============================================================

sort -k8,8nr 09_peaks/SIM/SIM_rep1_peaks.narrowPeak \
  > 10_idr/SIM/SIM_rep1.sorted.narrowPeak

sort -k8,8nr 09_peaks/SIM/SIM_rep2_peaks.narrowPeak \
  > 10_idr/SIM/SIM_rep2.sorted.narrowPeak

# ============================================================
# STAGE 3, IDR (rep1 vs rep2)
# ============================================================

conda activate idr

idr \
  --samples \
    10_idr/SIM/SIM_rep1.sorted.narrowPeak \
    10_idr/SIM/SIM_rep2.sorted.narrowPeak \
  --input-file-type narrowPeak \
  --rank p.value \
  --idr-threshold 0.05 \
  --output-file 10_idr/SIM/SIM_1v2.txt \
  --plot \
  --log-output-file 10_idr/SIM/SIM_1v2.log

conda deactivate

echo "=== IDR model summary ==="
grep "Final parameter values" 10_idr/SIM/SIM_1v2.log
grep "Number of peaks passing" 10_idr/SIM/SIM_1v2.log

# ============================================================
# STAGE 4, CONSENSUS PEAK SET (globalIDR <= 0.05)
# ============================================================

NCOLS=$(awk '{print NF; exit}' 10_idr/SIM/SIM_1v2.txt)

if [ "$NCOLS" -ge 12 ]; then
  echo "Using column 12 (globalIDR -log10 scale)"
  awk '$12 >= 1.301 {print $1"\t"$2"\t"$3}' \
    10_idr/SIM/SIM_1v2.txt \
    | sort -k1,1 -k2,2n \
    | bedtools merge -i stdin \
    > 11_final/SIM_IDR_consensus.bed
else
  echo "Column 12 absent — using column 5 (scaled score >= 540)"
  awk '$5 >= 540 {print $1"\t"$2"\t"$3}' \
    10_idr/SIM/SIM_1v2.txt \
    | sort -k1,1 -k2,2n \
    | bedtools merge -i stdin \
    > 11_final/SIM_IDR_consensus.bed
fi

echo -n "SIM consensus peaks (IDR <= 0.05): "
wc -l < 11_final/SIM_IDR_consensus.bed

# ============================================================
# STAGE 5, MEAN SIGNAL FROM MERGED BAM
# ============================================================

bedtools coverage \
  -a 11_final/SIM_IDR_consensus.bed \
  -b 07_blacklist/SIM_merged.clean.bam \
  -mean \
  > 11_final/SIM_IDR_with_signal.txt

echo "Head of signal file:"
head -3 11_final/SIM_IDR_with_signal.txt

# ============================================================
# STAGE 6, FORMAT AS NARROWPEAK
# ============================================================

awk 'BEGIN{OFS="\t"}{
  print $1, $2, $3, ".", 1000, ".", $4, -1, -1, int(($3-$2)/2)
}' 11_final/SIM_IDR_with_signal.txt \
  > 11_final/SIM_IDR0.05.clean.narrowPeak

# ============================================================
# STAGE 7, CANONICAL CHROMOSOME FILTER
# ============================================================

grep -E '^chr[234X][LR]?[[:space:]]' \
  11_final/SIM_IDR0.05.clean.narrowPeak \
  > 11_final/SIM_IDR0.05.final.narrowPeak

echo "=== Final peak counts ==="
echo -n "Before chr filter: "
wc -l < 11_final/SIM_IDR0.05.clean.narrowPeak
echo -n "After chr filter:  "
wc -l < 11_final/SIM_IDR0.05.final.narrowPeak

# ============================================================
# STAGE 8, FRIP ON FINAL PEAK SET
# ============================================================

for REP in 1 2
do
  echo -n "SIM_rep${REP} FRiP: "
  bedtools bamtobed \
    -i 07_blacklist/SIM_rep${REP}.clean.bam \
    > /tmp/SIM_frip_rep${REP}.bed
  total=$(wc -l < /tmp/SIM_frip_rep${REP}.bed)
  inpeaks=$(bedtools intersect -u \
    -a /tmp/SIM_frip_rep${REP}.bed \
    -b 11_final/SIM_IDR0.05.final.narrowPeak | wc -l)
  awk -v a="$inpeaks" -v b="$total" \
    'BEGIN { printf "%.4f\n", a/b }'
  rm /tmp/SIM_frip_rep${REP}.bed
done

# ============================================================
# STAGE 9, QC CHECKS
# ============================================================

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
echo "Consensus BED:      11_final/SIM_IDR_consensus.bed"
echo "Final narrowPeak:   11_final/SIM_IDR0.05.final.narrowPeak"
wc -l 11_final/SIM_IDR0.05.final.narrowPeak
