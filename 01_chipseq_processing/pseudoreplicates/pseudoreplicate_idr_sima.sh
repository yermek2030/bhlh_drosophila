#!/usr/bin/env bash
# ============================================================
# Pseudoreplicate generation and IDR, SIMA (single-end)
# Nt/Npr/N1/N2 = true-rep, pooled-pr, rep1-pr, rep2-pr IDR peak counts
# RR = max(Nt,Npr)/min(Nt,Npr); SCR = max(N1,N2)/min(N1,N2)
# ENCODE: PASS both<=2, BORDERLINE one>2, FAIL both>2
# ============================================================

mkdir -p 12_pseudoreps/SIMA
mkdir -p 10_idr/SIMA_pr

# ------------------------------------------------------------
# STEP 1: POOLED PSEUDO-REPLICATES
# Merge rep1 + rep2, then randomly split 50/50
# ------------------------------------------------------------

# Merge replicates (single-end, no pair preservation needed)
samtools merge -f 12_pseudoreps/SIMA/SIMA_pooled.bam \
07_blacklist/SIMA_rep1.clean.bam \
07_blacklist/SIMA_rep2.clean.bam

samtools index 12_pseudoreps/SIMA/SIMA_pooled.bam

# Count total reads
TOTAL=$(samtools view -c 12_pseudoreps/SIMA/SIMA_pooled.bam)
HALF=$(( TOTAL / 2 ))
echo "Total SIMA pooled reads: $TOTAL"
echo "Split size per pseudo-rep: $HALF"

# Random split using samtools view with subsampling
# -s seed.fraction where integer part is seed, decimal part is fraction
samtools view -h -s 42.5 12_pseudoreps/SIMA/SIMA_pooled.bam \
| samtools view -b - \
> 12_pseudoreps/SIMA/SIMA_pooled_pr1.bam

# Complement: all reads not in pr1 = pr2
# Extract read IDs from pr1
samtools view 12_pseudoreps/SIMA/SIMA_pooled_pr1.bam \
| awk '{print $1}' | sort -u \
> 12_pseudoreps/SIMA/SIMA_pr1_read_ids.txt

# Pipe pooled BAM through picard FilterSamReads to exclude pr1 reads
picard -Djava.io.tmpdir=$PROJECT_DIR/tmp \
FilterSamReads \
I=12_pseudoreps/SIMA/SIMA_pooled.bam \
O=12_pseudoreps/SIMA/SIMA_pooled_pr2.bam \
READ_LIST_FILE=12_pseudoreps/SIMA/SIMA_pr1_read_ids.txt \
FILTER=excludeReadList

samtools index 12_pseudoreps/SIMA/SIMA_pooled_pr1.bam
samtools index 12_pseudoreps/SIMA/SIMA_pooled_pr2.bam

# Verify counts are roughly 50/50
echo -n "SIMA pooled_pr1 reads: "
samtools view -c 12_pseudoreps/SIMA/SIMA_pooled_pr1.bam
echo -n "SIMA pooled_pr2 reads: "
samtools view -c 12_pseudoreps/SIMA/SIMA_pooled_pr2.bam

# ------------------------------------------------------------
# STEP 2: SELF-CONSISTENCY PSEUDO-REPLICATES (per biological replicate)
# Each replicate is split 50/50 independently
# ------------------------------------------------------------

# rep1 split
samtools view -h -s 42.5 07_blacklist/SIMA_rep1.clean.bam \
| samtools view -b - \
> 12_pseudoreps/SIMA/SIMA_rep1_pr1.bam

samtools view 12_pseudoreps/SIMA/SIMA_rep1_pr1.bam \
| awk '{print $1}' | sort -u \
> 12_pseudoreps/SIMA/SIMA_rep1_pr1_ids.txt

picard -Djava.io.tmpdir=$PROJECT_DIR/tmp \
FilterSamReads \
I=07_blacklist/SIMA_rep1.clean.bam \
O=12_pseudoreps/SIMA/SIMA_rep1_pr2.bam \
READ_LIST_FILE=12_pseudoreps/SIMA/SIMA_rep1_pr1_ids.txt \
FILTER=excludeReadList

samtools index 12_pseudoreps/SIMA/SIMA_rep1_pr1.bam
samtools index 12_pseudoreps/SIMA/SIMA_rep1_pr2.bam

# rep2 split
samtools view -h -s 42.5 07_blacklist/SIMA_rep2.clean.bam \
| samtools view -b - \
> 12_pseudoreps/SIMA/SIMA_rep2_pr1.bam

samtools view 12_pseudoreps/SIMA/SIMA_rep2_pr1.bam \
| awk '{print $1}' | sort -u \
> 12_pseudoreps/SIMA/SIMA_rep2_pr1_ids.txt

picard -Djava.io.tmpdir=$PROJECT_DIR/tmp \
FilterSamReads \
I=07_blacklist/SIMA_rep2.clean.bam \
O=12_pseudoreps/SIMA/SIMA_rep2_pr2.bam \
READ_LIST_FILE=12_pseudoreps/SIMA/SIMA_rep2_pr1_ids.txt \
FILTER=excludeReadList

samtools index 12_pseudoreps/SIMA/SIMA_rep2_pr1.bam
samtools index 12_pseudoreps/SIMA/SIMA_rep2_pr2.bam

# ------------------------------------------------------------
# STEP 3: MACS3 ON ALL PSEUDO-REPLICATES
# Same parameters as production: p=0.01, --nomodel --extsize 225 (SPP-derived)
# ------------------------------------------------------------

# Pooled pseudo-replicates
macs3 callpeak \
-t 12_pseudoreps/SIMA/SIMA_pooled_pr1.bam \
-c 07_blacklist/SIMA_control.clean.bam \
-f BAM -g dm \
-n SIMA_pooled_pr1 \
--outdir 12_pseudoreps/SIMA \
--keep-dup all --nomodel --extsize 225 \
-p 0.01 --call-summits \
2> 12_pseudoreps/SIMA/SIMA_pooled_pr1_macs3.log

macs3 callpeak \
-t 12_pseudoreps/SIMA/SIMA_pooled_pr2.bam \
-c 07_blacklist/SIMA_control.clean.bam \
-f BAM -g dm \
-n SIMA_pooled_pr2 \
--outdir 12_pseudoreps/SIMA \
--keep-dup all --nomodel --extsize 225 \
-p 0.01 --call-summits \
2> 12_pseudoreps/SIMA/SIMA_pooled_pr2_macs3.log

# Self-consistency: rep1 pseudo-replicates
macs3 callpeak \
-t 12_pseudoreps/SIMA/SIMA_rep1_pr1.bam \
-c 07_blacklist/SIMA_control.clean.bam \
-f BAM -g dm \
-n SIMA_rep1_pr1 \
--outdir 12_pseudoreps/SIMA \
--keep-dup all --nomodel --extsize 225 \
-p 0.01 --call-summits \
2> 12_pseudoreps/SIMA/SIMA_rep1_pr1_macs3.log

macs3 callpeak \
-t 12_pseudoreps/SIMA/SIMA_rep1_pr2.bam \
-c 07_blacklist/SIMA_control.clean.bam \
-f BAM -g dm \
-n SIMA_rep1_pr2 \
--outdir 12_pseudoreps/SIMA \
--keep-dup all --nomodel --extsize 225 \
-p 0.01 --call-summits \
2> 12_pseudoreps/SIMA/SIMA_rep1_pr2_macs3.log

# Self-consistency: rep2 pseudo-replicates
macs3 callpeak \
-t 12_pseudoreps/SIMA/SIMA_rep2_pr1.bam \
-c 07_blacklist/SIMA_control.clean.bam \
-f BAM -g dm \
-n SIMA_rep2_pr1 \
--outdir 12_pseudoreps/SIMA \
--keep-dup all --nomodel --extsize 225 \
-p 0.01 --call-summits \
2> 12_pseudoreps/SIMA/SIMA_rep2_pr1_macs3.log

macs3 callpeak \
-t 12_pseudoreps/SIMA/SIMA_rep2_pr2.bam \
-c 07_blacklist/SIMA_control.clean.bam \
-f BAM -g dm \
-n SIMA_rep2_pr2 \
--outdir 12_pseudoreps/SIMA \
--keep-dup all --nomodel --extsize 225 \
-p 0.01 --call-summits \
2> 12_pseudoreps/SIMA/SIMA_rep2_pr2_macs3.log

# ------------------------------------------------------------
# STEP 4: SORT ALL PSEUDO-REPLICATE PEAKS BY SIGNAL SCORE
# ------------------------------------------------------------

for NAME in SIMA_pooled_pr1 SIMA_pooled_pr2 \
SIMA_rep1_pr1 SIMA_rep1_pr2 \
SIMA_rep2_pr1 SIMA_rep2_pr2; do
sort -k8,8nr 12_pseudoreps/SIMA/${NAME}_peaks.narrowPeak \
> 12_pseudoreps/SIMA/${NAME}.sorted.narrowPeak
done

# ------------------------------------------------------------
# STEP 5: IDR ON PSEUDO-REPLICATES
# ------------------------------------------------------------

conda activate idr

# Pooled pseudo-rep IDR (Npr)
idr \
--samples \
12_pseudoreps/SIMA/SIMA_pooled_pr1.sorted.narrowPeak \
12_pseudoreps/SIMA/SIMA_pooled_pr2.sorted.narrowPeak \
--input-file-type narrowPeak \
--rank p.value --idr-threshold 0.05 \
--output-file 10_idr/SIMA_pr/SIMA_pooled_pr.txt \
--plot \
--log-output-file 10_idr/SIMA_pr/SIMA_pooled_pr.log

# Self-consistency rep1 IDR (N1)
idr \
--samples \
12_pseudoreps/SIMA/SIMA_rep1_pr1.sorted.narrowPeak \
12_pseudoreps/SIMA/SIMA_rep1_pr2.sorted.narrowPeak \
--input-file-type narrowPeak \
--rank p.value --idr-threshold 0.05 \
--output-file 10_idr/SIMA_pr/SIMA_rep1_selfpr.txt \
--plot \
--log-output-file 10_idr/SIMA_pr/SIMA_rep1_selfpr.log

# Self-consistency rep2 IDR (N2)
idr \
--samples \
12_pseudoreps/SIMA/SIMA_rep2_pr1.sorted.narrowPeak \
12_pseudoreps/SIMA/SIMA_rep2_pr2.sorted.narrowPeak \
--input-file-type narrowPeak \
--rank p.value --idr-threshold 0.05 \
--output-file 10_idr/SIMA_pr/SIMA_rep2_selfpr.txt \
--plot \
--log-output-file 10_idr/SIMA_pr/SIMA_rep2_selfpr.log

conda deactivate

# ------------------------------------------------------------
# STEP 6: COMPUTE RESCUE AND SELF-CONSISTENCY RATIOS
# ------------------------------------------------------------

# Detect column count
NCOLS=$(awk '{print NF; exit}' 10_idr/SIMA_pr/SIMA_pooled_pr.txt)

if [ "$NCOLS" -ge 12 ]; then
FILTER='$12 >= 1.301'
else
  FILTER='$5 >= 540'
fi

# Count peaks passing IDR in each analysis
Nt=$(awk "$FILTER" 10_idr/SIMA/SIMA_1v2.txt | wc -l)
Npr=$(awk "$FILTER" 10_idr/SIMA_pr/SIMA_pooled_pr.txt | wc -l)
N1=$(awk "$FILTER" 10_idr/SIMA_pr/SIMA_rep1_selfpr.txt | wc -l)
N2=$(awk "$FILTER" 10_idr/SIMA_pr/SIMA_rep2_selfpr.txt | wc -l)

# Compute ratios
RR=$(awk -v a="$Nt" -v b="$Npr" 'BEGIN{
  if (a > b) print a/b; else print b/a
}')

SCR=$(awk -v a="$N1" -v b="$N2" 'BEGIN{
  if (a > b) print a/b; else print b/a
}')

# Output summary
echo "=== ENCODE Pseudo-Replicate Diagnostics: SIMA ==="
echo "Nt  (true replicate IDR, rep1 vs rep2): $Nt"
echo "Npr (pooled pseudo-replicate IDR):      $Npr"
echo "N1  (rep1 self-consistency IDR):        $N1"
echo "N2  (rep2 self-consistency IDR):        $N2"
echo ""
printf "Rescue Ratio (RR)           = max(Nt,Npr)/min = %.3f\n" $RR
printf "Self-Consistency Ratio (SCR) = max(N1,N2)/min = %.3f\n" $SCR
echo ""
echo "PASS criteria: both ratios <= 2"
echo "BORDERLINE:    one ratio > 2"
echo "FAIL:          both ratios > 2"

# Save to file for chapter text
cat > 12_pseudoreps/SIMA/SIMA_pr_diagnostics.txt <<EOF
SIMA ENCODE Pseudo-Replicate Diagnostics
========================================
  Nt  = $Nt
Npr = $Npr
N1  = $N1
N2  = $N2
RR  = $RR
SCR = $SCR
EOF

