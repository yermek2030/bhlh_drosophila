#!/usr/bin/env bash
  

set -euo pipefail
mkdir -p 12_pseudoreps/SIM
mkdir -p 10_idr/SIM_pr

# Merge replicates first
samtools merge -f 12_pseudoreps/SIM/SIM_pooled.bam \
07_blacklist/SIM_rep1.clean.bam \
07_blacklist/SIM_rep2.clean.bam

samtools index 12_pseudoreps/SIM/SIM_pooled.bam

# Name-sort so paired reads are adjacent, essential for correct pair-aware
# splitting
samtools sort -n -@ 8 12_pseudoreps/SIM/SIM_pooled.bam \
-o 12_pseudoreps/SIM/SIM_pooled.namesorted.bam

# Pair-aware splitting: assign by read_id so both mates go to same pseudo-rep
samtools view -h 12_pseudoreps/SIM/SIM_pooled.namesorted.bam \
| awk 'BEGIN{srand(42)}
    /^@/ {
      print > "/tmp/SIM_pooled_pr1.sam"
      print > "/tmp/SIM_pooled_pr2.sam"
      next
    }
    {
      read_id = $1
      if (!(read_id in assigned)) {
        assigned[read_id] = (rand() < 0.5 ? 1 : 2)
      }
      if (assigned[read_id] == 1) print > "/tmp/SIM_pooled_pr1.sam"
      else                        print > "/tmp/SIM_pooled_pr2.sam"
    }'

samtools view -bS /tmp/SIM_pooled_pr1.sam \
| samtools sort -@ 8 -o 12_pseudoreps/SIM/SIM_pooled_pr1.bam
samtools view -bS /tmp/SIM_pooled_pr2.sam \
| samtools sort -@ 8 -o 12_pseudoreps/SIM/SIM_pooled_pr2.bam

samtools index 12_pseudoreps/SIM/SIM_pooled_pr1.bam
samtools index 12_pseudoreps/SIM/SIM_pooled_pr2.bam
rm /tmp/SIM_pooled_pr1.sam /tmp/SIM_pooled_pr2.sam

# Verify split is ~50/50 and pair preservation
echo -n "SIM pooled_pr1 reads: "
samtools view -c 12_pseudoreps/SIM/SIM_pooled_pr1.bam
echo -n "SIM pooled_pr2 reads: "
samtools view -c 12_pseudoreps/SIM/SIM_pooled_pr2.bam

# Check that paired reads are preserved
# flagstat should show 0 singletons and 100% properly paired
echo "--- pooled_pr1 flagstat ---"
samtools flagstat 12_pseudoreps/SIM/SIM_pooled_pr1.bam | head -14
echo "--- pooled_pr2 flagstat ---"
samtools flagstat 12_pseudoreps/SIM/SIM_pooled_pr2.bam | head -14

# Same pair-aware approach as pooled
samtools sort -n -@ 8 07_blacklist/SIM_rep1.clean.bam \
-o 12_pseudoreps/SIM/SIM_rep1.namesorted.bam

samtools view -h 12_pseudoreps/SIM/SIM_rep1.namesorted.bam \
| awk 'BEGIN{srand(123)}
    /^@/ {
      print > "/tmp/SIM_rep1_pr1.sam"
      print > "/tmp/SIM_rep1_pr2.sam"
      next
    }
    {
      read_id = $1
      if (!(read_id in assigned)) {
        assigned[read_id] = (rand() < 0.5 ? 1 : 2)
      }
      if (assigned[read_id] == 1) print > "/tmp/SIM_rep1_pr1.sam"
      else                        print > "/tmp/SIM_rep1_pr2.sam"
    }'

samtools view -bS /tmp/SIM_rep1_pr1.sam \
| samtools sort -@ 8 -o 12_pseudoreps/SIM/SIM_rep1_pr1.bam
samtools view -bS /tmp/SIM_rep1_pr2.sam \
| samtools sort -@ 8 -o 12_pseudoreps/SIM/SIM_rep1_pr2.bam
samtools index 12_pseudoreps/SIM/SIM_rep1_pr1.bam
samtools index 12_pseudoreps/SIM/SIM_rep1_pr2.bam
rm /tmp/SIM_rep1_pr1.sam /tmp/SIM_rep1_pr2.sam

echo "--- SIM_rep1_pr1 flagstat ---"
samtools flagstat 12_pseudoreps/SIM/SIM_rep1_pr1.bam | head -14

samtools sort -n -@ 8 07_blacklist/SIM_rep2.clean.bam \
-o 12_pseudoreps/SIM/SIM_rep2.namesorted.bam

samtools view -h 12_pseudoreps/SIM/SIM_rep2.namesorted.bam \
| awk 'BEGIN{srand(456)}
    /^@/ {
      print > "/tmp/SIM_rep2_pr1.sam"
      print > "/tmp/SIM_rep2_pr2.sam"
      next
    }
    {
      read_id = $1
      if (!(read_id in assigned)) {
        assigned[read_id] = (rand() < 0.5 ? 1 : 2)
      }
      if (assigned[read_id] == 1) print > "/tmp/SIM_rep2_pr1.sam"
      else                        print > "/tmp/SIM_rep2_pr2.sam"
    }'

samtools view -bS /tmp/SIM_rep2_pr1.sam \
| samtools sort -@ 8 -o 12_pseudoreps/SIM/SIM_rep2_pr1.bam
samtools view -bS /tmp/SIM_rep2_pr2.sam \
| samtools sort -@ 8 -o 12_pseudoreps/SIM/SIM_rep2_pr2.bam
samtools index 12_pseudoreps/SIM/SIM_rep2_pr1.bam
samtools index 12_pseudoreps/SIM/SIM_rep2_pr2.bam
rm /tmp/SIM_rep2_pr1.sam /tmp/SIM_rep2_pr2.sam

echo "--- SIM_rep2_pr1 flagstat ---"
samtools flagstat 12_pseudoreps/SIM/SIM_rep2_pr1.bam | head -14

# Pooled pseudo-replicates
macs3 callpeak \
-t 12_pseudoreps/SIM/SIM_pooled_pr1.bam \
-c 07_blacklist/SIM_control.clean.bam \
-f BAMPE -g dm \
-n SIM_pooled_pr1 \
--outdir 12_pseudoreps/SIM \
--keep-dup all -p 0.01 --call-summits \
2> 12_pseudoreps/SIM/SIM_pooled_pr1_macs3.log

macs3 callpeak \
-t 12_pseudoreps/SIM/SIM_pooled_pr2.bam \
-c 07_blacklist/SIM_control.clean.bam \
-f BAMPE -g dm \
-n SIM_pooled_pr2 \
--outdir 12_pseudoreps/SIM \
--keep-dup all -p 0.01 --call-summits \
2> 12_pseudoreps/SIM/SIM_pooled_pr2_macs3.log

# Self-consistency rep1 pseudo-replicates
macs3 callpeak \
-t 12_pseudoreps/SIM/SIM_rep1_pr1.bam \
-c 07_blacklist/SIM_control.clean.bam \
-f BAMPE -g dm \
-n SIM_rep1_pr1 \
--outdir 12_pseudoreps/SIM \
--keep-dup all -p 0.01 --call-summits \
2> 12_pseudoreps/SIM/SIM_rep1_pr1_macs3.log

macs3 callpeak \
-t 12_pseudoreps/SIM/SIM_rep1_pr2.bam \
-c 07_blacklist/SIM_control.clean.bam \
-f BAMPE -g dm \
-n SIM_rep1_pr2 \
--outdir 12_pseudoreps/SIM \
--keep-dup all -p 0.01 --call-summits \
2> 12_pseudoreps/SIM/SIM_rep1_pr2_macs3.log

# Self-consistency rep2 pseudo-replicates
macs3 callpeak \
-t 12_pseudoreps/SIM/SIM_rep2_pr1.bam \
-c 07_blacklist/SIM_control.clean.bam \
-f BAMPE -g dm \
-n SIM_rep2_pr1 \
--outdir 12_pseudoreps/SIM \
--keep-dup all -p 0.01 --call-summits \
2> 12_pseudoreps/SIM/SIM_rep2_pr1_macs3.log

macs3 callpeak \
-t 12_pseudoreps/SIM/SIM_rep2_pr2.bam \
-c 07_blacklist/SIM_control.clean.bam \
-f BAMPE -g dm \
-n SIM_rep2_pr2 \
--outdir 12_pseudoreps/SIM \
--keep-dup all -p 0.01 --call-summits \
2> 12_pseudoreps/SIM/SIM_rep2_pr2_macs3.log

for NAME in SIM_pooled_pr1 SIM_pooled_pr2 \
SIM_rep1_pr1 SIM_rep1_pr2 \
SIM_rep2_pr1 SIM_rep2_pr2; do
sort -k8,8nr 12_pseudoreps/SIM/${NAME}_peaks.narrowPeak \
> 12_pseudoreps/SIM/${NAME}.sorted.narrowPeak
done

conda activate idr

# Pooled pseudo-rep IDR (Npr)
idr \
--samples \
12_pseudoreps/SIM/SIM_pooled_pr1.sorted.narrowPeak \
12_pseudoreps/SIM/SIM_pooled_pr2.sorted.narrowPeak \
--input-file-type narrowPeak \
--rank p.value --idr-threshold 0.05 \
--output-file 10_idr/SIM_pr/SIM_pooled_pr.txt \
--plot \
--log-output-file 10_idr/SIM_pr/SIM_pooled_pr.log

# Self-consistency rep1 IDR (N1)
idr \
--samples \
12_pseudoreps/SIM/SIM_rep1_pr1.sorted.narrowPeak \
12_pseudoreps/SIM/SIM_rep1_pr2.sorted.narrowPeak \
--input-file-type narrowPeak \
--rank p.value --idr-threshold 0.05 \
--output-file 10_idr/SIM_pr/SIM_rep1_selfpr.txt \
--plot \
--log-output-file 10_idr/SIM_pr/SIM_rep1_selfpr.log

# Self-consistency rep2 IDR (N2)
idr \
--samples \
12_pseudoreps/SIM/SIM_rep2_pr1.sorted.narrowPeak \
12_pseudoreps/SIM/SIM_rep2_pr2.sorted.narrowPeak \
--input-file-type narrowPeak \
--rank p.value --idr-threshold 0.05 \
--output-file 10_idr/SIM_pr/SIM_rep2_selfpr.txt \
--plot \
--log-output-file 10_idr/SIM_pr/SIM_rep2_selfpr.log

conda deactivate

NCOLS=$(awk '{print NF; exit}' 10_idr/SIM_pr/SIM_pooled_pr.txt)

if [ "$NCOLS" -ge 12 ]; then
FILTER='$12 >= 1.301'
else
  FILTER='$5 >= 540'
fi

Nt=$(awk "$FILTER"  10_idr/SIM/SIM_1v2.txt            | wc -l)
Npr=$(awk "$FILTER" 10_idr/SIM_pr/SIM_pooled_pr.txt   | wc -l)
N1=$(awk "$FILTER"  10_idr/SIM_pr/SIM_rep1_selfpr.txt | wc -l)
N2=$(awk "$FILTER"  10_idr/SIM_pr/SIM_rep2_selfpr.txt | wc -l)

RR=$(awk -v a="$Nt" -v b="$Npr" 'BEGIN{
  if (a > b) print a/b; else print b/a
}')

SCR=$(awk -v a="$N1" -v b="$N2" 'BEGIN{
  if (a > b) print a/b; else print b/a
}')

echo "=== ENCODE Pseudo-Replicate Diagnostics: SIM ==="
echo "Nt  (true replicate IDR, rep1 vs rep2): $Nt"
echo "Npr (pooled pseudo-replicate IDR):      $Npr"
echo "N1  (rep1 self-consistency IDR):        $N1"
echo "N2  (rep2 self-consistency IDR):        $N2"
echo ""
printf "Rescue Ratio (RR)           = %.3f\n" $RR
printf "Self-Consistency Ratio (SCR) = %.3f\n" $SCR
echo ""
echo "PASS:       both ratios <= 2"
echo "BORDERLINE: one ratio > 2"
echo "FAIL:       both ratios > 2"

# Save diagnostics
cat > 12_pseudoreps/SIM/SIM_pr_diagnostics.txt <<EOF
SIM ENCODE Pseudo-Replicate Diagnostics
========================================
  Nt  = $Nt
Npr = $Npr
N1  = $N1
N2  = $N2
RR  = $RR
SCR = $SCR
EOF
