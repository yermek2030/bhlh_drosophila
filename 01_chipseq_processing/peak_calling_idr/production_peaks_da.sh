#!/usr/bin/env bash
# ============================================================
# production_peaks_da.sh
# da ChIP-seq: download -> align -> dedup -> blacklist -> SPP -> MACS3 -> IDR -> peak set -> FRiP
# ENCSR416MJQ target da (da-eGFP, anti-GFP): ENCFF303EVV rep2, ENCFF831FKP rep3, single-end 50bp
# Input ENCSR515BLI: ENCFF405BXW, single-end
#SBATCH --job-name=da_full
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=08:00:00
#SBATCH --nodes=1
#SBATCH --partition=<partition>
#SBATCH --output=logs/peaks_da_%j.out
#SBATCH --error=logs/peaks_da_%j.err
#SBATCH --no-requeue

set -euo pipefail
source "$(conda info --base)/etc/profile.d/conda.sh"

BASE="$PROJECT_DIR"
THREADS="${SLURM_CPUS_PER_TASK:-32}"
BLACKLIST="${BASE}/dm6-blacklist.v2.bed"
IDX="${BASE}/dm6"
MACS_PVAL="1e-2"
IDR_THRESH="0.05"
FACTOR="DA"
# run_spp.R from phantompeakqualtools 1.2.2; set RUN_SPP to its path if it is not on PATH.

FQ="${BASE}/01_fastq"; BAMD="${BASE}/03_bam"; QC="${BASE}/04_qc"
BDIR="${BASE}/07_blacklist"; PDIR="${BASE}/09_peaks/${FACTOR}"
IDIR="${BASE}/10_idr/${FACTOR}"; FDIR="${BASE}/11_final"; LOGD="${BASE}/logs"
mkdir -p "$FQ" "$BAMD" "$QC" "$BDIR" "$PDIR" "$IDIR" "$FDIR" "$LOGD"

declare -A ACC=( [rep1]=ENCFF303EVV [rep2]=ENCFF831FKP [control]=ENCFF405BXW )

echo "=========================================="
echo "Factor     : ${FACTOR} (Daughterless)"
echo "Reps       : 2 true biological (ENCODE bio reps 2 and 3)"
echo "Control    : ENCFF405BXW (ENCSR515BLI, input, unreplicated)"
echo "Library    : single-end 50 bp"
echo "Index      : ${IDX} (existing dm6, UCSC naming)"
echo "MACS p-val : ${MACS_PVAL}   IDR: ${IDR_THRESH}"
echo "Threads    : ${THREADS}"
echo "Started    : $(date)"
echo "=========================================="

# ---------------------------------------------------------------- STAGE 0: get
for K in rep1 rep2 control; do
  A="${ACC[$K]}"; DEST="${FQ}/${A}.fastq.gz"
  if [ -s "$DEST" ]; then echo "[$(date)] $A already present, skipping"; continue; fi
  echo "[$(date)] downloading $A ($K)..."
  curl -fsSL --retry 3 --retry-delay 5 \
    "https://www.encodeproject.org/files/${A}/@@download/${A}.fastq.gz" -o "${DEST}.part"
  gzip -t "${DEST}.part" || { echo "CORRUPT DOWNLOAD: $A"; exit 1; }
  mv "${DEST}.part" "$DEST"
  echo "[$(date)] $A  $(du -h "$DEST" | cut -f1)  reads=$(( $(zcat "$DEST" | wc -l) / 4 ))"
done

conda activate phd

# ------------------------------------------------- STAGE 1: align, dedup, clean
for K in rep1 rep2 control; do
  A="${ACC[$K]}"
  case "$K" in control) TAG="${FACTOR}_control";; *) TAG="${FACTOR}_${K}";; esac
  CLEAN="${BDIR}/${TAG}.clean.bam"
  if [ -s "$CLEAN" ]; then echo "[$(date)] $TAG clean BAM exists, skipping"; continue; fi

  echo "[$(date)] bowtie2 ${TAG}..."
  bowtie2 -x "$IDX" -U "${FQ}/${A}.fastq.gz" -p "$THREADS" \
    2> "${QC}/${TAG}.bowtie2.log" \
  | samtools view -b -q 30 -F 1804 -@ 4 - \
  | samtools sort -@ "$THREADS" -o "${BAMD}/${TAG}.sorted.bam" -
  samtools index "${BAMD}/${TAG}.sorted.bam"

  # picard's conda launcher defaults to a 2 GB heap. rep1 (8.3 M reads) fits and
  # passes; a larger library exceeds the in-memory SortingCollection, spills to
  # disk, and the spill path loads the Snappy native library, which fails on this
  # host (UnsatisfiedLinkError in SnappyLoader). Two independent guards: a heap
  # large enough that no spill happens, and Snappy disabled so that if one does
  # happen it falls back to the JDK deflater. TMP_DIR is pinned to the project scratch rather
  # than the node's 70 GB root filesystem.
  echo "[$(date)] picard MarkDuplicates ${TAG}..."
  mkdir -p "${BAMD}/tmp_${TAG}"
  picard -Xmx32g -Dsamjdk.snappy.disable=true MarkDuplicates \
    -I "${BAMD}/${TAG}.sorted.bam" \
    -O "${BAMD}/${TAG}.dedup.bam" -M "${QC}/${TAG}.dupmetrics.txt" \
    --REMOVE_DUPLICATES true --VALIDATION_STRINGENCY LENIENT \
    --TMP_DIR "${BAMD}/tmp_${TAG}" --MAX_RECORDS_IN_RAM 8000000 \
    2> "${QC}/${TAG}.picard.log"
  rm -rf "${BAMD}/tmp_${TAG}"
  samtools index "${BAMD}/${TAG}.dedup.bam"

  echo "[$(date)] blacklist filter ${TAG}..."
  bedtools intersect -v -abam "${BAMD}/${TAG}.dedup.bam" -b "$BLACKLIST" > "$CLEAN"
  samtools index "$CLEAN"
  echo "[$(date)] ${TAG}: $(samtools view -c "$CLEAN") reads retained"
done

# --------------------------------------------- STAGE 2: SPP fragment length
EXT_VALS=()
for REP in rep1 rep2; do
  TAG="${FACTOR}_${REP}"; OUT="${QC}/${TAG}.spp_metrics.txt"
  if [ ! -s "$OUT" ]; then
    echo "[$(date)] SPP cross-correlation ${TAG}..."
    Rscript "${RUN_SPP:-run_spp.R}" \
      -c="${BDIR}/${TAG}.clean.bam" -savp -out="$OUT" -p="$THREADS" \
      > "${QC}/${TAG}.spp.log" 2>&1 || echo "SPP failed for ${TAG}, see log"
  fi
  V=$(awk -F'\t' 'NR==1{split($3,a,","); print a[1]}' "$OUT" 2>/dev/null || echo "")
  echo "[$(date)] ${TAG} SPP fragment length: ${V:-UNAVAILABLE}"
  [ -n "${V:-}" ] && EXT_VALS+=("$V")
done

if [ "${#EXT_VALS[@]}" -eq 2 ]; then
  D=$(( EXT_VALS[0] > EXT_VALS[1] ? EXT_VALS[0]-EXT_VALS[1] : EXT_VALS[1]-EXT_VALS[0] ))
  EXTSIZE=$(( (EXT_VALS[0] + EXT_VALS[1]) / 2 ))
  echo "[$(date)] SPP estimates ${EXT_VALS[0]} and ${EXT_VALS[1]}, difference ${D} bp"
  if [ "$D" -gt 30 ]; then
    echo "WARNING: replicate SPP estimates differ by more than the 30 bp the project"
    echo "         convention allows for a single shared --extsize. Using the mean"
    echo "         ${EXTSIZE}; per-replicate values may be preferable."
  fi
else
  echo "ERROR: SPP did not yield two fragment lengths. Stopping before MACS3 rather"
  echo "       than substituting a guessed --extsize, which would break"
  echo "       cross-comparability with TWI/VVL/RIB."
  exit 1
fi
echo "[$(date)] EXTSIZE = ${EXTSIZE}"

# --------------------------------------------------------- STAGE 3: MACS3
for REP in rep1 rep2; do
  echo "[$(date)] MACS3 ${FACTOR} ${REP}..."
  macs3 callpeak \
    -t "${BDIR}/${FACTOR}_${REP}.clean.bam" \
    -c "${BDIR}/${FACTOR}_control.clean.bam" \
    -f BAM -g dm -n "${FACTOR}_${REP}" --outdir "$PDIR" \
    --keep-dup all --nomodel --extsize "$EXTSIZE" \
    -p "$MACS_PVAL" --call-summits \
    2> "${PDIR}/${FACTOR}_${REP}_macs3.log"
  sort -k8,8nr "${PDIR}/${FACTOR}_${REP}_peaks.narrowPeak" \
    > "${PDIR}/${FACTOR}_${REP}.sorted.narrowPeak"
  echo "[$(date)] ${FACTOR} ${REP}: $(wc -l < "${PDIR}/${FACTOR}_${REP}.sorted.narrowPeak") peaks"
done
conda deactivate

# ----------------------------------------------------------- STAGE 4: IDR
conda activate idr
echo "[$(date)] IDR rep1 vs rep2..."
idr --samples "${PDIR}/${FACTOR}_rep1.sorted.narrowPeak" \
              "${PDIR}/${FACTOR}_rep2.sorted.narrowPeak" \
    --input-file-type narrowPeak --rank p.value \
    --idr-threshold "$IDR_THRESH" \
    --output-file "${IDIR}/${FACTOR}_rep1_vs_rep2.txt" \
    --log-output-file "${IDIR}/${FACTOR}_rep1_vs_rep2.log" \
    --plot > "${IDIR}/${FACTOR}_idr.stdout" 2>&1
conda deactivate

NCOLS=$(awk '{print NF; exit}' "${IDIR}/${FACTOR}_rep1_vs_rep2.txt")
if [ "$NCOLS" -ge 12 ]; then IDR_FILTER='$12 >= 1.301'; else IDR_FILTER='$5 >= 540'; fi
awk "$IDR_FILTER" "${IDIR}/${FACTOR}_rep1_vs_rep2.txt" | sort -k7,7nr \
  > "${FDIR}/${FACTOR}_IDR${IDR_THRESH}.final.narrowPeak"
N_FINAL=$(wc -l < "${FDIR}/${FACTOR}_IDR${IDR_THRESH}.final.narrowPeak")
echo "[$(date)] FINAL: ${N_FINAL} peaks -> ${FDIR}/${FACTOR}_IDR${IDR_THRESH}.final.narrowPeak"
echo "  (for scale: TWI 4115, TRH 1033, VVL 854)"

# ---------------------------------------------------------- STAGE 5: FRiP
conda activate phd
for REP in rep1 rep2; do
  TOT=$(samtools view -c "${BDIR}/${FACTOR}_${REP}.clean.bam")
  INP=$(bedtools intersect -u -a "${BDIR}/${FACTOR}_${REP}.clean.bam" \
        -b "${FDIR}/${FACTOR}_IDR${IDR_THRESH}.final.narrowPeak" -ubam | samtools view -c)
  echo "[$(date)] FRiP ${REP}: $(awk -v a="$INP" -v b="$TOT" 'BEGIN{printf "%.4f", a/b}') (${INP}/${TOT})"
done
conda deactivate
echo "[$(date)] DONE"
