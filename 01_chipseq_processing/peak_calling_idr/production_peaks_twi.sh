#!/usr/bin/env bash
# ============================================================
# production_peaks_twi.sh: peak calling, IDR, consensus, FRiP for TWI
# Single-end ChIP, 2 replicates, pairwise IDR (rep1 vs rep2)
# MACS3: -f BAM --nomodel --extsize ${TWI_EXTSIZE}, -p 1e-2
# Submit: sbatch production_peaks_twi.sh
# ============================================================

#SBATCH --job-name=peaks_twi
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --partition=<partition>
#SBATCH --output=logs/peaks_twi_%j.out
#SBATCH --error=logs/peaks_twi_%j.err
#SBATCH --no-requeue

set -euo pipefail

# ============================================================
# CONDA INITIALISATION
# ============================================================
source "$(conda info --base)/etc/profile.d/conda.sh"

# ============================================================
# CONFIGURATION
# ============================================================
BASE_DIR="$PROJECT_DIR"
BLACKLIST="$PROJECT_DIR/dm6-blacklist.v2.bed"
MACS_PVAL="1e-2"
IDR_THRESH="0.05"

# SPP fragment-size estimate for TWI: mean of rep1 (210 bp) and rep2 (230 bp),
# from 04_qc/TWI_rep{1,2}.spp_metrics.txt (column 3, first value).
TWI_EXTSIZE=220

CONDA_ENV_PEAKS="phd"
CONDA_ENV_IDR="idr"

FACTOR="TWI"
CONTROL="TWI_control.clean.bam"

# ============================================================
# DERIVED PATHS
# ============================================================
BDIR="${BASE_DIR}/07_blacklist"
PDIR="${BASE_DIR}/09_peaks/${FACTOR}"
IDIR="${BASE_DIR}/10_idr/${FACTOR}"
FDIR="${BASE_DIR}/11_final"
LOG_DIR="${BASE_DIR}/logs"

mkdir -p "${PDIR}" "${IDIR}" "${FDIR}" "${LOG_DIR}"

# ============================================================
# PRE-FLIGHT CHECKS
# ============================================================
for REP in 1 2; do
  BAM="${BDIR}/${FACTOR}_rep${REP}.clean.bam"
  [ -f "${BAM}" ] || { echo "MISSING BAM: ${BAM}"; exit 1; }
done
[ -f "${BDIR}/${CONTROL}" ] || { echo "MISSING CONTROL BAM: ${BDIR}/${CONTROL}"; exit 1; }
[ -f "${BLACKLIST}" ]       || { echo "MISSING BLACKLIST: ${BLACKLIST}"; exit 1; }

echo "========================================"
echo "Factor      : ${FACTOR}"
echo "Control     : ${CONTROL}"
echo "Library     : single-end"
echo "EXTSIZE     : ${TWI_EXTSIZE} (from SPP)"
echo "MACS p-val  : ${MACS_PVAL}"
echo "IDR thresh  : ${IDR_THRESH}"
echo "Started     : $(date)"
echo "========================================"

# ============================================================
# STAGE 1, MACS3 ON INDIVIDUAL REPLICATES (single-end)
# ============================================================
conda activate "${CONDA_ENV_PEAKS}"

for REP in 1 2; do
  echo "[$(date)] MACS3 ${FACTOR} rep${REP}..."

  macs3 callpeak \
    -t "${BDIR}/${FACTOR}_rep${REP}.clean.bam" \
    -c "${BDIR}/${CONTROL}" \
    -f BAM \
    -g dm \
    -n "${FACTOR}_rep${REP}" \
    --outdir "${PDIR}" \
    --keep-dup all \
    --nomodel --extsize "${TWI_EXTSIZE}" \
    -p "${MACS_PVAL}" \
    --call-summits \
    2> "${PDIR}/${FACTOR}_rep${REP}_macs3.log"

  sort -k8,8nr \
    "${PDIR}/${FACTOR}_rep${REP}_peaks.narrowPeak" \
    > "${PDIR}/${FACTOR}_rep${REP}.sorted.narrowPeak"

  N=$(wc -l < "${PDIR}/${FACTOR}_rep${REP}.sorted.narrowPeak")
  echo "[$(date)] ${FACTOR} rep${REP}: ${N} peaks"
done

conda deactivate

# ============================================================
# STAGE 2, IDR ON BIOLOGICAL REPLICATES
# ============================================================
conda activate "${CONDA_ENV_IDR}"

echo "[$(date)] Running IDR on biological replicates..."

idr \
  --samples \
    "${PDIR}/${FACTOR}_rep1.sorted.narrowPeak" \
    "${PDIR}/${FACTOR}_rep2.sorted.narrowPeak" \
  --input-file-type narrowPeak \
  --rank p.value \
  --idr-threshold "${IDR_THRESH}" \
  --output-file "${IDIR}/${FACTOR}_rep1_vs_rep2.txt" \
  --log-output-file "${IDIR}/${FACTOR}_rep1_vs_rep2.log" \
  --plot \
  > "${IDIR}/${FACTOR}_idr.stdout" 2>&1

conda deactivate

NCOLS=$(awk '{print NF; exit}' "${IDIR}/${FACTOR}_rep1_vs_rep2.txt")
if [ "${NCOLS}" -ge 12 ]; then
  IDR_FILTER='$12 >= 1.301'
else
  IDR_FILTER='$5 >= 540'
fi

IDR_N=$(awk "${IDR_FILTER}" "${IDIR}/${FACTOR}_rep1_vs_rep2.txt" | wc -l)
echo "[$(date)] IDR peaks at threshold ${IDR_THRESH}: ${IDR_N}"

# ============================================================
# STAGE 3, PRODUCE FINAL PEAK SET
# ============================================================
awk "${IDR_FILTER}" \
  "${IDIR}/${FACTOR}_rep1_vs_rep2.txt" \
  | sort -k7,7nr \
  > "${FDIR}/${FACTOR}_IDR${IDR_THRESH}.final.narrowPeak"

N_FINAL=$(wc -l < "${FDIR}/${FACTOR}_IDR${IDR_THRESH}.final.narrowPeak")
echo "[$(date)] Final peak set: ${N_FINAL} peaks"
echo "[$(date)] Written to:    ${FDIR}/${FACTOR}_IDR${IDR_THRESH}.final.narrowPeak"

# ============================================================
# STAGE 4, FRiP ON FINAL PEAK SET (single-end)
#
# For SE data we count one record per read (no mate flag needed).
# -F 1804 excludes unmapped/mate-unmapped/secondary/QCfail/duplicate/supplementary.
# This matches the SE FRiP convention used elsewhere for single-end factors.
# ============================================================
echo "[$(date)] Computing FRiP..."
conda activate "${CONDA_ENV_PEAKS}"

for REP in 1 2; do
  TOTAL_READS=$(samtools view -c -F 1804 \
    "${BDIR}/${FACTOR}_rep${REP}.clean.bam")

  IN_READS=$(bedtools intersect -u \
      -abam "${BDIR}/${FACTOR}_rep${REP}.clean.bam" \
      -b "${FDIR}/${FACTOR}_IDR${IDR_THRESH}.final.narrowPeak" \
    | samtools view -c -F 1804 -)

  FRIP=$(awk -v a="${IN_READS}" -v b="${TOTAL_READS}" \
    'BEGIN{ if (b>0) printf "%.4f", a/b; else printf "NA" }')

  echo "FRiP ${FACTOR}_rep${REP}: ${FRIP}  (${IN_READS} / ${TOTAL_READS} reads)"
done

conda deactivate

# ============================================================
# STAGE 5, MEDIAN PEAK WIDTH
# ============================================================
MED_WIDTH=$(awk '{print $3-$2}' \
  "${FDIR}/${FACTOR}_IDR${IDR_THRESH}.final.narrowPeak" \
  | sort -n \
  | awk '{ a[NR]=$1 }
         END {
           if (NR == 0)         { print "NA" }
           else if (NR % 2 == 1){ print a[int(NR/2)+1] }
           else                 { printf "%.0f\n", (a[NR/2]+a[NR/2+1])/2
         }}')

echo "Median peak width ${FACTOR}: ${MED_WIDTH} bp"
echo "[$(date)] === DONE: ${FACTOR} ==="