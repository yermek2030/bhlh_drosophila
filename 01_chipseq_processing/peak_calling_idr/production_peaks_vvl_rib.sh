#!/usr/bin/env bash
# ============================================================
# production_peaks_vvl_rib.sh: peak calling, IDR, consensus, FRiP for VVL and RIB
# Paired-end ChIP, 2 replicates per factor, pairwise IDR rep1 vs rep2
# One job array element per factor runs MACS3, IDR, FRiP for both reps, one log
# Threshold -p 1e-2. Submit: sbatch production_peaks_vvl_rib.sh
# ============================================================

#SBATCH --job-name=peaks_vvl_rib
#SBATCH --array=0-1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --partition=<partition>
#SBATCH --output=logs/peaks_vvl_rib_%A_%a.out
#SBATCH --error=logs/peaks_vvl_rib_%A_%a.err
#SBATCH --no-requeue

set -euo pipefail

# ============================================================
# Source conda shell hook before activating envs; avoid `conda run`, it breaks pipes (e.g. bedtools | samtools view -c)
# ============================================================
source "$(conda info --base)/etc/profile.d/conda.sh"

# ============================================================
# CONFIGURATION, EDIT IF NEEDED
# ============================================================
BASE_DIR="$PROJECT_DIR"
BLACKLIST="$PROJECT_DIR/dm6-blacklist.v2.bed"
MACS_PVAL="1e-2"
IDR_THRESH="0.05"

CONDA_ENV_PEAKS="phd"   # contains macs3, samtools, bedtools
CONDA_ENV_IDR="idr"     # contains the ENCODE idr package

# ============================================================
# FACTOR TABLE
# Format: "FACTOR  CONTROL_BAM_BASENAME"
# REP1 and REP2 clean BAMs are inferred as ${FACTOR}_rep{1,2}.clean.bam
# ============================================================
FACTORS=(
  "VVL  VVL_control.clean.bam"
  "RIB  RIB_control.clean.bam"
)

# ============================================================
# DERIVED PATHS
# ============================================================
ENTRY="${FACTORS[$SLURM_ARRAY_TASK_ID]}"
FACTOR=$(echo  "${ENTRY}" | awk '{print $1}')
CONTROL=$(echo "${ENTRY}" | awk '{print $2}')

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
echo "MACS p-val  : ${MACS_PVAL}"
echo "IDR thresh  : ${IDR_THRESH}"
echo "Output peaks: ${PDIR}"
echo "Output IDR  : ${IDIR}"
echo "Output final: ${FDIR}"
echo "Started     : $(date)"
echo "========================================"

# ============================================================
# STAGE 1, MACS3 ON INDIVIDUAL REPLICATES (paired-end)
# ============================================================
conda activate "${CONDA_ENV_PEAKS}"

for REP in 1 2; do
  echo "[$(date)] MACS3 ${FACTOR} rep${REP}..."

  macs3 callpeak \
    -t "${BDIR}/${FACTOR}_rep${REP}.clean.bam" \
    -c "${BDIR}/${CONTROL}" \
    -f BAMPE \
    -g dm \
    -n "${FACTOR}_rep${REP}" \
    --outdir "${PDIR}" \
    --keep-dup all \
    -p "${MACS_PVAL}" \
    --call-summits \
    2> "${PDIR}/${FACTOR}_rep${REP}_macs3.log"

  # Sort by -log10(p-value) descending for IDR input (column 8 = -log10p)
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

# Detect column count to choose the correct IDR score column.
# Standard IDR output has column 12 = -log10(local IDR); some older builds
# truncate to 10 columns and store the score in column 5.
NCOLS=$(awk '{print NF; exit}' "${IDIR}/${FACTOR}_rep1_vs_rep2.txt")
if [ "${NCOLS}" -ge 12 ]; then
  IDR_FILTER='$12 >= 1.301'   # -log10(0.05) = 1.301
else
  IDR_FILTER='$5 >= 540'      # IDR column 5: scaled IDR score (540 corresponds to IDR = 0.05)
fi

IDR_N=$(awk "${IDR_FILTER}" "${IDIR}/${FACTOR}_rep1_vs_rep2.txt" | wc -l)
echo "[$(date)] IDR peaks at threshold ${IDR_THRESH}: ${IDR_N}"

# ============================================================
# STAGE 3, PRODUCE FINAL PEAK SET
# Filter IDR output to passing peaks, sort by signalValue (col 7) descending.
# This is the file that gets read by the R analysis chunk.
# ============================================================
awk "${IDR_FILTER}" \
  "${IDIR}/${FACTOR}_rep1_vs_rep2.txt" \
  | sort -k7,7nr \
  > "${FDIR}/${FACTOR}_IDR${IDR_THRESH}.final.narrowPeak"

N_FINAL=$(wc -l < "${FDIR}/${FACTOR}_IDR${IDR_THRESH}.final.narrowPeak")
echo "[$(date)] Final peak set: ${N_FINAL} peaks"
echo "[$(date)] Written to:    ${FDIR}/${FACTOR}_IDR${IDR_THRESH}.final.narrowPeak"

# ============================================================
# STAGE 4, FRiP ON FINAL PEAK SET
# Counts read pairs (one record per fragment) overlapping IDR peaks.
# -f 66: proper pair + first in pair, one line per fragment
# -F 1804: ENCODE exclusion mask (unmapped, mate unmapped, secondary, QC-fail, duplicate)
# ============================================================
echo "[$(date)] Computing FRiP..."
conda activate "${CONDA_ENV_PEAKS}"

for REP in 1 2; do
  TOTAL_FRAGS=$(samtools view -c -f 66 -F 1804 \
    "${BDIR}/${FACTOR}_rep${REP}.clean.bam")

  IN_FRAGS=$(bedtools intersect -u \
      -abam "${BDIR}/${FACTOR}_rep${REP}.clean.bam" \
      -b "${FDIR}/${FACTOR}_IDR${IDR_THRESH}.final.narrowPeak" \
    | samtools view -c -f 66 -F 1804 -)

  FRIP=$(awk -v a="${IN_FRAGS}" -v b="${TOTAL_FRAGS}" \
    'BEGIN{ if (b>0) printf "%.4f", a/b; else printf "NA" }')

  echo "FRiP ${FACTOR}_rep${REP}: ${FRIP}  (${IN_FRAGS} / ${TOTAL_FRAGS} fragments)"
done

conda deactivate

# ============================================================
# STAGE 5, MEDIAN PEAK WIDTH (for QC table)
# ============================================================
MED_WIDTH=$(awk '{print $3-$2}' \
  "${FDIR}/${FACTOR}_IDR${IDR_THRESH}.final.narrowPeak" \
  | sort -n \
  | awk '{ a[NR]=$1 }
         END {
           if (NR == 0)         { print "NA" }
           else if (NR % 2 == 1){ print a[int(NR/2)+1] }
           else                 { printf "%.0f\n", (a[NR/2]+a[NR/2+1])/2 }
         }')

echo "Median peak width ${FACTOR}: ${MED_WIDTH} bp"
echo "[$(date)] === DONE: ${FACTOR} ==="