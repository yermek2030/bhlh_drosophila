#!/usr/bin/env bash
# ============================================================
# spp_cross_correlation_vvl_rib.sh
# SPP cross-correlation QC on coordinate-sorted BAMs, VVL and RIB ChIP replicates.
# PE reads treated as SE for the calculation; NSC/RSC comparable to SE data.
# Array covers 4 ChIP replicates, controls excluded.
# Submit: sbatch spp_cross_correlation_vvl_rib.sh
# ============================================================

#SBATCH --job-name=spp_vvl_rib
#SBATCH --array=0-3
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=03:00:00
#SBATCH --nodes=1
#SBATCH --partition=<partition>
#SBATCH --output=logs/spp_vvl_rib_%A_%a.out
#SBATCH --error=logs/spp_vvl_rib_%A_%a.err
#SBATCH --no-requeue

set -euo pipefail

# ============================================================
# CONFIGURATION, EDIT THESE
# ============================================================
BASE_DIR="$PROJECT_DIR"
PICARD_TMP="$PROJECT_DIR/tmp"

# ABSOLUTE path to run_spp.R inside the phantompeakqualtools install.
# Find it with:
#   find $HOME $PROJECT_DIR \
#        -name run_spp.R 2>/dev/null
SPP_SCRIPT="${SPP_SCRIPT:-run_spp.R}"   # phantompeakqualtools run_spp.R (r_phylo environment)

CONDA_ENV_R="r_phylo"   # conda env that has run_spp.R's R dependencies

SAMPLES=(
  "VVL_rep1"
  "VVL_rep2"
  "RIB_rep1"
  "RIB_rep2"
)

# ============================================================
# DERIVED PATHS
# ============================================================
SAMPLE="${SAMPLES[$SLURM_ARRAY_TASK_ID]}"
QC_DIR="${BASE_DIR}/04_qc"
LOG_DIR="${BASE_DIR}/logs"
CLEAN_BAM="${BASE_DIR}/07_blacklist/${SAMPLE}.clean.bam"

mkdir -p "${QC_DIR}" "${LOG_DIR}" "${PICARD_TMP}"

# ============================================================
# PRE-FLIGHT CHECKS
# ============================================================
[ -f "${CLEAN_BAM}" ]   || { echo "MISSING BAM: ${CLEAN_BAM}";    exit 1; }
[ -f "${SPP_SCRIPT}" ]  || { echo "MISSING run_spp.R: ${SPP_SCRIPT}"; exit 1; }

echo "========================================"
echo "Sample      : ${SAMPLE}"
echo "BAM         : ${CLEAN_BAM}"
echo "Output dir  : ${QC_DIR}"
echo "SPP script  : ${SPP_SCRIPT}"
echo "Conda env   : ${CONDA_ENV_R}"
echo "Started     : $(date)"
echo "========================================"

# ============================================================
# RUN SPP
# ============================================================
conda run -n "${CONDA_ENV_R}" Rscript "${SPP_SCRIPT}" \
  -c="${CLEAN_BAM}" \
  -savp="${QC_DIR}/${SAMPLE}.spp_plot.pdf" \
  -out="${QC_DIR}/${SAMPLE}.spp_metrics.txt" \
  -tmpdir="${PICARD_TMP}" \
  -p="${SLURM_CPUS_PER_TASK}"

echo ""
echo "[$(date)] SPP complete for ${SAMPLE}"
echo "--- ${SAMPLE} SPP metrics ---"
cat "${QC_DIR}/${SAMPLE}.spp_metrics.txt"

# ============================================================
# Post-run one-liner compiles VVL_rep1/2, RIB_rep1/2 spp_metrics.txt into vvl_rib_spp_summary.tsv; ENCODE thresholds: NSC>=1.05, RSC>=0.8, quality tag>=0
# ============================================================