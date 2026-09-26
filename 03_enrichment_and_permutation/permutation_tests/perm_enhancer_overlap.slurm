#!/bin/bash
#SBATCH --job-name=perm_enh_canon
#SBATCH --partition=<partition>
#SBATCH --account=<account>
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=perm_enh_canon_%j.log
#
# Enhancer overlap permutations, unsigned proximal/distal convention:
# 4 heterodimers x 3 strata = 12 tests (Ch4), plus Trh/Tgo triple
# re-exported for Ch3. 5000 permutations each.
#
#   sbatch perm_enhancer_overlap.slurm
#
# Checkpoints per heterodimer; resubmit after timeout to resume.

set -eo pipefail

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate phd

export PROJECT_DIR="${PROJECT_DIR:-$(pwd)}"
echo "PROJECT_DIR=$PROJECT_DIR  cores=$SLURM_CPUS_PER_TASK"

Rscript --vanilla "$PROJECT_DIR/03_enrichment_and_permutation/permutation_tests/perm_enhancer_overlap.R"

echo "--- outputs ---"
ls -lh "$PROJECT_DIR/data/chapter3/perm_enhancers_trh_tgo.rds" \
       "$PROJECT_DIR/data/chapter4/perm_enhancers.rds"
