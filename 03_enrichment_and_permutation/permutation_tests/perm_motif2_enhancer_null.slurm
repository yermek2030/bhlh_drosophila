#!/bin/bash
#SBATCH --job-name=perm_motif2_enh
#SBATCH --partition=<partition>
#SBATCH --account=<account>
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --output=perm_motif2_enh_%j.log
#
# Location-matched permutation null, Trh/Tgo joint peaks vs Motif2-positive
# validated embryonic enhancers. 5000 permutations, one test. Run: sbatch perm_motif2_enhancer_null.slurm. 2 h on 32 cores.

set -eo pipefail

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate phd

export PROJECT_DIR="${PROJECT_DIR:-$(pwd)}"
echo "PROJECT_DIR=$PROJECT_DIR  cores=$SLURM_CPUS_PER_TASK"

Rscript --vanilla "$PROJECT_DIR/03_enrichment_and_permutation/permutation_tests/perm_motif2_enhancer_null.R"

echo "--- outputs ---"
ls -lh "$PROJECT_DIR/data/chapter3/perm_motif2_enh_or.rds" \
       "$PROJECT_DIR/data/chapter3/perm_motif2_enh_or_FULL.rds"
