#!/bin/bash
#SBATCH --job-name=perm_cme
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=06:00:00
#SBATCH --output=logs/perm_cme_%j.out
#SBATCH --error=logs/perm_cme_%j.err
#SBATCH --nodes=1
#SBATCH --partition=<partition>
#SBATCH --get-user-env
#SBATCH --no-requeue

source $(conda info --base)/etc/profile.d/conda.sh
conda activate phd

mkdir -p logs

export PROJECT_DIR="${PROJECT_DIR:?set PROJECT_DIR to the repository root}"
cd "$PROJECT_DIR"
Rscript --no-restore --no-save 03_enrichment_and_permutation/permutation_tests/perm_tissue_tss_proximity.R 5000
