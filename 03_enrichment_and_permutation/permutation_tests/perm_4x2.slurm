#!/bin/bash
#SBATCH --job-name=perm_4x2
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=25
#SBATCH --mem=100G
#SBATCH --time=08:00:00
#SBATCH --output=logs/perm_4x2_%j.out
#SBATCH --error=logs/perm_4x2_%j.err
#SBATCH --nodes=1
#SBATCH --partition=<partition>
#SBATCH --get-user-env
#SBATCH --no-requeue

source $(conda info --base)/etc/profile.d/conda.sh
conda activate phd

mkdir -p logs

export PROJECT_DIR="${PROJECT_DIR:?set PROJECT_DIR to the repository root}"
cd "$PROJECT_DIR"
Rscript --no-restore --no-save 03_enrichment_and_permutation/permutation_tests/perm_4x2.R 5000
