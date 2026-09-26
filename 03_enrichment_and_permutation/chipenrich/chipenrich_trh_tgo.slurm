#!/bin/bash
#SBATCH --job-name=chipenrich_ch3
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH --output=logs/chipenrich_ch3_%j.out
#SBATCH --error=logs/chipenrich_ch3_%j.err
#SBATCH --nodes=1
#SBATCH --partition=<partition>
#SBATCH --get-user-env
#SBATCH --no-requeue

source $(conda info --base)/etc/profile.d/conda.sh
conda activate phd

mkdir -p logs

export PROJECT_DIR="${PROJECT_DIR:?set PROJECT_DIR to the repository root}"
cd "$PROJECT_DIR"
Rscript --no-restore --no-save 03_enrichment_and_permutation/chipenrich/chipenrich_trh_tgo.R 
