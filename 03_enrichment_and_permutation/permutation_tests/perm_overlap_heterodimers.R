#!/usr/bin/env Rscript
# perm_overlap_heterodimers.R, permutation test for heterodimer overlap significance.
# Tests whether Sim/Tgo, Dys/Tgo, Sima/Tgo overlaps with Tgo are above chance.
# Uses the Chapter 3 approach: peakPermTest, ntimes=5000, seed=777.
#
# Usage:  Rscript --no-restore --no-save perm_overlap_heterodimers.R [ntimes]
# SLURM:  see perm_overlap_heterodimers.slurm

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(GenomeInfoDb)
  library(rtracklayer)
  library(ChIPpeakAnno)
  library(TxDb.Dmelanogaster.UCSC.dm6.ensGene)
  library(BSgenome.Dmelanogaster.UCSC.dm6)
})

# ---- Paths -------------------------------------------------------------------
## Repository root: the environment variable PROJECT_DIR, or the working directory.
PROJECT_DIR <- normalizePath(Sys.getenv("PROJECT_DIR", unset = "."), mustWork = TRUE)
if (!file.exists(file.path(PROJECT_DIR, "README.md")))
  stop("PROJECT_DIR does not point at the repository root: ", PROJECT_DIR)
FINAL_DIR <- file.path(PROJECT_DIR, "data/final_peaks")
EXPORT_DIR   <- file.path(PROJECT_DIR, "data/chapter4")
dir.create(EXPORT_DIR, showWarnings = FALSE, recursive = TRUE)

NTIMES <- as.integer(commandArgs(trailingOnly = TRUE)[1])
if (is.na(NTIMES)) NTIMES <- 5000L
NCORES <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "1"))
cat(sprintf("ntimes=%d | cores=%d\n", NTIMES, NCORES))

# ---- Constants ---------------------------------------------------------------
good_chr  <- c("chr2L", "chr2R", "chr3L", "chr3R", "chr4", "chrX")
txdb      <- TxDb.Dmelanogaster.UCSC.dm6.ensGene
extraCols <- c(signalValue = "numeric", pValue = "numeric",
               qValue = "numeric", peak = "integer")

# ---- Load peak sets ----------------------------------------------------------
load_peaks <- function(path) {
  gr <- import(path, format = "BED", extraCols = extraCols)
  gr <- gr[seqnames(gr) %in% good_chr]
  gr <- keepSeqlevels(gr, good_chr, pruning.mode = "coarse")
  gr[!duplicated(gr)]
}

cat("Loading peaks...\n")
Tgo_clean  <- load_peaks(file.path(FINAL_DIR, "TGO_IDR0.05.final.narrowPeak"))
Sim_clean  <- load_peaks(file.path(FINAL_DIR, "SIM_IDR0.05.final.narrowPeak"))
Dys_clean  <- load_peaks(file.path(FINAL_DIR, "DYS_IDR0.05.final.narrowPeak"))
Sima_clean <- load_peaks(file.path(FINAL_DIR, "SIMA_IDR0.05.final.narrowPeak"))
cat(sprintf("Tgo=%d Sim=%d Dys=%d Sima=%d\n",
            length(Tgo_clean), length(Sim_clean),
            length(Dys_clean), length(Sima_clean)))

# ---- Permutation -------------------------------------------------------------
run_peakPermTest <- function(tf_peaks, tgo_peaks, label, seed = 777) {
  set.seed(seed)
  pt <- peakPermTest(
    tf_peaks, tgo_peaks,
    ntimes         = NTIMES,
    mc.cores       = NCORES,
    TxDb           = txdb,
    seed           = seed,
    force.parallel = TRUE
  )
  list(
    label       = label,
    observed    = pt$cntOverlaps$observed,
    perm_counts = pt$cntOverlaps$permuted,
    pval        = pt$cntOverlaps$pval,
    zscore      = pt$cntOverlaps$zscore
  )
}

cat("Running permutation tests...\n")
pt_sim_tgo  <- run_peakPermTest(Sim_clean,  Tgo_clean, "Sim/Tgo")
pt_dys_tgo  <- run_peakPermTest(Dys_clean,  Tgo_clean, "Dys/Tgo")
pt_sima_tgo <- run_peakPermTest(Sima_clean, Tgo_clean, "Sima/Tgo")

# ---- Summary table -----------------------------------------------------------
perm_summary <- data.frame(
  Pair      = c("Sim/Tgo", "Dys/Tgo", "Sima/Tgo"),
  Observed  = c(pt_sim_tgo$observed,  pt_dys_tgo$observed,  pt_sima_tgo$observed),
  Null_mean = c(mean(pt_sim_tgo$perm_counts),
                mean(pt_dys_tgo$perm_counts),
                mean(pt_sima_tgo$perm_counts)),
  Null_SD   = c(sd(pt_sim_tgo$perm_counts),
                sd(pt_dys_tgo$perm_counts),
                sd(pt_sima_tgo$perm_counts)),
  Z_score   = c(pt_sim_tgo$zscore, pt_dys_tgo$zscore, pt_sima_tgo$zscore),
  p_value   = c(pt_sim_tgo$pval,   pt_dys_tgo$pval,   pt_sima_tgo$pval),
  N_perm    = NTIMES,
  stringsAsFactors = FALSE
)
print(perm_summary, row.names = FALSE)

# ---- Save --------------------------------------------------------------------
saveRDS(list(
  pt_sim_tgo  = pt_sim_tgo,
  pt_dys_tgo  = pt_dys_tgo,
  pt_sima_tgo = pt_sima_tgo,
  summary     = perm_summary
), file.path(EXPORT_DIR, "perm_overlaps.rds"))

write.csv(perm_summary,
          file.path(EXPORT_DIR, "perm_overlaps_summary.csv"), row.names = FALSE)
cat(sprintf("Saved: %s/perm_overlaps.rds\n", EXPORT_DIR))
cat("Done.\n")
