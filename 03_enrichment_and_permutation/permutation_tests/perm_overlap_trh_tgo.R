#!/usr/bin/env Rscript
# perm_overlap_trh_tgo.R, Trh/Tgo overlap significance test (peakPermTest).
# Self-contained: rebuilds Trh_ChIP_peaks and Tgo_ChIP_peaks from raw files.
#
# Usage:  Rscript --no-restore --no-save perm_overlap_trh_tgo.R [ntimes]
# SLURM:  see perm_overlap_trh_tgo.slurm

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
EXPORT_DIR   <- file.path(PROJECT_DIR, "data/chapter3")
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

cat("Loading Trh and Tgo peaks...\n")
Trh_ChIP_peaks <- load_peaks(file.path(FINAL_DIR, "TRH_IDR0.05.final.narrowPeak"))
Tgo_ChIP_peaks <- load_peaks(file.path(FINAL_DIR, "TGO_IDR0.05.final.narrowPeak"))
cat(sprintf("Trh=%d  Tgo=%d\n", length(Trh_ChIP_peaks), length(Tgo_ChIP_peaks)))

# ---- peakPermTest ------------------------------------------------------------
set.seed(777)
pt2_ <- peakPermTest(
  Trh_ChIP_peaks, Tgo_ChIP_peaks,
  ntimes         = NTIMES,
  mc.cores       = NCORES,
  TxDb           = txdb,
  seed           = 777,
  force.parallel = TRUE
)

# ---- Summary -----------------------------------------------------------------
perm_counts <- pt2_$cntOverlaps$permuted
obs_count   <- pt2_$cntOverlaps$observed
pval        <- pt2_$cntOverlaps$pval
zscore      <- pt2_$cntOverlaps$zscore
ev_perm_95  <- as.numeric(quantile(perm_counts, 0.95))

cat(sprintf("\nObserved overlap : %d\n", obs_count))
cat(sprintf("Null mean        : %.2f\n", mean(perm_counts)))
cat(sprintf("Null SD          : %.2f\n", sd(perm_counts)))
cat(sprintf("Z-score          : %.2f\n", zscore))
cat(sprintf("Empirical p      : %.3e\n", pval))
cat(sprintf("95th-pct of null : %.0f\n", ev_perm_95))

summary_df <- data.frame(
  Observed       = obs_count,
  Null_mean      = round(mean(perm_counts), 2),
  Null_SD        = round(sd(perm_counts), 2),
  Z_score        = round(zscore, 2),
  p_empirical    = pval,
  Ev_perm_95     = round(ev_perm_95),
  N_permutations = length(perm_counts),
  stringsAsFactors = FALSE
)

# ---- Save --------------------------------------------------------------------
saveRDS(list(pt2_ = pt2_, summary = summary_df),
        file.path(EXPORT_DIR, "perm_overlap_trh_tgo.rds"))
write.csv(summary_df,
          file.path(EXPORT_DIR, "perm_overlap_trh_tgo_summary.csv"),
          row.names = FALSE)
cat(sprintf("\nSaved: %s/perm_overlap_trh_tgo.rds\n", EXPORT_DIR))
cat("Done.\n")
