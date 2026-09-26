#!/usr/bin/env Rscript
# ============================================================
# chipenrich_trh_tgo_distal_split.R
# Promoter-proximal vs distal ChIP-Enrich on Trh/Tgo joint peaks (614).
# Checks whether cytoplasmic-translation GO signal is a promoter-proximal locus-density artefact.
# Proximal (<=1 kb): locusdef "1kb"; distal (>1 kb): locusdef "nearest_tss".
# Writes: data/chapter3/chipenrich_ch3_distalsplit.rds
# ============================================================

## Repository root: the environment variable PROJECT_DIR, or the working directory.
PROJECT_DIR <- normalizePath(Sys.getenv("PROJECT_DIR", unset = "."), mustWork = TRUE)
if (!file.exists(file.path(PROJECT_DIR, "README.md")))
  stop("PROJECT_DIR does not point at the repository root: ", PROJECT_DIR)
EXPORT_DIR <- file.path(PROJECT_DIR, "data/chapter3")
stopifnot(dir.exists(EXPORT_DIR))
suppressPackageStartupMessages({
  library(chipenrich)
  library(GenomicRanges)
})
# --- Rebuild the canonical 614-peak set with distanceToTSS -------------------
# This script is standalone, so it reconstructs joint_peaks_filt by loading the
# saved peak object (joint_peaks_filt.rds) rather than repeating the peak-build
# steps here.
jp_path <- file.path(EXPORT_DIR, "joint_peaks_filt.rds")
if (!file.exists(jp_path)) {
  stop("joint_peaks_filt.rds not found in chapter3/. ",
       "Add a saveRDS(joint_peaks_filt, jp_path) line to the Rmd peak-build ",
       "chunk and knit once, or point jp_path at your existing object.")
}
joint_peaks_filt <- readRDS(jp_path)
stopifnot(length(joint_peaks_filt) == 614)
stopifnot("distanceToTSS" %in% names(mcols(joint_peaks_filt)))

# --- Split (identical thresholds to Rmd lines 849-850) -----------------------
promoter_peaks <- joint_peaks_filt[joint_peaks_filt$distanceToTSS <= 1000]
distal_peaks   <- joint_peaks_filt[joint_peaks_filt$distanceToTSS >  1000]
cat("Promoter peaks (<=1 kb):", length(promoter_peaks), "\n")
cat("Distal peaks   (>1 kb):",  length(distal_peaks),   "\n")
stopifnot(length(promoter_peaks) + length(distal_peaks) == 614)
to_bed <- function(gr) data.frame(
  chrom = as.character(seqnames(gr)),  # already chr-prefixed UCSC
  start = start(gr), end = end(gr),
  stringsAsFactors = FALSE
)
run_bp <- function(bed, ldef) {
  chipenrich(peaks = bed, genome = "dm6", genesets = "GOBP",
             locusdef = ldef, qc_plots = FALSE, out_name = NULL, n_cores = 1)$results
}
res_prox_bp   <- run_bp(to_bed(promoter_peaks), "nearest_tss")
res_distal_bp <- run_bp(to_bed(distal_peaks),   "nearest_tss")
# Significant subsets (FDR < 0.05), ordered
sig <- function(df) df[df$FDR < 0.05, ][order(df[df$FDR < 0.05, ]$FDR), ]
res_prox_bp_sig   <- sig(res_prox_bp)
res_distal_bp_sig <- sig(res_distal_bp)

cat("\n== PROXIMAL top terms ==\n")
print(head(res_prox_bp_sig[, c("Description","Odds.Ratio","P.value","FDR","N.Geneset.Peak.Genes")], 15))
cat("\n== DISTAL top terms ==\n")
print(head(res_distal_bp_sig[, c("Description","Odds.Ratio","P.value","FDR","N.Geneset.Peak.Genes")], 15))

saveRDS(list(
  res_prox_bp        = res_prox_bp,
  res_distal_bp      = res_distal_bp,
  res_prox_bp_sig    = res_prox_bp_sig,
  res_distal_bp_sig  = res_distal_bp_sig,
  n_promoter         = length(promoter_peaks),
  n_distal           = length(distal_peaks)
), file.path(EXPORT_DIR, "chipenrich_ch3_distalsplit.rds"))
cat("\nSaved: chipenrich_ch3_distalsplit.rds\n")