#!/usr/bin/env Rscript
# perm_4x2.R, permutation test: 4 heterodimers x 2 tissue gene lists.
# Self-contained: rebuilds every input from raw files. No RDS input.
#
# Usage:  Rscript --no-restore --no-save perm_4x2.R [ntimes]
# SLURM:  see perm_4x2.slurm

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(GenomeInfoDb)
  library(rtracklayer)
  library(regioneR)
  library(BSgenome.Dmelanogaster.UCSC.dm6)
  library(TxDb.Dmelanogaster.UCSC.dm6.ensGene)
  library(AnnotationDbi)
  library(org.Dm.eg.db)
})

# ---- Paths -------------------------------------------------------------------
## Repository root: the environment variable PROJECT_DIR, or the working directory.
PROJECT_DIR <- normalizePath(Sys.getenv("PROJECT_DIR", unset = "."), mustWork = TRUE)
if (!file.exists(file.path(PROJECT_DIR, "README.md")))
  stop("PROJECT_DIR does not point at the repository root: ", PROJECT_DIR)
REF_DIR    <- file.path(PROJECT_DIR, "data/reference_gene_lists")
FINAL_DIR  <- file.path(PROJECT_DIR, "data/final_peaks")
EXPORT_DIR    <- file.path(PROJECT_DIR, "data/chapter4")        # outputs land here
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

# ---- Helpers -----------------------------------------------------------------
load_peaks <- function(path) {
  gr <- import(path, format = "BED", extraCols = extraCols)
  gr <- gr[seqnames(gr) %in% good_chr]
  gr <- keepSeqlevels(gr, good_chr, pruning.mode = "coarse")
  gr[!duplicated(gr)]
}

# ---- Load peak sets ----------------------------------------------------------
cat("Loading peaks...\n")
Trh_clean  <- load_peaks(file.path(FINAL_DIR, "TRH_IDR0.05.final.narrowPeak"))
Tgo_clean  <- load_peaks(file.path(FINAL_DIR, "TGO_IDR0.05.final.narrowPeak"))
Sim_clean  <- load_peaks(file.path(FINAL_DIR, "SIM_IDR0.05.final.narrowPeak"))
Dys_clean  <- load_peaks(file.path(FINAL_DIR, "DYS_IDR0.05.final.narrowPeak"))
Sima_clean <- load_peaks(file.path(FINAL_DIR, "SIMA_IDR0.05.final.narrowPeak"))

# ---- Trh/Tgo joint (from Chapter 3 BED) --------------------------------------
trh_bed <- file.path(PROJECT_DIR, "data/chapter3",
                     "Trh_Tgo_joint_peaks_native_width_for_MEME_ChIP.bed")
stopifnot(file.exists(trh_bed))
trh_df <- read.table(trh_bed, header = FALSE, stringsAsFactors = FALSE,
                     col.names = c("chr","start","end","name","score","strand"))
Trh_Tgo_joint <- GRanges(
  seqnames = trh_df$chr,
  ranges   = IRanges(start = trh_df$start + 1L, end = trh_df$end),
  name     = trh_df$name
)
Trh_Tgo_joint <- keepSeqlevels(Trh_Tgo_joint, good_chr, pruning.mode = "coarse")
stopifnot(length(Trh_Tgo_joint) == 614)

# ---- Sim/Tgo, Dys/Tgo, Sima/Tgo joint peaks ----------------------------------
make_joint <- function(tf, tgo) {
  ov  <- findOverlaps(tf, tgo, minoverlap = 1, ignore.strand = TRUE)
  gr  <- pintersect(tf[queryHits(ov)], tgo[subjectHits(ov)])
  mcols(gr) <- NULL
  gr  <- keepSeqlevels(gr, good_chr, pruning.mode = "coarse")
  gr[!duplicated(gr)]
}
Sim_Tgo_joint_filt  <- make_joint(Sim_clean,  Tgo_clean)
Dys_Tgo_joint_filt  <- make_joint(Dys_clean,  Tgo_clean)
Sima_Tgo_joint_filt <- make_joint(Sima_clean, Tgo_clean)

# ---- tss_gr with FlyBase IDs -------------------------------------------------
tss_gr <- suppressMessages(promoters(genes(txdb), upstream = 0, downstream = 1))
tss_gr <- keepSeqlevels(tss_gr, good_chr, pruning.mode = "coarse")
tss_gr <- trim(tss_gr)
flybase_list <- mapIds(org.Dm.eg.db, keys = names(tss_gr),
                       keytype = "ENSEMBL", column = "FLYBASE",
                       multiVals = "list")
ambiguous   <- lengths(flybase_list) > 1
flybase_vec <- mapIds(org.Dm.eg.db, keys = names(tss_gr),
                      keytype = "ENSEMBL", column = "FLYBASE",
                      multiVals = "first")
flybase_vec[ambiguous] <- NA
tss_gr$flybase <- flybase_vec

# ---- Gene lists --------------------------------------------------------------
read_fbgn <- function(path) {
  df <- read.csv(path, stringsAsFactors = FALSE)
  v  <- if ("genes" %in% names(df)) df[["genes"]] else df[[1]]
  v  <- unique(trimws(as.character(v)))
  v[grepl("^FBgn\\d+$", v)]
}
tracheal_sg_genes <- read_fbgn(file.path(REF_DIR,
                                "trh_ts_sg_tissues_reference.csv"))
cns_midline_genes <- read_fbgn(file.path(REF_DIR,
                                "sim_cns_tissues_DEGs.csv"))
cat(sprintf("tracheal_sg=%d | cns_midline=%d\n",
            length(tracheal_sg_genes), length(cns_midline_genes)))

tss_gr$in_tracheal_sg <- !is.na(tss_gr$flybase) &
                          tss_gr$flybase %in% tracheal_sg_genes
tss_gr$in_cns_midline <- !is.na(tss_gr$flybase) &
                          tss_gr$flybase %in% cns_midline_genes
tss_tracheal <- tss_gr[tss_gr$in_tracheal_sg]
tss_cns      <- tss_gr[tss_gr$in_cns_midline]
cat(sprintf("TSSs tracheal+SG: %d | CNS+midline: %d\n",
            length(tss_tracheal), length(tss_cns)))

# ---- Permutation -------------------------------------------------------------
perm_eval_tissue <- function(A, B, ..., tissue_tss) {
  if (length(tissue_tss) == 0L) return(0L)
  hits <- findOverlaps(A, tissue_tss, maxgap = 1000L, ignore.strand = TRUE)
  length(unique(subjectHits(hits)))
}

run_tissue_permtest <- function(peaks, tissue_tss, het_label, tissue_label,
                                n = NTIMES, seed = 777) {
  set.seed(seed)
  pt <- permTest(
    A                  = peaks,
    randomize.function = randomizeRegions,
    evaluate.function  = perm_eval_tissue,
    ntimes             = n,
    per.chromosome     = TRUE,
    genome             = "dm6",
    mc.cores           = NCORES,
    tissue_tss         = tissue_tss
  )
  obs  <- pt$perm_eval_tissue$observed
  null <- pt$perm_eval_tissue$permuted
  list(
    Heterodimer = het_label,
    Tissue_list = tissue_label,
    Observed    = obs,
    Null        = null,
    Null_mean   = mean(null),
    Null_SD     = sd(null),
    Z_score     = (obs - mean(null)) / sd(null),
    p_empirical = (sum(null >= obs) + 1) / (length(null) + 1),
    N_perm      = length(null)
  )
}

cat("Running permutation tests...\n")
perm_results_raw <- list(
  run_tissue_permtest(Trh_Tgo_joint,       tss_tracheal, "Trh/Tgo",  "tracheal+SG"),
  run_tissue_permtest(Trh_Tgo_joint,       tss_cns,      "Trh/Tgo",  "CNS+midline"),
  run_tissue_permtest(Sim_Tgo_joint_filt,  tss_tracheal, "Sim/Tgo",  "tracheal+SG"),
  run_tissue_permtest(Sim_Tgo_joint_filt,  tss_cns,      "Sim/Tgo",  "CNS+midline"),
  run_tissue_permtest(Dys_Tgo_joint_filt,  tss_tracheal, "Dys/Tgo",  "tracheal+SG"),
  run_tissue_permtest(Dys_Tgo_joint_filt,  tss_cns,      "Dys/Tgo",  "CNS+midline"),
  run_tissue_permtest(Sima_Tgo_joint_filt, tss_tracheal, "Sima/Tgo", "tracheal+SG"),
  run_tissue_permtest(Sima_Tgo_joint_filt, tss_cns,      "Sima/Tgo", "CNS+midline")
)

# ---- Summary table -----------------------------------------------------------
perm_matrix <- do.call(rbind, lapply(perm_results_raw, function(r) {
  data.frame(
    Heterodimer = r$Heterodimer, Tissue_list = r$Tissue_list,
    Observed = r$Observed,
    Null_mean = round(r$Null_mean, 2), Null_SD = round(r$Null_SD, 2),
    Z_score  = round(r$Z_score, 3),
    p_empirical = r$p_empirical, N_perm = r$N_perm,
    stringsAsFactors = FALSE
  )
}))
print(perm_matrix, row.names = FALSE)

# ---- Save: full RDS (with null distributions) + CSV summary ------------------
saveRDS(list(raw = perm_results_raw, summary = perm_matrix),
        file.path(EXPORT_DIR, "perm_4x2.rds"))
write.csv(perm_matrix,
          file.path(EXPORT_DIR, "perm_4x2_summary.csv"), row.names = FALSE)
cat(sprintf("Saved: %s/perm_4x2.rds\n", EXPORT_DIR))
cat("Done.\n")
