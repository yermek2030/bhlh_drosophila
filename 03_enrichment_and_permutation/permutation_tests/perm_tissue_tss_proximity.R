#!/usr/bin/env Rscript
# perm_tissue_tss_proximity.R, HPC permutation test for Trh/Tgo proximity to TS/SG TSSs.
# Single-core forced (mc.cores = 1) because multi-core mclapply triggers regioneR
# serialisation failures.
#
# Usage:  Rscript --no-restore --no-save perm_tissue_tss_proximity.R [ntimes]
# SLURM:  see perm_tissue_tss_proximity.slurm   (1 CPU / 8 GB)

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
REF_DIR   <- file.path(PROJECT_DIR, "data/reference_gene_lists")
EXPORT_DIR   <- file.path(PROJECT_DIR, "data/chapter3")
dir.create(EXPORT_DIR, showWarnings = FALSE, recursive = TRUE)

NTIMES <- as.integer(commandArgs(trailingOnly = TRUE)[1])
if (is.na(NTIMES)) NTIMES <- 5000L
cat(sprintf("ntimes=%d | cores=1 (forced)\n", NTIMES))
options(mc.cores = 1)

good_chr <- c("chr2L", "chr2R", "chr3L", "chr3R", "chr4", "chrX")
txdb     <- TxDb.Dmelanogaster.UCSC.dm6.ensGene

# ---- joint_peaks_filt -------------------------------------------------------
trh_bed <- file.path(PROJECT_DIR, "data/chapter3",
                     "Trh_Tgo_joint_peaks_native_width_for_MEME_ChIP.bed")
stopifnot(file.exists(trh_bed))
trh_df <- read.table(trh_bed, header = FALSE, stringsAsFactors = FALSE,
                     col.names = c("chr","start","end","name","score","strand"))
joint_peaks_filt <- GRanges(
  seqnames = trh_df$chr,
  ranges   = IRanges(start = trh_df$start + 1L, end = trh_df$end)
)
joint_peaks_filt <- keepSeqlevels(joint_peaks_filt, good_chr,
                                  pruning.mode = "coarse")
stopifnot(length(joint_peaks_filt) == 614)

# ---- tss_gr with FlyBase IDs ------------------------------------------------
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

# ---- TS/SG gene list ---------------------------------------------------------
read_fbgn <- function(path) {
  df <- read.csv(path, stringsAsFactors = FALSE)
  v  <- if ("genes" %in% names(df)) df[["genes"]] else df[[1]]
  v  <- unique(trimws(as.character(v)))
  v[grepl("^FBgn\\d+$", v)]
}
ts_sg_genes <- read_fbgn(file.path(REF_DIR, "trh_ts_sg_tissues_reference.csv"))
cat(sprintf("ts_sg_genes: %d\n", length(ts_sg_genes)))

tss_gr$is_ts_sg <- !is.na(tss_gr$flybase) & tss_gr$flybase %in% ts_sg_genes
tss_ts_sg <- tss_gr[tss_gr$is_ts_sg]
cat(sprintf("TS/SG TSSs for permutation: %d\n", length(tss_ts_sg)))

# ---- Permutation -------------------------------------------------------------
perm_eval <- function(A, B, ..., tss_ts_sg) {
  if (length(tss_ts_sg) == 0L) return(0L)
  hits <- findOverlaps(A, tss_ts_sg, maxgap = 1000L, ignore.strand = TRUE)
  length(unique(subjectHits(hits)))
}

set.seed(777)
pt <- permTest(
  A                  = joint_peaks_filt,
  randomize.function = randomizeRegions,
  evaluate.function  = perm_eval,
  ntimes             = NTIMES,
  per.chromosome     = TRUE,
  genome             = "dm6",
  tss_ts_sg          = tss_ts_sg
)

# ---- Summary -----------------------------------------------------------------
obs  <- pt$perm_eval$observed
null <- pt$perm_eval$permuted
tbl_perm <- data.frame(
  Context        = "Trh/Tgo peaks within 1 kb of TS/SG TSS",
  Observed       = obs,
  Null_mean      = round(mean(null), 2),
  Null_SD        = round(sd(null),   2),
  Z_score        = round((obs - mean(null)) / sd(null), 3),
  p_empirical    = (sum(null >= obs) + 1) / (length(null) + 1),  # unbiased
  N_permutations = length(null),
  stringsAsFactors = FALSE
)
print(tbl_perm, row.names = FALSE)

# ---- Save --------------------------------------------------------------------
saveRDS(list(pt = pt, summary = tbl_perm),
        file.path(EXPORT_DIR, "perm_cme.rds"))
write.csv(tbl_perm,
          file.path(EXPORT_DIR, "perm_cme_summary.csv"), row.names = FALSE)
cat(sprintf("Saved: %s/perm_cme.rds\n", EXPORT_DIR))
cat("Done.\n")
