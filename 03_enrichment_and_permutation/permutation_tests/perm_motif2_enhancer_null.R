#!/usr/bin/env Rscript
# perm_motif2_enhancer_null.R: location-matched permutation null for
# co-localisation of 614 Trh/Tgo joint peaks with Motif2-positive
# enhancers, peaks scattered at random across dm6.
# Output: observed = 109, z = 4.48, p = 2e-04 (5000 perms), n_peaks = 614.
# Usage: Rscript --vanilla perm_motif2_enhancer_null.R [ntimes]
# SLURM: sbatch perm_motif2_enhancer_null.slurm

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(GenomeInfoDb)
  library(regioneR)
})

# ---- Paths -------------------------------------------------------------------
# PROJECT_DIR is set as an environment variable; if not set, defaults to the
# current directory so the script can be run locally (observed value only).
## Repository root: the environment variable PROJECT_DIR, or the working directory.
PROJECT_DIR <- normalizePath(Sys.getenv("PROJECT_DIR", unset = "."), mustWork = TRUE)
if (!file.exists(file.path(PROJECT_DIR, "README.md")))
  stop("PROJECT_DIR does not point at the repository root: ", PROJECT_DIR)
if (!nzchar(PROJECT_DIR)) PROJECT_DIR <- PROJECT_DIR

AUX_DIR    <- file.path(PROJECT_DIR, "data/chapter3")
REF_DIR    <- file.path(PROJECT_DIR, "data/reference_gene_lists")
EXPORT_DIR <- AUX_DIR
dir.create(EXPORT_DIR, showWarnings = FALSE, recursive = TRUE)

NTIMES <- as.integer(commandArgs(trailingOnly = TRUE)[1])
if (is.na(NTIMES)) NTIMES <- 5000L
NCORES <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "1"))
SEED   <- 777L
cat(sprintf("PROJECT_DIR=%s\nntimes=%d | cores=%d | seed=%d\n",
            PROJECT_DIR, NTIMES, NCORES, SEED))

good_chr <- c("chr2L", "chr2R", "chr3L", "chr3R", "chr4", "chrX")

# ---- 1. Trh/Tgo joint peaks (614) -------------------------------------------
# Same BED as perm_enhancer_overlap.R, so the peak set is byte-identical to
# the one used for the enhancer-stratum permutations and for MEME-ChIP.
trh_bed <- file.path(AUX_DIR, "Trh_Tgo_joint_peaks_native_width_for_MEME_ChIP.bed")
stopifnot(file.exists(trh_bed))
trh_df <- read.table(trh_bed, header = FALSE, stringsAsFactors = FALSE,
                     col.names = c("chr", "start", "end", "name", "score", "strand"))
joint_peaks_filt <- GRanges(
  seqnames = trh_df$chr,
  ranges   = IRanges(start = trh_df$start + 1L, end = trh_df$end),  # BED 0-based -> 1-based
  name     = trh_df$name
)
joint_peaks_filt <- keepSeqlevels(joint_peaks_filt, good_chr, pruning.mode = "coarse")
stopifnot(length(joint_peaks_filt) == 614L)

# ---- 2. Annotated embryonic enhancers ---------------------------------------
enh_bed <- file.path(REF_DIR, "merged_embryo_5-20_enhancers_dm6.bed")
stopifnot(file.exists(enh_bed))
enh_df <- read.table(enh_bed, header = FALSE, stringsAsFactors = FALSE)
colnames(enh_df)[1:3] <- c("seqnames", "start", "end")
enhancers_filt <- GRanges(
  seqnames = enh_df$seqnames,
  ranges   = IRanges(start = enh_df$start + 1L, end = enh_df$end)
)
enhancers_filt <- keepSeqlevels(enhancers_filt, good_chr, pruning.mode = "coarse")

# fimo.tsv: genome-wide scan of Motif2 PWM (STRCACTGARARAAA.meme) vs dm6; sequence_name/start/stop already chromosomal, no conversion needed
# threshold p < 1e-4 (FIMO default, not q <= 0.05)
# gives 109, recorded in results/reported_values.tsv, Chapter 3
fimo_tsv <- file.path(AUX_DIR, "fimo_enhancers_genomewide/fimo.tsv")
stopifnot(file.exists(fimo_tsv))
fimo <- read.delim(fimo_tsv, comment.char = "#", stringsAsFactors = FALSE)
fimo <- fimo[!is.na(fimo$start) & nzchar(fimo$sequence_name), ]
fimo <- fimo[fimo$motif_id == "STRCACTGARARAAA", ]   # explicit: file may gain motifs
fimo <- fimo[fimo$p.value < 1e-4, ]
stopifnot(nrow(fimo) > 0L)

motif2_hits <- GRanges(
  seqnames = fimo$sequence_name,
  ranges   = IRanges(start = fimo$start, end = fimo$stop)   # FIMO is 1-based inclusive
)
motif2_hits <- keepSeqlevels(motif2_hits, good_chr, pruning.mode = "coarse")

enh_motif2_pos <- enhancers_filt[countOverlaps(enhancers_filt, motif2_hits) > 0L]

cat(sprintf("joint peaks=%d | enhancers=%d | motif2 hits=%d | motif2+ enhancers=%d (%.1f%%)\n",
            length(joint_peaks_filt), length(enhancers_filt),
            length(motif2_hits), length(enh_motif2_pos),
            100 * length(enh_motif2_pos) / length(enhancers_filt)))

# ---- 4. Permutation test -----------------------------------------------------
# randomizeRegions conserves peak NUMBER and interval WIDTH, and (with
# per.chromosome = TRUE) the per-chromosome peak count, but not chromatin
# state. Positions are drawn from the whole genome with no accessibility mask.
# This caveat is addressed by the control factors described in Chapter 4.
# count.once = TRUE => an enhancer overlapped twice still counts once, so the
# statistic is "number of peaks overlapping >= 1 Motif2+ enhancer".
set.seed(SEED)
pt_motif2_enh <- permTest(
  A                  = joint_peaks_filt,
  B                  = enh_motif2_pos,
  randomize.function = randomizeRegions,
  evaluate.function  = numOverlaps,
  ntimes             = NTIMES,
  genome             = "dm6",
  per.chromosome     = TRUE,
  count.once         = TRUE,
  mc.cores           = NCORES
)

obs <- pt_motif2_enh$numOverlaps$observed
nv  <- pt_motif2_enh$numOverlaps$permuted
z   <- (obs - mean(nv)) / sd(nv)

cat(sprintf("\nobserved=%d  null mean=%.2f  null sd=%.2f  z=%.2f  p=%.2e  (n=%d)\n",
            obs, mean(nv), sd(nv), z, pt_motif2_enh$numOverlaps$pval, length(nv)))

# Regression check against the expected values reported in results/reported_values.tsv.
# The observed overlap is deterministic given the inputs; the null is seeded. If either of
# these fires, an input file has changed.
if (obs != 109L)
  warning(sprintf("observed = %d, expected 109 -- inputs have changed", obs))
if (abs(z - 4.48) > 0.5)
  warning(sprintf("z = %.2f, expected ~4.48 -- inputs have changed", z))

# ---- 5. Save -----------------------------------------------------------------
# Flat list schema, exactly as consumed by the Fig 3.19 plotting code below
# and by the values recorded in results/reported_values.tsv for Chapter 3.
res <- list(
  observed = obs,
  permuted = nv,
  zscore   = round(z, 2),
  pval     = pt_motif2_enh$numOverlaps$pval,
  ntimes   = length(nv),
  n_peaks  = length(joint_peaks_filt)
)
saveRDS(res, file.path(EXPORT_DIR, "perm_motif2_enh_or.rds"))
saveRDS(pt_motif2_enh, file.path(EXPORT_DIR, "perm_motif2_enh_or_FULL.rds"))
cat(sprintf("Saved: %s/perm_motif2_enh_or.rds\n", EXPORT_DIR))
cat("Done.\n")
