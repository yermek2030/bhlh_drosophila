#!/usr/bin/env Rscript
# ====
# perm_enhancers_tss_matched_all.R: extends trh_tgo version to all four heterodimers, tests promoter-proximal depletion of enhancers against a TSS-proximity-matched null.
#

suppressPackageStartupMessages({
  library(GenomicRanges); library(GenomeInfoDb); library(rtracklayer)
  library(regioneR); library(TxDb.Dmelanogaster.UCSC.dm6.ensGene)
})

## Repository root: the environment variable PROJECT_DIR, or the working directory.
PROJECT_DIR <- normalizePath(Sys.getenv("PROJECT_DIR", unset = "."), mustWork = TRUE)
if (!file.exists(file.path(PROJECT_DIR, "README.md")))
  stop("PROJECT_DIR does not point at the repository root: ", PROJECT_DIR)

REF_DIR   <- file.path(PROJECT_DIR, "data/reference_gene_lists")
FINAL_DIR <- file.path(PROJECT_DIR, "data/final_peaks")
EXP3      <- file.path(PROJECT_DIR, "data", "chapter3")
EXP4      <- file.path(PROJECT_DIR, "data", "chapter4")
dir.create(EXP4, showWarnings = FALSE, recursive = TRUE)

args      <- commandArgs(trailingOnly = TRUE)
NTIMES    <- if (length(args) >= 1 && !is.na(as.integer(args[1]))) as.integer(args[1]) else 5000L
POOL_MULT <- if (length(args) >= 2 && !is.na(as.integer(args[2]))) as.integer(args[2]) else 300L
cat(sprintf("ntimes=%d pool_mult=%d\n", NTIMES, POOL_MULT))

good_chr  <- c("chr2L", "chr2R", "chr3L", "chr3R", "chr4", "chrX")
txdb      <- TxDb.Dmelanogaster.UCSC.dm6.ensGene
extraCols <- c(signalValue = "numeric", pValue = "numeric",
               qValue = "numeric", peak = "integer")

## ---- canonical TSS points ---------------------------------------------------
tss_pt <- resize(transcripts(txdb), width = 1L, fix = "start")
tss_pt <- keepSeqlevels(tss_pt, good_chr, pruning.mode = "coarse")

annotate_tss <- function(gr) {
  gr <- keepSeqlevels(gr, good_chr, pruning.mode = "coarse")
  ov <- distanceToNearest(gr, tss_pt, ignore.strand = TRUE)
  gr$distanceToTSS <- NA_integer_
  gr$distanceToTSS[queryHits(ov)] <- mcols(ov)$distance
  stopifnot(!anyNA(gr$distanceToTSS), all(gr$distanceToTSS >= 0))
  gr
}

## ---- peak sets --------------------------------------------------------------
load_peaks <- function(path) {
  gr <- import(path, format = "BED", extraCols = extraCols)
  gr <- gr[seqnames(gr) %in% good_chr]
  gr <- keepSeqlevels(gr, good_chr, pruning.mode = "coarse")
  gr[!duplicated(gr)]
}
Tgo_clean  <- load_peaks(file.path(FINAL_DIR, "TGO_IDR0.05.final.narrowPeak"))
Sim_clean  <- load_peaks(file.path(FINAL_DIR, "SIM_IDR0.05.final.narrowPeak"))
Dys_clean  <- load_peaks(file.path(FINAL_DIR, "DYS_IDR0.05.final.narrowPeak"))
Sima_clean <- load_peaks(file.path(FINAL_DIR, "SIMA_IDR0.05.final.narrowPeak"))

make_joint <- function(tf, tgo) {
  ov <- findOverlaps(tf, tgo, minoverlap = 1, ignore.strand = TRUE)
  gr <- pintersect(tf[queryHits(ov)], tgo[subjectHits(ov)])
  mcols(gr) <- NULL
  gr <- keepSeqlevels(gr, good_chr, pruning.mode = "coarse")
  gr[!duplicated(gr)]
}

# Trh/Tgo comes from the exported BED, exactly as in Ch3 (614 peaks)
trh_df <- read.table(file.path(EXP3, "Trh_Tgo_joint_peaks_native_width_for_MEME_ChIP.bed"),
                     header = FALSE, stringsAsFactors = FALSE,
                     col.names = c("chr", "start", "end", "name", "score", "strand"))
joint <- list(
  Trh  = GRanges(trh_df$chr, IRanges(trh_df$start + 1L, trh_df$end)),
  Sim  = make_joint(Sim_clean,  Tgo_clean),
  Dys  = make_joint(Dys_clean,  Tgo_clean),
  Sima = make_joint(Sima_clean, Tgo_clean)
)
joint <- lapply(joint, annotate_tss)
stopifnot(length(joint$Trh) == 614L, length(joint$Sim)  == 355L,
          length(joint$Dys) == 870L, length(joint$Sima) == 296L)
# canonical distal counts reported in Ch3/Ch4
stopifnot(sum(joint$Trh$distanceToTSS  > 1000) == 146L,
          sum(joint$Sim$distanceToTSS  > 1000) == 112L,
          sum(joint$Dys$distanceToTSS  > 1000) == 300L,
          sum(joint$Sima$distanceToTSS > 1000) ==  78L)

## ---- enhancers + genome -----------------------------------------------------
enh_df <- read.table(file.path(REF_DIR, "merged_embryo_5-20_enhancers_dm6.bed"),
                     header = FALSE, stringsAsFactors = FALSE)
colnames(enh_df)[1:3] <- c("seqnames", "start", "end")
enhancers <- GRanges(enh_df$seqnames, IRanges(enh_df$start + 1L, enh_df$end))
enhancers <- keepSeqlevels(enhancers, good_chr, pruning.mode = "coarse")
cat(sprintf("enhancers: %d\n", length(enhancers)))

genome_gr <- getGenomeAndMask("dm6", mask = NA)$genome
genome_gr <- keepSeqlevels(genome_gr, good_chr, pruning.mode = "coarse")

## ---- one heterodimer --------------------------------------------------------
run_het <- function(gr, label) {
  set.seed(777)                            # per-heterodimer: Trh must reproduce
  d_obs <- gr$distanceToTSS
  stat  <- function(x) sum(overlapsAny(x, enhancers, ignore.strand = TRUE))
  obs   <- c(all = stat(gr), promoter = stat(gr[d_obs <= 1000]),
             distal = stat(gr[d_obs > 1000]))

  CHUNK <- 10L
  pool <- do.call(c, lapply(seq_len(ceiling(POOL_MULT / CHUNK)), function(i)
    randomizeRegions(rep(gr, CHUNK), genome = genome_gr,
                     per.chromosome = TRUE, allow.overlaps = TRUE)))
  mcols(pool) <- NULL
  pool     <- annotate_tss(pool)
  pool_hit <- overlapsAny(pool, enhancers, ignore.strand = TRUE)
  pool_d   <- pool$distanceToTSS
  cat(sprintf("\n=== %s === n=%d pool=%d (%.1f%% enhancer-overlapping)\n",
              label, length(gr), length(pool), 100 * mean(pool_hit)))

  run_arm <- function(dvec) {
    n    <- length(dvec)
    qs   <- quantile(dvec, probs = seq(0, 1, by = 0.1))
    brks <- c(-Inf, unique(qs[-c(1, length(qs))]), Inf)
    dec_obs  <- cut(dvec,   breaks = brks, labels = FALSE, include.lowest = TRUE)
    dec_pool <- cut(pool_d, breaks = brks, labels = FALSE, include.lowest = TRUE)
    tab  <- table(factor(dec_obs, levels = seq_len(max(dec_obs))))
    by_k <- split(seq_along(pool), dec_pool)
    stopifnot(all(vapply(names(tab), function(k)
      length(by_k[[k]]) >= 10L * tab[[k]], logical(1))))
    list(n = n,
         null_gen = vapply(seq_len(NTIMES), function(i)
           sum(pool_hit[sample.int(length(pool), n)]), numeric(1)),
         null_mat = vapply(seq_len(NTIMES), function(i) {
           idx <- unlist(lapply(names(tab), function(k) sample(by_k[[k]], tab[[k]])))
           sum(pool_hit[idx])
         }, numeric(1)))
  }
  arms <- list(all      = run_arm(d_obs),
               promoter = run_arm(d_obs[d_obs <= 1000]),
               distal   = run_arm(d_obs[d_obs  > 1000]))

  do.call(rbind, lapply(names(arms), function(a) {
    A <- arms[[a]]; o <- obs[[a]]
    f <- function(nv, tag) data.frame(
      heterodimer = label, arm = a, null = tag, n_peaks = A$n, observed = o,
      null_mean = mean(nv), null_sd = sd(nv), z = (o - mean(nv)) / sd(nv),
      p_emp = (sum(if (o >= mean(nv)) nv >= o else nv <= o) + 1) / (NTIMES + 1),
      fold = o / mean(nv), n_perm = NTIMES, row.names = NULL)
    rbind(f(A$null_gen, "genomic"), f(A$null_mat, "tss_matched"))
  }))
}

summ <- do.call(rbind, lapply(names(joint), function(h) run_het(joint[[h]], h)))
cat("\n")
print(summ, digits = 4)

write.csv(summ, file.path(EXP4, "perm_enhancers_distmatched_all.csv"), row.names = FALSE)
saveRDS(list(summ = summ, ntimes = NTIMES, pool_mult = POOL_MULT, seed = 777L),
        file.path(EXP4, "perm_enhancers_distmatched_all.rds"))
cat("\nwrote perm_enhancers_distmatched_all.{csv,rds} to", EXP4, "\n")
