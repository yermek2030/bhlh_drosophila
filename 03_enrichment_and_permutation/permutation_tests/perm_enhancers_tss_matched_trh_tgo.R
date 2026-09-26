#!/usr/bin/env Rscript
# ============================================================
# perm_enhancers_tss_matched_trh_tgo.R: enhancer overlap permutation for the three Trh/Tgo strata, null matched on unsigned distance to nearest TSS (TxDb.Dmelanogaster.UCSC.dm6.ensGene); proximal <=1000, distal >1000.
# Statistic: peaks overlapping >=1 enhancer (

suppressPackageStartupMessages({
  library(GenomicRanges); library(GenomeInfoDb)
  library(regioneR); library(TxDb.Dmelanogaster.UCSC.dm6.ensGene)
})

## Repository root: the environment variable PROJECT_DIR, or the working directory.
PROJECT_DIR <- normalizePath(Sys.getenv("PROJECT_DIR", unset = "."), mustWork = TRUE)
if (!file.exists(file.path(PROJECT_DIR, "README.md")))
  stop("PROJECT_DIR does not point at the repository root: ", PROJECT_DIR)

REF_DIR <- file.path(PROJECT_DIR, "data/reference_gene_lists")
EXP3    <- file.path(PROJECT_DIR, "data", "chapter3")
dir.create(EXP3, showWarnings = FALSE, recursive = TRUE)

args      <- commandArgs(trailingOnly = TRUE)
NTIMES    <- if (length(args) >= 1 && !is.na(as.integer(args[1]))) as.integer(args[1]) else 5000L
POOL_MULT <- if (length(args) >= 2 && !is.na(as.integer(args[2]))) as.integer(args[2]) else 300L
set.seed(777)
cat(sprintf("ntimes=%d pool_mult=%d\n", NTIMES, POOL_MULT))

good_chr <- c("chr2L", "chr2R", "chr3L", "chr3R", "chr4", "chrX")
txdb     <- TxDb.Dmelanogaster.UCSC.dm6.ensGene

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

## ---- joint peaks ------------------------------------------------------------
trh_df <- read.table(file.path(EXP3, "Trh_Tgo_joint_peaks_native_width_for_MEME_ChIP.bed"),
                     header = FALSE, stringsAsFactors = FALSE,
                     col.names = c("chr", "start", "end", "name", "score", "strand"))
joint <- GRanges(trh_df$chr, IRanges(trh_df$start + 1L, trh_df$end))
joint <- keepSeqlevels(joint, good_chr, pruning.mode = "coarse")
stopifnot(length(joint) == 614L)
joint <- annotate_tss(joint)
d_obs <- joint$distanceToTSS
stopifnot(sum(d_obs <= 1000) == 468L, sum(d_obs > 1000) == 146L)

## ---- enhancers --------------------------------------------------------------
enh_df <- read.table(file.path(REF_DIR, "merged_embryo_5-20_enhancers_dm6.bed"),
                     header = FALSE, stringsAsFactors = FALSE)
colnames(enh_df)[1:3] <- c("seqnames", "start", "end")
enhancers <- GRanges(enh_df$seqnames, IRanges(enh_df$start + 1L, enh_df$end))
enhancers <- keepSeqlevels(enhancers, good_chr, pruning.mode = "coarse")
cat(sprintf("enhancers: %d\n", length(enhancers)))

stat <- function(gr) sum(overlapsAny(gr, enhancers, ignore.strand = TRUE))
obs <- c(all      = stat(joint),
         promoter = stat(joint[d_obs <= 1000]),
         distal   = stat(joint[d_obs  > 1000]))
stopifnot(obs[["all"]] == 169L, obs[["promoter"]] == 55L, obs[["distal"]] == 114L)

## ---- one candidate pool, shared by both nulls -------------------------------
genome_gr <- getGenomeAndMask("dm6", mask = NA)$genome
genome_gr <- keepSeqlevels(genome_gr, good_chr, pruning.mode = "coarse")

CHUNK <- 10L
pool <- do.call(c, lapply(seq_len(ceiling(POOL_MULT / CHUNK)), function(i)
  randomizeRegions(rep(joint, CHUNK), genome = genome_gr,
                   per.chromosome = TRUE, allow.overlaps = TRUE)))
mcols(pool) <- NULL
pool     <- annotate_tss(pool)
pool_hit <- overlapsAny(pool, enhancers, ignore.strand = TRUE)   # precomputed
pool_d   <- pool$distanceToTSS
cat(sprintf("pool: %d candidates (%.1f%% enhancer-overlapping)\n",
            length(pool), 100 * mean(pool_hit)))

## ---- nulls ------------------------------------------------------------------
# genomic  : sample the pool ignoring TSS distance
# matched  : sample within TSS-distance deciles of the OBSERVED arm, so the
#            permuted set reproduces that arm's promoter-proximity profile
run_arm <- function(dvec) {
  n <- length(dvec)
  qs   <- quantile(dvec, probs = seq(0, 1, by = 0.1))
  brks <- c(-Inf, unique(qs[-c(1, length(qs))]), Inf)
  dec_obs  <- cut(dvec,   breaks = brks, labels = FALSE, include.lowest = TRUE)
  dec_pool <- cut(pool_d, breaks = brks, labels = FALSE, include.lowest = TRUE)
  tab   <- table(factor(dec_obs, levels = seq_len(max(dec_obs))))
  by_k  <- split(seq_along(pool), dec_pool)
  stopifnot(all(vapply(names(tab), function(k)
    length(by_k[[k]]) >= 10L * tab[[k]], logical(1))))

  null_gen <- vapply(seq_len(NTIMES), function(i)
    sum(pool_hit[sample.int(length(pool), n)]), numeric(1))
  null_mat <- vapply(seq_len(NTIMES), function(i) {
    idx <- unlist(lapply(names(tab), function(k) sample(by_k[[k]], tab[[k]])))
    sum(pool_hit[idx])
  }, numeric(1))
  list(n = n, null_gen = null_gen, null_mat = null_mat)
}

arms <- list(all      = run_arm(d_obs),
             promoter = run_arm(d_obs[d_obs <= 1000]),
             distal   = run_arm(d_obs[d_obs  > 1000]))

summ <- do.call(rbind, lapply(names(arms), function(a) {
  A <- arms[[a]]; o <- obs[[a]]
  f <- function(nv, tag) data.frame(
    arm = a, null = tag, n_peaks = A$n, observed = o,
    null_mean = mean(nv), null_sd = sd(nv),
    z = (o - mean(nv)) / sd(nv),
    p_emp = (sum(if (o >= mean(nv)) nv >= o else nv <= o) + 1) / (NTIMES + 1),
    fold = o / mean(nv), n_perm = NTIMES, row.names = NULL)
  rbind(f(A$null_gen, "genomic"), f(A$null_mat, "tss_matched"))
}))
print(summ, digits = 4)

write.csv(summ, file.path(EXP3, "perm_enhancers_distmatched.csv"), row.names = FALSE)
saveRDS(list(summ = summ, arms = arms, obs = obs, ntimes = NTIMES,
             pool_mult = POOL_MULT, pool_n = length(pool), seed = 777L),
        file.path(EXP3, "perm_enhancers_distmatched.rds"))
cat("\nwrote perm_enhancers_distmatched.{csv,rds} to", EXP3, "\n")
