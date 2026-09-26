#!/usr/bin/env Rscript
# ============================================================
# perm_overlap_trh_tgo_tss_matched.R
# TSS-distance-decile-matched null for Trh x Tgo peak overlap (IDR peaks, 2000 permutations); writes results/reported_

## Repository root: the environment variable PROJECT_DIR, or the working directory.
PROJECT_DIR <- normalizePath(Sys.getenv("PROJECT_DIR", unset = "."), mustWork = TRUE)
if (!file.exists(file.path(PROJECT_DIR, "README.md")))
  stop("PROJECT_DIR does not point at the repository root: ", PROJECT_DIR)

## Values quoted in the dissertation are kept in results/reported_values.tsv
## (columns section, name, value). Writing a name that is already there
## replaces its value, so scripts can be re-run singly and in any order.
VALUES_FILE <- file.path(PROJECT_DIR, "results", "reported_values.tsv")
save_values <- function(section, values) {
  dir.create(dirname(VALUES_FILE), recursive = TRUE, showWarnings = FALSE)
  tab <- if (file.exists(VALUES_FILE)) {
    read.delim(VALUES_FILE, colClasses = "character", quote = "", comment.char = "")
  } else {
    data.frame(section = character(), name = character(), value = character())
  }
  for (nm in names(values)) {
    v   <- gsub("[\t\r\n]+", " ", paste(as.character(values[[nm]]), collapse = " "))
    tab <- tab[tab$name != nm, , drop = FALSE]
    tab <- rbind(tab, data.frame(section = section, name = nm, value = v))
  }
  tab <- tab[order(tab$section, tab$name), , drop = FALSE]
  write.table(tab, VALUES_FILE, sep = "\t", quote = FALSE, row.names = FALSE)
  cat(sprintf("%d values written to %s\n", length(values), VALUES_FILE))
}

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(GenomeInfoDb)
  library(regioneR)
  library(TxDb.Dmelanogaster.UCSC.dm6.ensGene)
})

args      <- commandArgs(trailingOnly = TRUE)
NPERM     <- if (length(args) >= 1 && !is.na(as.integer(args[1]))) as.integer(args[1]) else 2000L
POOL_MULT <- if (length(args) >= 2 && !is.na(as.integer(args[2]))) as.integer(args[2]) else 300L
set.seed(777)

good_chr <- c("chr2L", "chr2R", "chr3L", "chr3R", "chr4", "chrX")
PEAK_DIR <- file.path(PROJECT_DIR, "data/final_peaks")
OUT_DIR  <- file.path(PROJECT_DIR, "data", "chapter3")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ---- 1. peak sets -----------------------------------------------------------
read_np <- function(f) {
  stopifnot(file.exists(f))
  df <- read.table(f, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
  gr <- GRanges(df$V1, IRanges(df$V2 + 1L, df$V3))          # BED -> 1-based
  sort(keepSeqlevels(gr, good_chr, pruning.mode = "coarse"))
}
trh <- read_np(file.path(PEAK_DIR, "TRH_IDR0.05.clean.narrowPeak"))
tgo <- read_np(file.path(PEAK_DIR, "TGO_IDR0.05.clean.narrowPeak"))
stopifnot(length(trh) == 1033L, length(tgo) == 1403L)

stat <- function(gr) sum(overlapsAny(gr, tgo))
obs  <- stat(trh)
stopifnot(obs == 614L)                                       # == Ev_obs, Fig 3.5B

# ---- 2. TSS and decile definition -------------------------------------------
txdb <- TxDb.Dmelanogaster.UCSC.dm6.ensGene
tss  <- suppressMessages(promoters(genes(txdb), upstream = 0, downstream = 1))
tss  <- trim(keepSeqlevels(tss, good_chr, pruning.mode = "coarse"))

# unsigned distance from peak centre to nearest TSS: signed distanceToTSS
# must never be used for a distance-band definition
dist_tss <- function(gr) {
  m <- distanceToNearest(resize(gr, width = 1L, fix = "center"), tss,
                         ignore.strand = TRUE)
  d <- rep(NA_integer_, length(gr))
  d[queryHits(m)] <- mcols(m)$distance
  d
}

d_obs <- dist_tss(trh)
stopifnot(!anyNA(d_obs))
# decile cut points from the OBSERVED Trh distance distribution
qs   <- quantile(d_obs, probs = seq(0, 1, by = 0.1))
brks <- c(-Inf, unique(qs[-c(1, length(qs))]), Inf)
to_dec  <- function(d) cut(d, breaks = brks, labels = FALSE, include.lowest = TRUE)
dec_obs <- to_dec(d_obs)
n_dec   <- max(dec_obs)
tab_obs <- table(factor(dec_obs, levels = seq_len(n_dec)))
stopifnot(sum(tab_obs) == length(trh))

# ---- 3. one candidate pool, shared by both nulls -----------------------------
sl  <- seqlengths(txdb)[good_chr]
gen <- GRanges(names(sl), IRanges(1L, sl))

CHUNK <- 10L
pool <- do.call(c, lapply(seq_len(ceiling(POOL_MULT / CHUNK)), function(i) {
  randomizeRegions(rep(trh, CHUNK), genome = gen,
                   per.chromosome = TRUE, allow.overlaps = TRUE)
}))
pool_dec <- to_dec(dist_tss(pool))
keep     <- !is.na(pool_dec)
pool     <- pool[keep]; pool_dec <- pool_dec[keep]
pool_by  <- split(seq_along(pool), pool_dec)
# every observed decile must be over-represented in the pool for sampling
stopifnot(all(vapply(seq_len(n_dec),
                     function(k) length(pool_by[[as.character(k)]]) >= 10L * tab_obs[[k]],
                     logical(1))))
cat(sprintf("pool = %d candidate intervals across %d deciles\n", length(pool), n_dec))

# ---- 4. the two nulls --------------------------------------------------------
# genomic: sample the pool ignoring decile strata
null_gen <- vapply(seq_len(NPERM), function(i)
  stat(pool[sample.int(length(pool), length(trh))]), numeric(1))

# matched: sample within each decile, preserving the observed decile counts
null_mat <- vapply(seq_len(NPERM), function(i) {
  idx <- unlist(lapply(seq_len(n_dec), function(k)
    sample(pool_by[[as.character(k)]], tab_obs[[k]])))
  stat(pool[idx])
}, numeric(1))

emp_p <- function(nullv) (sum(nullv >= obs) + 1) / (length(nullv) + 1)
zsc   <- function(nullv) (obs - mean(nullv)) / sd(nullv)

res <- data.frame(
  null        = c("genomic", "tss_decile_matched"),
  observed    = obs,
  null_mean   = c(mean(null_gen), mean(null_mat)),
  null_sd     = c(sd(null_gen),   sd(null_mat)),
  null_q95    = c(quantile(null_gen, 0.95), quantile(null_mat, 0.95)),
  z           = c(zsc(null_gen),  zsc(null_mat)),
  p_emp       = c(emp_p(null_gen), emp_p(null_mat)),
  fold        = obs / c(mean(null_gen), mean(null_mat)),
  n_perm      = NPERM,
  row.names   = NULL
)
inflation  <- mean(null_mat) / mean(null_gen)   # the "2.4-fold" quantity
prom_share <- 100 * (1 - 1 / inflation)         # share of genomic fold from promoter bias

print(res, digits = 4)
cat(sprintf("\nexpected-overlap inflation (matched / genomic) = %.2f-fold\n", inflation))
cat(sprintf("promoter-preference share of the genomic fold  = %.1f%%\n", prom_share))

write.csv(res, file.path(OUT_DIR, "perm_trhtgo_tssmatched.csv"), row.names = FALSE)
saveRDS(list(res = res, inflation = inflation, prom_share = prom_share,
             null_gen = null_gen, null_mat = null_mat, brks = brks,
             tab_obs = tab_obs, seed = 777L, nperm = NPERM,
             pool_mult = POOL_MULT, pool_n = length(pool)),
        file.path(OUT_DIR, "perm_trhtgo_tssmatched.rds"))

# ---- 5. reported values -------------------------------------------------------

fmt_p <- function(p, n) if (p <= 1 / (n + 1) + 1e-12)
  sprintf("< %.2fe-4", 1e4 / (n + 1)) else signif(p, 3)

V <- list(
  permMatchedN         = format(NPERM, scientific = FALSE),
  permMatchedPoolN     = format(length(pool), scientific = FALSE),
  permGenomicNullMean  = sprintf("%.2f", res$null_mean[1]),
  permGenomicFold      = sprintf("%.1f", res$fold[1]),
  permMatchedNullMean  = sprintf("%.2f", res$null_mean[2]),
  permMatchedNullSd    = sprintf("%.2f", res$null_sd[2]),
  permMatchedZ         = sprintf("%.2f", res$z[2]),
  permMatchedFold      = sprintf("%.1f", res$fold[2]),
  permMatchedExpFold   = sprintf("%.1f", inflation),
  permPromoterSharePct = sprintf("%.0f", prom_share)
)
save_values("ch3_permutations", V)
print(unlist(V))
