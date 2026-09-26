#!/usr/bin/env Rscript
# ============================================================
# Site-level register test: mean-profile Z-scores not comparable (n=165 in-peak vs n=22925 background).
# Effect size = register-0 score minus site's own flank, measured per site, compared between sets.
# Reuses cached genome-wide FIMO from 02_split_scan_register.R. ~3 min.
# ============================================================
## Repository root: environment variable PROJECT_DIR, or working directory.
PROJECT_DIR <- normalizePath(Sys.getenv("PROJECT_DIR", unset = "."), mustWork = TRUE)
if (!file.exists(file.path(PROJECT_DIR, "README.md")))
  stop("PROJECT_DIR does not point at the repository root: ", PROJECT_DIR)
setwd(PROJECT_DIR)

suppressPackageStartupMessages({
  library(Biostrings); library(GenomicRanges); library(IRanges)
})

OUT  <- file.path(PROJECT_DIR, "data", "motif2_concat_audit")
GENOME  <- file.path(PROJECT_DIR, "reference", "dm6_hardmasked.fa")
PEAKBED <- file.path(PROJECT_DIR, "data", "chapter3",
                     "Trh_Tgo_joint_peaks_filt.bed")
FIMO <- file.path(OUT, "fimo_core9_genome", "fimo.tsv")
stopifnot(file.exists(FIMO), file.exists(GENOME), file.exists(PEAKBED))

BG   <- c(A = .306, C = .194, G = .194, T = .306)
MAIN <- c("chr2L", "chr2R", "chr3L", "chr3R", "chr4", "chrX")
W    <- 30L; LEN <- 2L * W + 6L

x <- readLines(file.path(PROJECT_DIR, "data", "chapter3",
                         "STRCACTGARARAAA.meme"), warn = FALSE)
i <- grep("^letter-probability matrix", x)[1]
p <- do.call(rbind, lapply(strsplit(trimws(x[(i + 1):(i + 15)]),
                                    "[[:space:]]+"), as.numeric))
lo_tail <- log2(((p[10:15, ] + 0.01 * matrix(BG, 6, 4, byrow = TRUE)) / 1.01) /
                  matrix(BG, 6, 4, byrow = TRUE))
max_tail <- sum(apply(lo_tail, 1, max))     # best attainable tail score

gnm <- readDNAStringSet(GENOME); names(gnm) <- sub(" .*", "", names(gnm))
gnm <- gnm[MAIN]; chrlen <- setNames(width(gnm), names(gnm))

fi <- read.table(FIMO, header = TRUE, sep = "\t", comment.char = "#",
                 quote = "", stringsAsFactors = FALSE)
names(fi)[names(fi) == "p.value"] <- "pval"
fi <- fi[fi$sequence_name %in% MAIN & fi$pval < 1e-4, ]
core_gr <- GRanges(fi$sequence_name, IRanges(fi$start, fi$stop), strand = fi$strand)

pk       <- read.table(PEAKBED, header = FALSE, stringsAsFactors = FALSE)[, 1:3]
peaks_gr <- GRanges(pk[[1]], IRanges(pk[[2]] + 1L, pk[[3]]))

ctx_for <- function(gr) {
  st <- as.character(strand(gr))
  s  <- ifelse(st == "+", end(gr) + 1L - W, start(gr) - 6L - W)
  e  <- s + LEN - 1L
  ok <- s >= 1L & e <= chrlen[as.character(seqnames(gr))]
  list(ok = ok, seq = local({
    gr <- gr[ok]; s <- s[ok]; e <- e[ok]; st <- st[ok]
    chr <- as.character(seqnames(gr))
    out <- DNAStringSet(character(0)); ord <- integer(0)
    for (cc in unique(chr)) {
      k <- which(chr == cc)
      v <- as(Views(gnm[[cc]], start = s[k], end = e[k]), "DNAStringSet")
      neg <- st[k] == "-"; if (any(neg)) v[neg] <- reverseComplement(v[neg])
      out <- c(out, v); ord <- c(ord, k)
    }
    out[order(ord)]
  }))
}
score_windows <- function(dss) {
  mat <- as.matrix(dss); idx <- match(mat, c("A", "C", "G", "T")); dim(idx) <- dim(mat)
  s <- numeric(nrow(mat))
  for (j in seq_len(ncol(mat)))
    s <- s + ifelse(is.na(idx[, j]), -10, lo_tail[j, pmax(idx[, j], 1L)])
  s
}

cx  <- ctx_for(core_gr)
gr2 <- core_gr[cx$ok]; ctx <- cx$seq
inpk <- countOverlaps(gr2, peaks_gr, ignore.strand = TRUE) > 0

s0    <- score_windows(subseq(ctx, start = W + 1L, width = 6L))      # register 0
flank_off <- c(-30:-15, 15:30)
fl <- vapply(flank_off, function(d)
  score_windows(subseq(ctx, start = W + d + 1L, width = 6L)), numeric(length(ctx)))
flmean <- rowMeans(fl)
delta  <- s0 - flmean                       # per-site effect size

site <- data.frame(set = ifelse(inpk, "in_joint_peaks", "genome_background"),
                   s0 = s0, flank = flmean, delta = delta)
write.csv(site, file.path(OUT, "T2c_site_scores.csv"), row.names = FALSE)

# ---- 1. effect size per site, comparable across sets -------------------------
agg <- do.call(rbind, lapply(split(site, site$set), function(d) data.frame(
  set = d$set[1], n = nrow(d),
  median_register0 = round(median(d$s0), 2),
  median_flank     = round(median(d$flank), 2),
  mean_delta       = round(mean(d$delta), 2),
  median_delta     = round(median(d$delta), 2))))
cat("\n=== per-site register-0 effect (log2 units above own flank) ===\n")
print(agg, row.names = FALSE)

wt <- wilcox.test(delta ~ set, data = site)
cat(sprintf("\nMann-Whitney on per-site delta: W = %.0f, p = %.3g\n",
            unname(wt$statistic), wt$p.value))
cat(sprintf("Flank composition matched? in-peak median flank %.2f vs background %.2f\n",
            agg$median_flank[agg$set == "in_joint_peaks"],
            agg$median_flank[agg$set == "genome_background"]))

# ---- 2. is it a subset of sites, or a global shift? -------------------------
# "Tail-positive" = register-0 tail within 2 log2 units of the best possible.
THR <- max_tail - 2
site$tail_pos <- site$s0 >= THR
tb <- table(factor(site$set, c("in_joint_peaks", "genome_background")), site$tail_pos)
ft <- fisher.test(tb[, c("TRUE", "FALSE")])
cat(sprintf("\nTail-positive (register-0 score >= %.2f, i.e. within 2 of max %.2f):\n",
            THR, max_tail))
cat(sprintf("  in peaks   : %d/%d (%.1f%%)\n", tb[1, "TRUE"], sum(tb[1, ]),
            100 * tb[1, "TRUE"] / sum(tb[1, ])))
cat(sprintf("  background : %d/%d (%.1f%%)\n", tb[2, "TRUE"], sum(tb[2, ]),
            100 * tb[2, "TRUE"] / sum(tb[2, ])))
cat(sprintf("  OR = %.2f [%.2f-%.2f], p = %.3g\n", unname(ft$estimate),
            ft$conf.int[1], ft$conf.int[2], ft$p.value))

# ---- 3. off-register control: does the same excess appear at offset +-5? ----
ctrl <- vapply(c(-5L, 5L), function(d)
  score_windows(subseq(ctx, start = W + d + 1L, width = 6L)), numeric(length(ctx)))
ctrl_pos <- (ctrl[, 1] >= THR) | (ctrl[, 2] >= THR)
tbc <- table(factor(site$set, c("in_joint_peaks", "genome_background")), ctrl_pos)
ftc <- fisher.test(tbc[, c("TRUE", "FALSE")])
cat(sprintf("\nOff-register control (offset -5 or +5): OR = %.2f, p = %.3g\n",
            unname(ftc$estimate), ftc$p.value))
cat("  A specific register effect requires OR(register 0) >> OR(off-register).\n")

res <- data.frame(
  n_inpeak = sum(inpk), n_background = sum(!inpk),
  delta_inpeak = agg$mean_delta[agg$set == "in_joint_peaks"],
  delta_background = agg$mean_delta[agg$set == "genome_background"],
  delta_ratio = round(agg$mean_delta[agg$set == "in_joint_peaks"] /
                        agg$mean_delta[agg$set == "genome_background"], 2),
  wilcox_p = signif(wt$p.value, 3),
  tailpos_OR = round(unname(ft$estimate), 2), tailpos_p = signif(ft$p.value, 3),
  offreg_OR = round(unname(ftc$estimate), 2), offreg_p = signif(ftc$p.value, 3))
write.csv(res, file.path(OUT, "T2c_register_sitelevel.csv"), row.names = FALSE)

cat("\nVERDICT RULE (site level):\n",
    " tailpos_OR >> offreg_OR and wilcox p < 0.05 -> register is real AND peak-specific.\n",
    " tailpos_OR ~= offreg_OR                     -> in-peak sequence is just A-rich; not a register.\n",
    " tailpos_OR ~= 1                             -> the 15-bp call is a concatenation.\n")
cat("\nWrote T2c_* to", OUT, "\n")
