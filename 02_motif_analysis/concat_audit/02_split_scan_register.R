#!/usr/bin/env Rscript
# ====
# Split Motif2 PWM at CACTG|ARARAAA boundary; check genome-wide frequency of each half, tail register given a core hit, and whether register is peak-specific or genome-wide.
# Test (c): composite element gives sharp offset-0 spike in peaks and flat profile outside; AT-rich concatenation artefact gives the same profile in both.
# Depends: fimo via conda run -n phd; Biostrings, GenomicRanges. Runtime: FIMO ~5-15 min (cached on rerun), R part ~2 min.
# ====
## Repository root: PROJECT_DIR environment variable, or working directory.
PROJECT_DIR <- normalizePath(Sys.getenv("PROJECT_DIR", unset = "."), mustWork = TRUE)
if (!file.exists(file.path(PROJECT_DIR, "README.md")))
  stop("PROJECT_DIR does not point at the repository root: ", PROJECT_DIR)
setwd(PROJECT_DIR)

suppressPackageStartupMessages({
  library(Biostrings); library(GenomicRanges); library(IRanges)
})

OUT <- file.path(PROJECT_DIR, "data", "motif2_concat_audit")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

PWM_FILE <- file.path(PROJECT_DIR, "data", "chapter3",
                      "STRCACTGARARAAA.meme")
GENOME   <- file.path(PROJECT_DIR, "reference", "dm6_hardmasked.fa")
PEAKBED  <- file.path(PROJECT_DIR, "data", "chapter3",
                      "Trh_Tgo_joint_peaks_filt.bed")
stopifnot(file.exists(PWM_FILE), file.exists(GENOME), file.exists(PEAKBED))

BG   <- c(A = .306, C = .194, G = .194, T = .306)
CORE <- 1:9    # STRCACTGA
TAIL <- 10:15  # RARAAA
MAIN <- c("chr2L", "chr2R", "chr3L", "chr3R", "chr4", "chrX")

# ---------------------------------------------------------------- 1. split ---
read_meme_lpm <- function(f) {
  x <- readLines(f); i <- grep("^letter-probability matrix", x)[1]
  w <- as.integer(sub(".*w= *([0-9]+).*", "\\1", x[i]))
  m <- do.call(rbind, lapply(strsplit(trimws(x[(i + 1):(i + w)]),
                                      "[[:space:]]+"), as.numeric))
  colnames(m) <- c("A", "C", "G", "T"); m
}
write_meme <- function(m, id, f) {
  writeLines(c("MEME version 5", "", "ALPHABET= ACGT", "", "strands: + -", "",
    "Background letter frequencies (from unknown source):",
    sprintf("A %.3f C %.3f G %.3f T %.3f", BG[1], BG[2], BG[3], BG[4]),
    "", paste("MOTIF", id),
    sprintf("letter-probability matrix: alength= 4 w= %d nsites= 36 E= 0", nrow(m)),
    apply(m, 1, function(r) paste(sprintf("%.6f", r), collapse = " "))), f)
}
p <- read_meme_lpm(PWM_FILE); stopifnot(nrow(p) == 15L)

f_core <- file.path(OUT, "m2_core9.meme")
f_tail <- file.path(OUT, "m2_tail6.meme")
write_meme(p[CORE, ], "M2_CORE9_STRCACTGA", f_core)
write_meme(p[TAIL, ], "M2_TAIL6_RARAAA",    f_tail)

# --- minimum attainable FIMO p-value: a SHORT PWM CANNOT reach p < 1e-4 ------
min_attainable_p <- function(m, bg, pseudo = 0.01, scale = 100) {
  q  <- (m + pseudo * matrix(bg, nrow(m), 4, byrow = TRUE)) / (1 + pseudo)
  S  <- log2(sweep(q, 2, bg, "/"))
  Si <- round((S - apply(S, 1, min)) * scale)          # non-negative integers
  span <- sum(apply(Si, 1, max))
  pmf <- numeric(span + 1); pmf[1] <- 1                # index = score + 1
  for (j in seq_len(nrow(m))) {
    nxt <- numeric(span + 1)
    nz  <- which(pmf > 0)
    for (b in 1:4) {
      k <- nz + Si[j, b]
      nxt[k] <- nxt[k] + pmf[nz] * bg[b]
    }
    pmf <- nxt
  }
  sum(pmf[max(which(pmf > 0))])
}
mp <- data.frame(
  pwm  = c("full_15", "core_9", "tail_6"),
  minp = c(min_attainable_p(p, BG), min_attainable_p(p[CORE, ], BG),
           min_attainable_p(p[TAIL, ], BG)))
mp$usable_at_1e4 <- ifelse(mp$minp < 1e-4, "YES", "NO - returns 0 hits by construction")
cat("\n=== minimum attainable FIMO p-value per PWM ===\n"); print(mp, row.names = FALSE)
write.csv(mp, file.path(OUT, "T2_min_attainable_p.csv"), row.names = FALSE)

# ------------------------------------------------- 2. genome-wide core scan ---
FIMO_DIR <- file.path(OUT, "fimo_core9_genome")
if (!file.exists(file.path(FIMO_DIR, "fimo.tsv"))) {
  cat("\nRunning genome-wide FIMO (core 9-mer, p < 1e-4) -- several minutes...\n")
  system2("conda", c("run", "-n", "phd", "fimo",
                     "--oc", shQuote(FIMO_DIR), "--verbosity", "1",
                     "--bgfile", "--nrdb--", "--thresh", "1.0E-4",
                     "--max-stored-scores", "10000000",
                     shQuote(f_core), shQuote(GENOME)))
}
stopifnot(file.exists(file.path(FIMO_DIR, "fimo.tsv")))

fi <- read.table(file.path(FIMO_DIR, "fimo.tsv"), header = TRUE, sep = "\t",
                 comment.char = "#", quote = "", stringsAsFactors = FALSE)
names(fi)[names(fi) == "p.value"] <- "pval"
fi <- fi[fi$sequence_name %in% MAIN & fi$pval < 1e-4, ]

gnm <- readDNAStringSet(GENOME); names(gnm) <- sub(" .*", "", names(gnm))
gnm <- gnm[MAIN]
chrlen <- setNames(width(gnm), names(gnm))
eff_mb <- sum(as.numeric(letterFrequency(gnm, "ACGT"))) / 1e6

cat(sprintf("\nEffective (non-N) main-chromosome length: %.1f Mb\n", eff_mb))
cat(sprintf("Core 9-mer genome-wide hits (p<1e-4): %d  (%.1f per Mb)\n",
            nrow(fi), nrow(fi) / eff_mb))
cat("The 6-bp tail is not scanned genome-wide: min attainable p =",
    signif(mp$minp[3], 3), "- p<1e-4 returns 0 hits by construction.\n",
    "Its occupancy is measured by direct log-odds scoring below instead.\n")

# --------------------------------------- 3. register profile around the core ---
core_gr <- GRanges(fi$sequence_name, IRanges(fi$start, fi$stop), strand = fi$strand)
pk       <- read.table(PEAKBED, header = FALSE, stringsAsFactors = FALSE)[, 1:3]
peaks_gr <- GRanges(pk[[1]], IRanges(pk[[2]] + 1L, pk[[3]]))
in_peak  <- countOverlaps(core_gr, peaks_gr, ignore.strand = TRUE) > 0
cat(sprintf("Core hits inside Trh/Tgo joint peaks: %d ; outside: %d\n",
            sum(in_peak), sum(!in_peak)))

lo_tail <- log2(((p[TAIL, ] + 0.01 * matrix(BG, 6, 4, byrow = TRUE)) / 1.01) /
                  matrix(BG, 6, 4, byrow = TRUE))

W  <- 30L
LEN <- 2L * W + 6L

# memory-safe context extraction: per-chromosome Views, never name-indexing gnm
ctx_for <- function(gr) {
  st <- as.character(strand(gr))
  s  <- ifelse(st == "+", end(gr) + 1L - W, start(gr) - 6L - W)
  e  <- s + LEN - 1L
  ok <- s >= 1L & e <= chrlen[as.character(seqnames(gr))]
  gr <- gr[ok]; s <- s[ok]; e <- e[ok]; st <- st[ok]
  chr <- as.character(seqnames(gr))
  out <- DNAStringSet(character(0)); ord <- integer(0)
  for (cc in unique(chr)) {
    k <- which(chr == cc)
    v <- as(Views(gnm[[cc]], start = s[k], end = e[k]), "DNAStringSet")
    neg <- st[k] == "-"
    if (any(neg)) v[neg] <- reverseComplement(v[neg])
    out <- c(out, v); ord <- c(ord, k)
  }
  out[order(ord)]
}

score_windows <- function(dss) {
  mat <- as.matrix(dss)
  idx <- match(mat, c("A", "C", "G", "T")); dim(idx) <- dim(mat)
  s <- numeric(nrow(mat))
  for (j in seq_len(ncol(mat)))
    s <- s + ifelse(is.na(idx[, j]), -10, lo_tail[j, pmax(idx[, j], 1L)])
  s
}

offsets <- -W:W
profile_of <- function(gr, label) {
  if (length(gr) == 0L) return(NULL)
  ctx <- ctx_for(gr)
  m <- vapply(offsets, function(d)
    mean(score_windows(subseq(ctx, start = W + d + 1L, width = 6L))), numeric(1))
  data.frame(offset = offsets, mean_tail_score = m, set = label, n = length(ctx))
}
prof <- rbind(profile_of(core_gr[in_peak],  "in_joint_peaks"),
              profile_of(core_gr[!in_peak], "genome_background"))
write.csv(prof, file.path(OUT, "T2_register_profile.csv"), row.names = FALSE)

sharp <- do.call(rbind, lapply(split(prof, prof$set), function(d) {
  d0  <- d$mean_tail_score[d$offset == 0]
  fl  <- mean(d$mean_tail_score[abs(d$offset) >= 15])
  sdf <- sd(d$mean_tail_score[abs(d$offset) >= 15])
  data.frame(set = d$set[1], n = d$n[1], score_at_register0 = round(d0, 3),
             flank_mean = round(fl, 3), flank_sd = round(sdf, 3),
             Z_at_register0 = round((d0 - fl) / sdf, 1))
}))
cat("\n=== register-0 sharpness (tail log-odds, motif orientation) ===\n")
print(sharp, row.names = FALSE)
write.csv(sharp, file.path(OUT, "T2_register_sharpness.csv"), row.names = FALSE)

png(file.path(OUT, "T2_register_profile.png"), width = 1700, height = 950, res = 180)
op <- par(mar = c(4.2, 4.4, 2.6, 1))
plot(NA, xlim = c(-W, W), ylim = range(prof$mean_tail_score),
     xlab = "Offset of 6-bp tail from 3' end of core (bp)",
     ylab = "Mean tail log-odds score",
     main = "Where does the A-tail sit relative to the CACTG core?")
cols <- c(in_joint_peaks = "black", genome_background = "grey60")
for (s in names(cols)) {
  d <- prof[prof$set == s, ]; lines(d$offset, d$mean_tail_score, col = cols[[s]], lwd = 2)
}
abline(v = 0, lty = 3)
legend("topright", names(cols), col = unname(cols), lwd = 2, bty = "n")
par(op); dev.off()

cat("\nVERDICT RULE:\n",
    " sharp spike at 0 in peaks AND flat genome-wide -> real composite element.\n",
    " same spike in both sets                        -> genomic composition artefact.\n",
    " no spike anywhere                              -> the 15-bp call is a concatenation.\n")
cat("\nWrote T2_* to", OUT, "\n")
