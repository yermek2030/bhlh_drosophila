# ============================================================
# motif2_robustness_c_halves.R
#   R7 CME string audit; R8 conditional register test; R9 AUROC full vs core vs tail;
#   R10 repeat overlap; R11 CentriMo summit centrality; R12 empirical FIMO FDR
# Run: MEME_BIN=... Rscript motif2_robustness_c_halves.R
# ============================================================

## Repository root: the environment variable PROJECT_DIR, or the working directory.
PROJECT_DIR <- normalizePath(Sys.getenv("PROJECT_DIR", unset = "."), mustWork = TRUE)
if (!file.exists(file.path(PROJECT_DIR, "README.md")))
  stop("PROJECT_DIR does not point at the repository root: ", PROJECT_DIR)
setwd(PROJECT_DIR)

suppressPackageStartupMessages({
  library(GenomicRanges); library(Biostrings); library(rtracklayer)
  library(BSgenome.Dmelanogaster.UCSC.dm6)
})

OUT <- file.path(PROJECT_DIR, "data", "motif2_robustness")
MEME_BIN <- Sys.getenv("MEME_BIN", unset = "")
if (nzchar(MEME_BIN)) Sys.setenv(PATH = paste(MEME_BIN, Sys.getenv("PATH"), sep = ":"))

JOINT_FA  <- file.path(PROJECT_DIR, "data/chapter3/joint.fa")
JOINT_BED <- file.path(PROJECT_DIR, "data/chapter3/Trh_Tgo_joint_peaks_filt.bed")
PWM_FILE  <- file.path(PROJECT_DIR, "data/chapter3/STRCACTGARARAAA.meme")
SHUF_FA   <- file.path(OUT, "joint_dishuf2.fa")
resC <- list()

lp <- readLines(PWM_FILE, warn = FALSE)
i0 <- grep("^letter-probability matrix", lp)
pm <- do.call(rbind, lapply(strsplit(trimws(lp[(i0 + 1):(i0 + 15)]), "\\s+"), as.numeric))
colnames(pm) <- c("A", "C", "G", "T")
best <- paste(colnames(pm)[apply(pm, 1, which.max)], collapse = "")

# ------------------------------------------------------------
# TEST R7, is the CME (ACGTG) inside Motif2 at all?
# ------------------------------------------------------------
d  <- DNAString(best); rcd <- reverseComplement(d)
cme_f <- length(matchPattern("ACGTG", d)); cme_r <- length(matchPattern("ACGTG", rcd))
# highest probability with which the PWM can spell ACGTG at any offset
best_cme_p <- max(vapply(0:10, function(o)
  prod(vapply(1:5, function(k) pm[o + k, substring("ACGTG", k, k)], numeric(1))), numeric(1)))
eb_f <- length(matchPattern("CANNTG", d, fixed = FALSE))
eb_r <- length(matchPattern("CANNTG", rcd, fixed = FALSE))
resC$R7 <- data.frame(modal_consensus = best, CME_ACGTG_fwd = cme_f, CME_ACGTG_rev = cme_r,
                      max_PWM_prob_of_ACGTG = signif(best_cme_p, 3),
                      CANNTG_fwd = eb_f, CANNTG_rev = eb_r)
cat("\n== TEST R7  CME / E-box string audit ==\n"); print(resC$R7, row.names = FALSE)

# ------------------------------------------------------------
# TEST R8, conditional register test (the real decomposition test)
# ------------------------------------------------------------
fimo_txt <- function(meme, fa_path, tag, thresh = "1.0E-3") {
  f <- file.path(OUT, sprintf("fimoC_%s.tsv", tag))
  system(sprintf("fimo --verbosity 1 --bgfile --nrdb-- --thresh %s --text %s %s > %s 2>/dev/null",
                 thresh, shQuote(meme), shQuote(fa_path), shQuote(f)))
  x <- tryCatch(read.delim(f, stringsAsFactors = FALSE), error = function(e) NULL)
  if (is.null(x) || !nrow(x)) return(NULL)
  x <- x[!startsWith(as.character(x[[1]]), "#") & nzchar(as.character(x[[1]])), , drop = FALSE]
  x$start <- as.integer(x$start); x$stop <- as.integer(x$stop); x
}
LM <- file.path(OUT, "m2_left.meme"); RM <- file.path(OUT, "m2_right.meme")

register_stats <- function(fa_path, tag) {
  L <- fimo_txt(LM, fa_path, paste0("L_", tag)); R <- fimo_txt(RM, fa_path, paste0("R_", tag))
  if (is.null(L)) return(list(nL = 0L, n0 = 0L, frac = NA_real_))
  key <- paste(R$sequence_name, R$strand)
  Rby <- split(R, key)
  n0 <- 0L
  for (i in seq_len(nrow(L))) {
    rr <- Rby[[paste(L$sequence_name[i], L$strand[i])]]
    if (is.null(rr)) next
    d <- if (L$strand[i] == "+") rr$start - L$stop[i] - 1L else L$start[i] - rr$stop[i] - 1L
    if (any(d == 0L, na.rm = TRUE)) n0 <- n0 + 1L
  }
  list(nL = nrow(L), n0 = n0, frac = n0 / nrow(L))
}
ro <- register_stats(JOINT_FA, "obs")
rs <- register_stats(SHUF_FA, "shuf")
ft8 <- fisher.test(matrix(c(ro$n0, ro$nL - ro$n0, rs$n0, rs$nL - rs$n0), nrow = 2))
resC$R8 <- data.frame(
  real_core_hits = ro$nL, real_with_tail_at_register0 = ro$n0,
  real_pct = round(100 * ro$frac, 2),
  shuf_core_hits = rs$nL, shuf_with_tail_at_register0 = rs$n0,
  shuf_pct = round(100 * rs$frac, 2),
  OR = round(unname(ft8$estimate), 2), p = signif(ft8$p.value, 3))
cat("\n== TEST R8  conditional register (decomposition test) ==\n")
print(resC$R8, row.names = FALSE)

# ------------------------------------------------------------
# TEST R9, threshold-free discrimination (max log-odds per sequence, AUROC)
# ------------------------------------------------------------
bg <- c(A = 0.306, C = 0.194, G = 0.194, T = 0.306)
mkpwm <- function(m) log2((m + 1e-3) / matrix(bg[colnames(m)], nrow(m), 4, byrow = TRUE))
scan_max <- function(seqs, pw) {
  hits <- Biostrings::PWMscoreStartingAt
  vapply(seq_along(seqs), function(i) {
    s <- seqs[[i]]
    f <- max(hits(pw, s, starting.at = seq_len(max(1, length(s) - ncol(pw) + 1))))
    r <- max(hits(pw, reverseComplement(s), starting.at = seq_len(max(1, length(s) - ncol(pw) + 1))))
    max(f, r)
  }, numeric(1))
}
auroc <- function(pos, neg) {
  r <- rank(c(pos, neg)); (sum(r[seq_along(pos)]) - length(pos) * (length(pos) + 1) / 2) /
    (length(pos) * length(neg))
}
real <- readDNAStringSet(JOINT_FA)
shuf <- readDNAStringSet(SHUF_FA)
set.seed(11); shuf <- shuf[sample(length(shuf), 6140)]
r9 <- do.call(rbind, lapply(list(full = 1:15, core_ebox = 1:9, tail_arich = 10:15),
  function(ix) NULL))
r9 <- do.call(rbind, lapply(names(list(full = 1:15, core_ebox = 1:9, tail_arich = 10:15)),
  function(nm) {
    ix <- list(full = 1:15, core_ebox = 1:9, tail_arich = 10:15)[[nm]]
    pw <- t(mkpwm(pm[ix, , drop = FALSE]))
    a <- auroc(scan_max(real, pw), scan_max(shuf, pw))
    data.frame(submotif = nm, width = length(ix), AUROC = round(a, 4))
  }))
resC$R9 <- r9
cat("\n== TEST R9  threshold-free AUROC (real vs dinucleotide-shuffled) ==\n")
print(r9, row.names = FALSE)

# ------------------------------------------------------------
# TEST R10, are Motif2 instances simple repeats?
# ------------------------------------------------------------
rep_bed <- file.path(PROJECT_DIR, "reference", "dm6_repeats.bed")
peaks <- import(JOINT_BED, format = "BED")
ob <- fimo_txt(PWM_FILE, JOINT_FA, "m2_1e4", thresh = "1.0E-4")
pk <- peaks[match(ob$sequence_name, peaks$name)]
inst <- GRanges(seqnames(pk), IRanges(start(pk) + ob$start - 1L, start(pk) + ob$stop - 1L))
if (file.exists(rep_bed)) {
  rp <- import(rep_bed, format = "BED")
  frac_rep <- mean(overlapsAny(inst, rp))
  ctrl <- GRanges(seqnames(peaks), IRanges(start(peaks) + floor((width(peaks) - 15) / 2),
                                           width = 15L))
  frac_ctrl <- mean(overlapsAny(ctrl, rp))
  resC$R10 <- data.frame(n_instances = length(inst),
                         pct_in_repeats = round(100 * frac_rep, 1),
                         pct_peak_centres_in_repeats = round(100 * frac_ctrl, 1))
} else resC$R10 <- data.frame(note = "dm6_repeats.bed missing")
cat("\n== TEST R10  repeat overlap ==\n"); print(resC$R10, row.names = FALSE)

# ------------------------------------------------------------
# TEST R11, summit centrality (CentriMo at Trh summits), Motif2 vs CME
# ------------------------------------------------------------
trh <- read.delim(file.path(PROJECT_DIR, "data/final_peaks/TRH_IDR0.05.clean.narrowPeak"),
                  header = FALSE, stringsAsFactors = FALSE)
trh_gr <- GRanges(trh$V1, IRanges(trh$V2 + 1L, trh$V3), summit = trh$V10)
h <- findOverlaps(peaks, trh_gr)
sm <- tapply(trh_gr$summit[subjectHits(h)], queryHits(h), function(x) x[1])
qi <- as.integer(names(sm))
cen <- GRanges(seqnames(peaks)[qi], IRanges(as.integer(sm) - 249L, as.integer(sm) + 250L))
cen <- trim(cen); cen <- cen[width(cen) == 500L]
sq <- getSeq(BSgenome.Dmelanogaster.UCSC.dm6, cen)
names(sq) <- paste0("summit_", seq_along(sq))
writeXStringSet(sq, file.path(OUT, "trh_summits_500bp.fa"))
cme_meme <- file.path(OUT, "cme.meme")
writeLines(c("MEME version 5", "", "ALPHABET= ACGT", "", "strands: + -", "",
             "Background letter frequencies (from unknown source):",
             "A 0.306 C 0.194 G 0.194 T 0.306", "", "MOTIF ACGTG CME",
             "letter-probability matrix: alength= 4 w= 5 nsites= 100 E= 0",
             "1 0 0 0", "0 1 0 0", "0 0 1 0", "0 0 0 1", "0 0 1 0"), cme_meme)
system(sprintf("centrimo --verbosity 1 --oc %s --local %s %s %s >/dev/null 2>&1",
               shQuote(file.path(OUT, "centrimo_summits")),
               shQuote(file.path(OUT, "trh_summits_500bp.fa")), shQuote(PWM_FILE),
               shQuote(cme_meme)))
ctsv <- file.path(OUT, "centrimo_summits", "centrimo.tsv")
if (file.exists(ctsv)) {
  ct <- read.delim(ctsv, stringsAsFactors = FALSE, comment.char = "#")
  ct <- ct[nzchar(as.character(ct[[1]])), , drop = FALSE]
  resC$R11 <- ct[, intersect(c("motif_id", "consensus", "E.value", "adj_p.value",
                               "bin_location", "bin_width", "sites_in_bin", "total_sites"),
                             names(ct))]
} else resC$R11 <- data.frame(note = "centrimo produced no tsv")
cat("\n== TEST R11  summit centrality (n=", length(sq), " summits) ==\n", sep = "")
print(resC$R11, row.names = FALSE)

# ------------------------------------------------------------
# TEST R12, empirical FIMO FDR from the shuffled scan
# ------------------------------------------------------------
obs_hits <- nrow(ob)
sh <- fimo_txt(PWM_FILE, SHUF_FA, "m2_shuf_1e4", thresh = "1.0E-4")
exp_hits <- nrow(sh) / 50
resC$R12 <- data.frame(observed_hits = obs_hits,
                       expected_hits_dinuc_null = round(exp_hits, 1),
                       empirical_FDR_pct = round(100 * exp_hits / obs_hits, 1),
                       fold = round(obs_hits / exp_hits, 2))
cat("\n== TEST R12  empirical FDR ==\n"); print(resC$R12, row.names = FALSE)

saveRDS(resC, file.path(OUT, "motif2_robustness_partC.rds"))
cat("\nPart C written to ", OUT, "\n", sep = "")