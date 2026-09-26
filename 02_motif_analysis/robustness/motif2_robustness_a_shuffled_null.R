# ============================================================
# motif2_robustness_a_shuffled_null.R
# Tests Motif2 (STRCACTGARARAAA): not an E-box plus A/T tract, not noise.
# R1 E-box audit, R2 internal spacing, R3 half vs whole, R4 split-half validation, R5 HOT loci check.
# Run: Rscript motif2_robustness_a_shuffled_null.R
# ============================================================

## Repository root: the environment variable PROJECT_DIR, or the working directory.
PROJECT_DIR <- normalizePath(Sys.getenv("PROJECT_DIR", unset = "."), mustWork = TRUE)
if (!file.exists(file.path(PROJECT_DIR, "README.md")))
  stop("PROJECT_DIR does not point at the repository root: ", PROJECT_DIR)
setwd(PROJECT_DIR)

suppressPackageStartupMessages({
  library(GenomicRanges); library(Biostrings); library(rtracklayer)
})

OUT <- file.path(PROJECT_DIR, "data", "motif2_robustness")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

MEME_BIN <- Sys.getenv("MEME_BIN", unset = "")
if (nzchar(MEME_BIN)) Sys.setenv(PATH = paste(MEME_BIN, Sys.getenv("PATH"), sep = ":"))

JOINT_FA  <- file.path(PROJECT_DIR, "data/chapter3/joint.fa")
JOINT_BED <- file.path(PROJECT_DIR, "data/chapter3/Trh_Tgo_joint_peaks_filt.bed")
PWM_FILE  <- file.path(PROJECT_DIR, "data/chapter3/STRCACTGARARAAA.meme")
stopifnot(file.exists(JOINT_FA), file.exists(JOINT_BED), file.exists(PWM_FILE))

res <- list()

# ------------------------------------------------------------
# E-box string audit
# ------------------------------------------------------------
lp <- readLines(PWM_FILE)
i0 <- grep("^letter-probability matrix", lp)
w  <- as.integer(sub(".*w= *([0-9]+).*", "\\1", lp[i0]))
pm <- do.call(rbind, lapply(strsplit(trimws(lp[(i0 + 1):(i0 + w)]), "\\s+"), as.numeric))
colnames(pm) <- c("A", "C", "G", "T")
stopifnot(nrow(pm) == 15L, all(abs(rowSums(pm) - 1) < 1e-6))

iupac <- function(p, thr = 0.25) {
  b <- colnames(pm)[p >= thr]
  if (length(b) == 1L) return(b)
  key <- c(AC = "M", AG = "R", AT = "W", CG = "S", CT = "Y", GT = "K",
           ACG = "V", ACT = "H", AGT = "D", CGT = "B", ACGT = "N")
  unname(key[paste(b, collapse = "")])
}
cons <- paste(apply(pm, 1, iupac), collapse = "")
best <- paste(colnames(pm)[apply(pm, 1, which.max)], collapse = "")

ebox_hits <- function(s) {
  d <- DNAString(s)
  fw <- matchPattern("CANNTG", d, fixed = FALSE)
  rv <- matchPattern("CANNTG", reverseComplement(d), fixed = FALSE)
  list(fwd = length(fw), rev = length(rv))
}
eb <- ebox_hits(best)
res$R1 <- data.frame(consensus_iupac = cons, consensus_best = best,
                     CANNTG_fwd = eb$fwd, CANNTG_rev = eb$rev)
cat("\n== TEST R1  E-box string audit ==\n")
cat(sprintf("IUPAC consensus : %s\nModal consensus : %s\n", cons, best))
cat(sprintf("CANNTG matches  : fwd=%d rev=%d  ->  %s\n", eb$fwd, eb$rev,
            ifelse(eb$fwd + eb$rev > 0, "canonical E-box PRESENT",
                   "NO canonical CANNTG E-box in the consensus")))

# ------------------------------------------------------------
# Test R2: is Motif2 one unit or two loose features. L=pos1-9 (STRCACTGA) E-box-like, R=pos10-15 (RARAAA) A-rich. Fixed offset means unit, broad spread means concatenation.
# ------------------------------------------------------------
write_meme <- function(mat, name, path, nsites = 36) {
  con <- file(path, "w")
  writeLines(c("MEME version 5", "", "ALPHABET= ACGT", "", "strands: + -", "",
               "Background letter frequencies (from unknown source):",
               "A 0.306 C 0.194 G 0.194 T 0.306", "",
               sprintf("MOTIF %s", name),
               sprintf("letter-probability matrix: alength= 4 w= %d nsites= %d E= 0",
                       nrow(mat), nsites),
               apply(mat, 1, function(r) paste(sprintf("%.6f", r), collapse = " "))), con)
  close(con)
}
L_IDX <- 1:9; R_IDX <- 10:15
write_meme(pm[L_IDX, , drop = FALSE], "M2_LEFT_EBOX",  file.path(OUT, "m2_left.meme"))
write_meme(pm[R_IDX, , drop = FALSE], "M2_RIGHT_ATAIL", file.path(OUT, "m2_right.meme"))

run_fimo <- function(meme, fa, tag, thresh = "1.0E-3") {
  f <- file.path(OUT, sprintf("fimo_%s.tsv", tag))
  cmd <- sprintf("fimo --verbosity 1 --bgfile --nrdb-- --thresh %s --text %s %s > %s 2>/dev/null",
                 thresh, shQuote(meme), shQuote(fa), shQuote(f))
  system(cmd)
  d <- tryCatch(read.delim(f, stringsAsFactors = FALSE), error = function(e) NULL)
  if (is.null(d) || !nrow(d)) return(NULL)
  d <- d[!startsWith(as.character(d[[1]]), "#") & nzchar(as.character(d[[1]])), ]
  d$start <- as.integer(d$start); d$stop <- as.integer(d$stop)
  d
}
Lh <- run_fimo(file.path(OUT, "m2_left.meme"),  JOINT_FA, "left")
Rh <- run_fimo(file.path(OUT, "m2_right.meme"), JOINT_FA, "right")
stopifnot(!is.null(Lh), !is.null(Rh))

# offset of every R hit relative to every L hit on the same strand & sequence
off <- do.call(rbind, lapply(split(seq_len(nrow(Lh)), Lh$sequence_name), function(ix) {
  s  <- Lh$sequence_name[ix[1]]
  rr <- Rh[Rh$sequence_name == s, , drop = FALSE]
  if (!nrow(rr)) return(NULL)
  do.call(rbind, lapply(ix, function(i) {
    same <- rr[rr$strand == Lh$strand[i], , drop = FALSE]
    if (!nrow(same)) return(NULL)
    d <- if (Lh$strand[i] == "+") same$start - Lh$stop[i] - 1L
         else Lh$start[i] - same$stop[i] - 1L
    d <- d[abs(d) <= 60]
    if (!length(d)) return(NULL)
    data.frame(seq = s, offset = d)
  }))
}))
tab <- table(factor(off$offset, levels = -60:60))
obs0 <- as.integer(tab[as.character(0)])
flank <- as.integer(tab[as.character(setdiff(-60:60, -3:3))])
res$R2 <- data.frame(
  n_L_hits = nrow(Lh), n_R_hits = nrow(Rh),
  pairs_within_60bp = nrow(off),
  at_register_0 = obs0,
  offcentre_mean = mean(flank), offcentre_sd = sd(flank),
  register_enrichment = obs0 / mean(flank),
  register_Z = (obs0 - mean(flank)) / sd(flank))
cat("\n== TEST R2  internal register ==\n"); print(res$R2, row.names = FALSE)
write.csv(data.frame(offset = -60:60, n = as.integer(tab)),
          file.path(OUT, "R2_register_profile.csv"), row.names = FALSE)

saveRDS(res, file.path(OUT, "motif2_robustness_partA.rds"))
cat("\nPart A written to ", OUT, "\n", sep = "")