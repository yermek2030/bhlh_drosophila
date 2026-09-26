## ====
##  motif2_summit_distances.R - writes R11b_summit_distances.csv
##  dist = (start+stop)/2 - 250 from fimo_truesummits.tsv, 380 rows
##  Upstream motif2_summit_centrality.R; cons

## Repository root: the environment variable PROJECT_DIR, or the working directory.
PROJECT_DIR <- normalizePath(Sys.getenv("PROJECT_DIR", unset = "."), mustWork = TRUE)
if (!file.exists(file.path(PROJECT_DIR, "README.md")))
  stop("PROJECT_DIR does not point at the repository root: ", PROJECT_DIR)

RB <- file.path(PROJECT_DIR, "data", "motif2_robustness")
stopifnot(dir.exists(RB))

FA        <- file.path(RB, "trh_TRUEsummits_500bp.fa")
FIMO_TSV  <- file.path(RB, "fimo_truesummits.tsv")
OUT_CSV   <- file.path(RB, "R11b_summit_distances.csv")
MEME_M2   <- file.path(PROJECT_DIR, "data", "chapter3",
                       "STRCACTGARARAAA.meme")

## Window geometry, as defined in motif2_summit_centrality.R:
##   cen <- GRanges(seqnames, IRanges(ctr - 249L, ctr + 250L))
## so the summit sits at offset 250 in 1-based within-sequence coordinates.
WINDOW_HALF   <- 250L
SUMMIT_OFFSET <- 250L
FIMO_THRESH   <- 1e-4          # FIMO default; observed max p = 9.65e-05

RERUN_FIMO <- FALSE            # see origin note above

## ---- 1. FIMO hits in the summit windows ----
if (RERUN_FIMO || !file.exists(FIMO_TSV)) {
  stopifnot(file.exists(FA), file.exists(MEME_M2))
  MEME_BIN <- Sys.getenv("MEME_BIN", "")
  if (nzchar(MEME_BIN))
    Sys.setenv(PATH = paste(MEME_BIN, Sys.getenv("PATH"), sep = ":"))
  if (nzchar(Sys.which("fimo")) == FALSE)
    stop("fimo not on PATH; set MEME_BIN to the MEME 5.5.5 bin directory")
  tmp <- file.path(RB, "fimo_truesummits_rerun")
  system(sprintf("fimo --verbosity 1 --thresh %g --oc %s %s %s",
                 FIMO_THRESH, shQuote(tmp), shQuote(MEME_M2), shQuote(FA)))
  file.copy(file.path(tmp, "fimo.tsv"), FIMO_TSV, overwrite = TRUE)
  message("FIMO re-run: ", FIMO_TSV,
          " -- re-check \\motifTwoSummitN before shipping.")
}
stopifnot(file.exists(FIMO_TSV))

fi <- read.delim(FIMO_TSV, comment.char = "#", stringsAsFactors = FALSE)
fi <- fi[nzchar(as.character(fi$sequence_name)), ]          # drop FIMO trailer
stopifnot(nrow(fi) > 0L, all(fi$motif_id == "STRCACTGARARAAA"))

## Both strands are retained. Motif2 is scored as a PWM on
## double-stranded genomic DNA; restricting to the plus strand would halve
## the sample and impose an orientation the element does not have.
stopifnot(all(fi$strand %in% c("+", "-")))

## ---- 2. Signed distance from motif centre to the ChIP summit ----
## Centre, not start: the element is 15 bp and its start is strand-dependent,
## so a start-referenced distance would displace minus-strand hits by the
## motif width and smear the annulus by 15 bp.
mid  <- (fi$start + fi$stop) / 2
dist <- as.integer(round(mid - SUMMIT_OFFSET))

out <- data.frame(dist = dist)

## ---- 3. Verify against the shipped file before overwriting ----
## The shipped CSV is the object the reported values were computed from.
## If this script cannot reproduce it exactly, the discrepancy is resolved
## before the file is overwritten.
if (file.exists(OUT_CSV)) {
  ref <- read.csv(OUT_CSV)$dist
  if (!identical(as.integer(ref), dist)) {
    stop("recomputed distances differ from the shipped ",
         basename(OUT_CSV), " (n_ref = ", length(ref),
         ", n_new = ", length(dist), "). Resolve before overwriting.")
  }
  message("verified: recomputed distances identical to shipped CSV (n = ",
          length(dist), ")")
}

write.csv(out, OUT_CSV, row.names = FALSE)

## ---- 4. Report the values recorded in results/reported_values.tsv ----
W        <- max(abs(dist))
core_pct <- 100 * sum(abs(dist) <= 25) / length(dist)
core_nul <- 100 * 51 / (2 * W + 1)

cat(sprintf("\nwritten: %s\n", OUT_CSV))
cat(sprintf("  n instances            = %d      (\\motifTwoSummitN)\n",
            length(dist)))
cat(sprintf("  scanned half-window    = %d bp   (\\motifTwoSummitWindow)\n", W))
cat(sprintf("  within +/-25 bp        = %.1f%%   (\\motifTwoSummitCorePct)\n",
            core_pct))
cat(sprintf("  uniform null expects   = %.1f%%   (\\motifTwoSummitCoreNull)\n",
            core_nul))
