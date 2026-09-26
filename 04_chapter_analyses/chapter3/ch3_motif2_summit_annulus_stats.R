## ============================================================
##  ch3_motif2_summit_annulus_stats.R
##  Summit core (|d| <= 25 bp) and annulus (50 <= |d| <= 100 bp) occupancy
##  vs uniform null, plus left/right symmetry, for Motif2 instances.
##  Input: data/motif2_robustness/R11b_summit_distances.csv
##  Output: 14 values in results/reported_values.tsv
## ============================================================

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


RB <- file.path(PROJECT_DIR, "data", "motif2_robustness")
d  <- read.csv(file.path(RB, "R11b_summit_distances.csv"))$dist

N <- length(d)
W <- max(abs(d))            # scanned half-window, +/- bp

## ---- Band definitions ----
CORE_HALF <- 25L            # summit core, |d| <= 25
RING_LO   <- 50L            # annulus, closed: RING_LO <= |d| <= RING_HI
RING_HI   <- 100L

in_core <- abs(d) <= CORE_HALF
in_ring <- abs(d) >= RING_LO & abs(d) <= RING_HI

## Uniform null over the scanned window [-W, W] (2W+1 positions).
core_null <- (2 * CORE_HALF + 1) / (2 * W + 1)
ring_null <- 2 * (RING_HI - RING_LO + 1) / (2 * W + 1)

odds <- function(p) p / (1 - p)
core_pct <- mean(in_core); ring_pct <- mean(in_ring)
core_or  <- odds(core_pct) / odds(core_null)
ring_or  <- odds(ring_pct) / odds(ring_null)
core_p   <- binom.test(sum(in_core), N, core_null)$p.value
ring_p   <- binom.test(sum(in_ring), N, ring_null)$p.value

## Symmetry within the same band.
symL <- sum(in_ring & d < 0)
symR <- sum(in_ring & d > 0)
stopifnot(symL + symR == sum(in_ring))          # band has no d == 0 member
sym_p <- binom.test(symL, symL + symR, 0.5)$p.value

## Symmetry across all instances, not just the band.
all_p <- binom.test(sum(d < 0), sum(d != 0), 0.5)$p.value

## Check the direction asserted in Section 3.7: core depletion and ring enrichment relative to the null.
stopifnot(core_pct < core_null, ring_pct > ring_null)

fmt_p <- function(p, digits = 2) {
  if (is.na(p)) return("NA")
  s <- formatC(p, format = "e", digits = digits)
  q <- strsplit(s, "e")[[1]]
  paste0(q[1], "e", as.integer(q[2]))
}

vals <- list(
  motifTwoSummitN        = N,
  motifTwoSummitWindow   = W,
  motifTwoSummitCorePct  = round(100 * core_pct,  1),
  motifTwoSummitCoreNull = round(100 * core_null, 1),
  motifTwoSummitCoreOR   = round(core_or, 2),
  motifTwoSummitCoreP    = fmt_p(core_p),
  motifTwoRingPct        = round(100 * ring_pct,  1),
  motifTwoRingNull       = round(100 * ring_null, 1),
  motifTwoRingOR         = round(ring_or, 2),
  motifTwoRingP          = fmt_p(ring_p),
  motifTwoRingSymL       = symL,
  motifTwoRingSymR       = symR,
  motifTwoRingSymP       = sprintf("%.3f", sym_p),
  motifTwoAllSymP        = sprintf("%.3f", all_p)
)

save_values("ch3_summit_centrality", vals)

cat(sprintf("\nband: %d <= |d| <= %d  (closed)\n", RING_LO, RING_HI))
cat(sprintf("n=%d  window=+/-%d\n", N, W))
cat(sprintf("core : %.1f%% vs %.1f%% null  OR=%.2f  p=%s\n",
            100*core_pct, 100*core_null, core_or, fmt_p(core_p)))
cat(sprintf("ring : %.1f%% vs %.1f%% null  OR=%.2f  p=%s\n",
            100*ring_pct, 100*ring_null, ring_or, fmt_p(ring_p)))
cat(sprintf("sym  : L=%d R=%d (sum %d, = %.1f%% of n)  p=%.3f\n",
            symL, symR, symL + symR, 100*(symL + symR)/N, sym_p))
cat(sprintf("all  : %d neg / %d pos  p=%.3f\n", sum(d < 0), sum(d > 0), all_p))

