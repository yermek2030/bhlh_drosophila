#!/usr/bin/env Rscript
# ============================================================
# icistarget_inputs_all_heterodimers.R
#
# Builds i-cisTarget region files for all four heterodimers (Trh/Tgo,
# Sim/Tgo, Dys/Tgo, Sima/Tgo), each split into Motif2-positive and
# Motif2-negative peaks. Six columns, column 6 in {+,-,?} (unstranded peaks use "?").
# Usage: Rscript 02_motif_analysis/icistarget/icistarget_inputs_all_heterodimers.R
# ============================================================

## Repository root: the environment variable PROJECT_DIR, or the working directory.
PROJECT_DIR <- normalizePath(Sys.getenv("PROJECT_DIR", unset = "."), mustWork = TRUE)
if (!file.exists(file.path(PROJECT_DIR, "README.md")))
  stop("PROJECT_DIR does not point at the repository root: ", PROJECT_DIR)

suppressPackageStartupMessages(library(GenomicRanges))

OUT <- file.path(PROJECT_DIR, "data", "icistarget", "all_heterodimers")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
A3 <- file.path(PROJECT_DIR, "data", "chapter3")
A4 <- file.path(PROJECT_DIR, "data", "chapter4")

MOTIF <- "STRCACTGARARAAA"
QMAX  <- 0.05

# peaks: BED (0-based) ; fimo: peak-local, sequence_name = "peak_N" indexing
# the same bed in row order. Trh/Tgo is taken from the canonical RDS so the
# 614 invariant is asserted rather than re-derived from a file.
het <- list(
  Trh_Tgo  = list(bed = NULL, rds = file.path(A3, "joint_peaks_filt.rds"),
                  fimo = file.path(A3, "fimo_joint", "fimo.tsv"), n = 614L),
  Sim_Tgo  = list(bed = file.path(A4, "Sim_Tgo_joint_peaks_filt.bed"),  rds = NULL,
                  fimo = file.path(A4, "fimo_Sim_Tgo", "fimo.tsv"),  n = 355L),
  Dys_Tgo  = list(bed = file.path(A4, "Dys_Tgo_joint_peaks_filt.bed"),  rds = NULL,
                  fimo = file.path(A4, "fimo_Dys_Tgo", "fimo.tsv"),  n = 870L),
  Sima_Tgo = list(bed = file.path(A4, "Sima_Tgo_joint_peaks_filt.bed"), rds = NULL,
                  fimo = file.path(A4, "fimo_Sima_Tgo", "fimo.tsv"), n = 296L)
)

write6 <- function(chrom, start0, end, tag, file) {
  d <- data.frame(chrom, start0, end,
                  name = sprintf("%s_%04d", tag, seq_along(chrom)),
                  score = 0L, strand = "?")
  stopifnot(all(d$start0 >= 0), all(d$end > d$start0))
  write.table(d, file.path(OUT, file), sep = "\t", quote = FALSE,
              row.names = FALSE, col.names = FALSE)
  nrow(d)
}

summary_rows <- list()
for (h in names(het)) {
  s <- het[[h]]
  if (!is.null(s$rds)) {
    gr <- readRDS(s$rds)
    chrom <- as.character(seqnames(gr)); st0 <- start(gr) - 1L; en <- end(gr)
  } else {
    b <- read.table(s$bed, sep = "\t", stringsAsFactors = FALSE)
    chrom <- b$V1; st0 <- b$V2; en <- b$V3
  }
  stopifnot(length(chrom) == s$n)                       # peak-count invariant

  fm <- read.delim(s$fimo, comment.char = "#", stringsAsFactors = FALSE)
  fm <- fm[!is.na(fm$start) & fm$motif_id == MOTIF & fm$q.value < QMAX, ]
  idx <- as.integer(sub("peak_", "", fm$sequence_name))
  stopifnot(!any(is.na(idx)), max(idx) <= s$n)          # index frame matches
  pos <- sort(unique(idx)); neg <- setdiff(seq_along(chrom), pos)
  stopifnot(length(pos) + length(neg) == s$n)

  np <- write6(chrom[pos], st0[pos], en[pos], paste0(h, "_m2pos"),
               sprintf("%s_motif2_positive_icistarget.bed", h))
  nn <- write6(chrom[neg], st0[neg], en[neg], paste0(h, "_m2neg"),
               sprintf("%s_motif2_negative_icistarget.bed", h))
  summary_rows[[h]] <- data.frame(heterodimer = h, peaks = s$n,
                                  motif2_hits = nrow(fm),
                                  positive = np, negative = nn,
                                  pct_positive = round(100 * np / s$n, 1))
}
sm <- do.call(rbind, summary_rows)
write.table(sm, file.path(OUT, "input_set_summary.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)
print(sm, row.names = FALSE)

# read back and assert the server's format rule on every file written
for (f in list.files(OUT, pattern = "\\.bed$", full.names = TRUE)) {
  x <- read.table(f, sep = "\t", stringsAsFactors = FALSE)
  stopifnot(ncol(x) == 6L, all(x$V6 == "?"))
}
cat("\nformat verified on", length(list.files(OUT, pattern = "\\.bed$")),
    "files: 6 columns, strand \"?\".\n")
cat("written to data/icistarget/all_heterodimers/\n")
