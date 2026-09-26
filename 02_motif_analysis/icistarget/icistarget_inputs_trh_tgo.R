#!/usr/bin/env Rscript
# ============================================================
# icistarget_inputs_trh_tgo.R
# Writes i-cisTarget region files, splitting 614 Trh/Tgo joint peaks into
# Motif2-positive (178) and Motif2-negative (436) sets.
# Six-column BED, strand column "?" (rtracklayer export gives "*").
# Usage: Rscript 02_motif_analysis/icistarget/icistarget_inputs_trh_tgo.R
# ============================================================

## Repository root: the environment variable PROJECT_DIR, or the working directory.
PROJECT_DIR <- normalizePath(Sys.getenv("PROJECT_DIR", unset = "."), mustWork = TRUE)
if (!file.exists(file.path(PROJECT_DIR, "README.md")))
  stop("PROJECT_DIR does not point at the repository root: ", PROJECT_DIR)

suppressPackageStartupMessages(library(GenomicRanges))

OUT <- file.path(PROJECT_DIR, "data", "icistarget", "trh_tgo")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

# ---- canonical peak set: the 614 invariant --------------------------------
jp <- readRDS(file.path(PROJECT_DIR, "data/chapter3/joint_peaks_filt.rds"))
stopifnot(length(jp) == 614)

# ---- Motif2 status. fimo_joint/fimo.tsv is PEAK-LOCAL (sequence_name =
#      "peak_N"); only the peak INDEX is needed here, not the coordinates.
fm <- read.table(file.path(PROJECT_DIR, "data/chapter3/fimo_joint/fimo.tsv"),
                 header = TRUE, sep = "\t", comment.char = "#")
fm <- fm[fm$motif_id == "STRCACTGARARAAA" & fm$q.value < 0.05, ]
stopifnot(nrow(fm) > 0)
i   <- as.integer(sub("peak_", "", fm$sequence_name))
stopifnot(!any(is.na(i)), max(i) <= length(jp))
pos <- sort(unique(i))
neg <- setdiff(seq_along(jp), pos)
stopifnot(length(pos) + length(neg) == 614L, length(pos) == 178L)

# ---- writer ----------------------------------------------------------------
write_icistarget_bed <- function(idx, tag, file) {
  d <- data.frame(chrom  = as.character(seqnames(jp))[idx],
                  start  = start(jp)[idx] - 1L,        # GRanges 1-based -> BED 0-based
                  end    = end(jp)[idx],
                  name   = sprintf("%s_%03d", tag, seq_along(idx)),
                  score  = 0L,
                  strand = "?")
  stopifnot(all(d$start >= 0), all(d$end > d$start))
  write.table(d, file.path(OUT, file), sep = "\t", quote = FALSE,
              row.names = FALSE, col.names = FALSE)
  cat(sprintf("%-44s %3d regions  %6d bp\n", file, nrow(d), sum(d$end - d$start)))
  invisible(d)
}
a <- write_icistarget_bed(pos, "m2pos", "motif2_positive_peaks_icistarget.bed")
b <- write_icistarget_bed(neg, "m2neg", "motif2_negative_peaks_icistarget.bed")

# ---- read back and assert the server's format rule -------------------------
for (f in c("motif2_positive_peaks_icistarget.bed",
            "motif2_negative_peaks_icistarget.bed")) {
  x <- read.table(file.path(OUT, f), sep = "\t", stringsAsFactors = FALSE)
  stopifnot(ncol(x) == 6L, all(x$V6 %in% c("+", "-", "?")))
}
cat("\nformat verified: 6 columns, strand in {+,-,?} -- i-cisTarget will accept these.\n")
cat("submit:", file.path("data/icistarget/trh_tgo",
                         "motif2_positive_peaks_icistarget.bed"), "\n")
