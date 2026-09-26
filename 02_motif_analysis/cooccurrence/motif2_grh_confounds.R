#!/usr/bin/env Rscript
# ============================================================
# motif2_grh_confounds.R: tests three explanations for Grh motif enrichment in Motif2 peaks (Trh/Tgo).
# C1 proximity: checks Motif2 vs promoter-proximity association, repeats Motif2/Grh contrast in distal peaks only.
# C2 sequence: checks whether Motif2 hits overlap Grh hits at instance level

## Repository root: the environment variable PROJECT_DIR, or the working directory.
PROJECT_DIR <- normalizePath(Sys.getenv("PROJECT_DIR", unset = "."), mustWork = TRUE)
if (!file.exists(file.path(PROJECT_DIR, "README.md")))
  stop("PROJECT_DIR does not point at the repository root: ", PROJECT_DIR)

suppressPackageStartupMessages({
  library(GenomicRanges); library(ChIPseeker)
  library(TxDb.Dmelanogaster.UCSC.dm6.ensGene)
})
OUT <- file.path(PROJECT_DIR, "data", "motif2_grh_confounds")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
A3 <- file.path(PROJECT_DIR, "data", "chapter3")

jp <- readRDS(file.path(A3, "joint_peaks_filt.rds"))
stopifnot(length(jp) == 614)                       # canonical invariant

read_fimo <- function(f, motif = NULL, q = 0.05) {
  x <- read.table(f, header = TRUE, sep = "\t", comment.char = "#",
                  stringsAsFactors = FALSE)
  x <- x[!is.na(x$start), ]
  if (!is.null(motif)) x <- x[x$motif_id %in% motif, ]
  x[x$q.value < q, ]
}
m2 <- read_fimo(file.path(A3, "fimo_joint", "fimo.tsv"), "STRCACTGARARAAA")
m2$peak <- as.integer(sub("peak_", "", m2$sequence_name))
jp$M2 <- seq_along(jp) %in% sort(unique(m2$peak))

# ---- proximity -------------------------------------------------------------
an <- as.data.frame(annotatePeak(jp, TxDb = TxDb.Dmelanogaster.UCSC.dm6.ensGene,
                                 tssRegion = c(-1000, 1000), verbose = FALSE))
jp$promoter <- grepl("^Promoter", an$annotation)
t1 <- table(Motif2 = jp$M2, Promoter = jp$promoter)
f1 <- fisher.test(t1)
cat("=== C1  Motif2 versus promoter-proximity ===\n"); print(t1)
cat(sprintf("  OR = %.2f (%.2f-%.2f), p = %.3g\n",
            f1$estimate, f1$conf.int[1], f1$conf.int[2], f1$p.value))
cat(sprintf("  distal fraction: Motif2+ %.1f%% vs Motif2- %.1f%%\n\n",
            100 * mean(!jp$promoter[jp$M2]), 100 * mean(!jp$promoter[!jp$M2])))

# ---- does Motif2 physically contain a Grh site? ------------------------
grh_f <- file.path(OUT, "fimo_grh", "fimo.tsv")
if (file.exists(grh_f)) {
  gh <- read_fimo(grh_f)
  gh$peak <- as.integer(sub("peak_", "", gh$sequence_name))
  jp$GRH <- seq_along(jp) %in% sort(unique(gh$peak))
  gr_m2 <- GRanges(m2$sequence_name, IRanges(m2$start, m2$stop))
  gr_gh <- GRanges(gh$sequence_name, IRanges(gh$start, gh$stop))
  ov <- sum(countOverlaps(gr_m2, gr_gh) > 0)
  cat("=== C2  instance-level overlap of Motif2 with Grh matches ===\n")
  cat(sprintf("  Motif2 instances: %d | Grh instances: %d\n", nrow(m2), nrow(gh)))
  cat(sprintf("  Motif2 instances overlapping a Grh match: %d (%.1f%%)\n",
              ov, 100 * ov / nrow(m2)))
  t2 <- table(Motif2 = jp$M2, Grh = jp$GRH); f2 <- fisher.test(t2)
  cat("  peak level:\n"); print(t2)
  cat(sprintf("  OR = %.2f, p = %.3g\n", f2$estimate, f2$p.value))
  d <- jp[!jp$promoter]
  t3 <- table(Motif2 = d$M2, Grh = d$GRH); f3 <- fisher.test(t3)
  cat(sprintf("\n  DISTAL peaks only (n = %d): OR = %.2f, p = %.3g\n",
              length(d), f3$estimate, f3$p.value))
  res <- data.frame(
    quantity = c("peaks", "Motif2+ peaks", "distal OR Motif2 x promoter",
                 "distal % Motif2+", "distal % Motif2-",
                 "Motif2 instances", "Grh instances",
                 "Motif2 instances overlapping Grh", "pct overlapping",
                 "peak-level OR Motif2 x Grh", "peak-level p",
                 "distal-only OR", "distal-only p"),
    value = c(length(jp), sum(jp$M2), round(unname(f1$estimate), 3),
              round(100 * mean(!jp$promoter[jp$M2]), 1),
              round(100 * mean(!jp$promoter[!jp$M2]), 1),
              nrow(m2), nrow(gh), ov, round(100 * ov / nrow(m2), 1),
              round(unname(f2$estimate), 3), signif(f2$p.value, 4),
              round(unname(f3$estimate), 3), signif(f3$p.value, 4)))
  write.table(res, file.path(OUT, "motif2_grh_confound_results.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  cat("\nwritten:", file.path(OUT, "motif2_grh_confound_results.tsv"), "\n")
} else {
  cat("=== C2 skipped: run the FIMO step first ===\n  expected:", grh_f, "\n")
}
