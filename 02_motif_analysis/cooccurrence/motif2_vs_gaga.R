#!/usr/bin/env Rscript
# ============================================================
# motif2_vs_gaga.R
# Tests whether Motif2 is distinct from GAGA/Trl motifs.
# T1 peak-level association; T2 instance-level overlap; T3 occurrence in GAGA-negative peaks.
# Usage: Rscript motif2_vs_gaga.R
# ============================================================

## Repository root: the environment variable PROJECT_DIR, or the working directory.
PROJECT_DIR <- normalizePath(Sys.getenv("PROJECT_DIR", unset = "."), mustWork = TRUE)
if (!file.exists(file.path(PROJECT_DIR, "README.md")))
  stop("PROJECT_DIR does not point at the repository root: ", PROJECT_DIR)

suppressPackageStartupMessages({
  library(GenomicRanges); library(Biostrings)
  library(BSgenome.Dmelanogaster.UCSC.dm6)
})
GENOME  <- BSgenome.Dmelanogaster.UCSC.dm6
OUT_DIR <- file.path(PROJECT_DIR, "data", "motif2_gaga")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ---- 1. canonical peaks + Motif2 instances --------------------------------
jp <- readRDS(file.path(PROJECT_DIR, "data/chapter3/joint_peaks_filt.rds"))
stopifnot(length(jp) == 614)

fimo <- read.table(file.path(PROJECT_DIR, "data/chapter3/fimo_joint/fimo.tsv"),
                   header = TRUE, sep = "\t", comment.char = "#")
fimo <- fimo[fimo$motif_id == "STRCACTGARARAAA" & fimo$q.value < 0.05, ]
stopifnot(nrow(fimo) > 0)
pidx <- as.integer(sub("peak_", "", fimo$sequence_name))
stopifnot(!any(is.na(pidx)), max(pidx) <= length(jp))
m2 <- GRanges(as.character(seqnames(jp))[pidx],
              IRanges(start(jp)[pidx] + fimo$start - 1L,
                      start(jp)[pidx] + fimo$stop  - 1L))
jp$Motif2 <- countOverlaps(jp, m2, ignore.strand = TRUE) > 0

# ---- 2. GAGA instances, via the matrices that actually scored --------------
# icistarget_289.meme carries the exact GAGA-type matrices enriched in these
# peaks, so the test uses those rather than a hand-picked substitute.
memef <- file.path(PROJECT_DIR, "data/icistarget/trh_tgo",
                   "icistarget_289.meme")
stopifnot(file.exists(memef))
want <- c("transfac_pro__M00723",                      # top Trl hit, NES 5.50
          "flyfactorsurvey__Trl_FlyReg_FBgn0013263",
          "homer__CTYTCTYTCTCTCTC_GAGA-repeat")
ml   <- readLines(memef, warn = FALSE)
keep <- grep("^MOTIF", ml, value = TRUE)
present <- want[vapply(want, function(w) any(grepl(w, keep, fixed = TRUE)), TRUE)]
stopifnot(length(present) > 0)
cat("GAGA matrices used:", paste(present, collapse = ", "), "\n")

# subset the MEME file to those motifs
hdr_end <- max(grep("^Background|^A 0\\.25", ml))
sub <- ml[1:hdr_end]
starts <- grep("^MOTIF", ml)
for (s in starts) {
  nm <- sub("^MOTIF\\s+", "", ml[s])
  if (!any(vapply(present, function(w) grepl(w, nm, fixed = TRUE), TRUE))) next
  e <- if (any(starts > s)) min(starts[starts > s]) - 1L else length(ml)
  sub <- c(sub, "", ml[s:e])
}
gaga_meme <- file.path(OUT_DIR, "gaga_matrices.meme")
writeLines(sub, gaga_meme)

# peak FASTA, named peak_N so FIMO coordinates convert the same way
seqs <- getSeq(GENOME, jp); names(seqs) <- paste0("peak_", seq_along(jp))
fa <- file.path(OUT_DIR, "joint_peaks.fa"); writeXStringSet(seqs, fa)

fimo_out <- file.path(OUT_DIR, "fimo_gaga")
st <- system2("fimo", c("--verbosity", "1", "--thresh", "1e-4",
                        "--oc", fimo_out, gaga_meme, fa),
              stdout = FALSE, stderr = FALSE)
stopifnot(st == 0, file.exists(file.path(fimo_out, "fimo.tsv")))
gf <- read.delim(file.path(fimo_out, "fimo.tsv"), comment.char = "#",
                 stringsAsFactors = FALSE)
gf <- gf[!is.na(gf$start) & gf$q.value < 0.05, ]
gidx <- as.integer(sub("peak_", "", gf$sequence_name))
gg <- GRanges(as.character(seqnames(jp))[gidx],
              IRanges(start(jp)[gidx] + gf$start - 1L,
                      start(jp)[gidx] + gf$stop  - 1L))
jp$GAGA <- countOverlaps(jp, gg, ignore.strand = TRUE) > 0
cat(sprintf("Motif2-positive peaks: %d/%d | GAGA-positive peaks: %d/%d\n",
            sum(jp$Motif2), length(jp), sum(jp$GAGA), length(jp)))

# ---- peak-level association ---------------------------------------------
tab <- table(Motif2 = factor(jp$Motif2, c(FALSE, TRUE)),
             GAGA   = factor(jp$GAGA,   c(FALSE, TRUE)))
ft  <- fisher.test(tab)
cat("\nT1 peak-level Motif2 x GAGA\n"); print(tab)
cat(sprintf("   OR = %.3f (%.3f-%.3f), p = %.4g\n", ft$estimate,
            ft$conf.int[1], ft$conf.int[2], ft$p.value))

# ---- Do the instances coincide? --------------------------------------------
ov <- sum(countOverlaps(m2, gg, ignore.strand = TRUE) > 0)
same <- findOverlaps(m2, jp, select = "first")
gpk  <- findOverlaps(gg, jp, select = "first")
d <- unlist(lapply(seq_along(m2), function(i) {
  j <- which(gpk == same[i]); if (!length(j)) return(NULL)
  min(abs(start(gg)[j] - start(m2)[i]))
}))
cat(sprintf("\nT2 instance level: %d Motif2 hits; %d (%.1f%%) directly overlap a GAGA hit\n",
            length(m2), ov, 100 * ov / length(m2)))
if (length(d)) cat(sprintf("   of %d co-occurring in one peak: median separation %.0f bp (range %.0f-%.0f)\n",
                           length(d), median(d), min(d), max(d)))

# ---- does Motif2 survive in GAGA-negative peaks? --------------------------
n_neg <- sum(!jp$GAGA); n_m2_neg <- sum(jp$Motif2 & !jp$GAGA)
cat(sprintf("\nT3 conditioning: %d/%d (%.1f%%) GAGA-NEGATIVE peaks carry Motif2\n",
            n_m2_neg, n_neg, 100 * n_m2_neg / max(n_neg, 1)))
cat(sprintf("   for comparison, %.1f%% of GAGA-positive peaks carry Motif2\n",
            100 * sum(jp$Motif2 & jp$GAGA) / max(sum(jp$GAGA), 1)))

res <- data.frame(
  test = c("peaks total","Motif2+ peaks","GAGA+ peaks","T1 OR","T1 p",
           "T2 Motif2 instances","T2 overlapping a GAGA hit","T2 pct overlapping",
           "T3 GAGA- peaks","T3 GAGA- peaks with Motif2","T3 pct",
           "T3 pct of GAGA+ peaks with Motif2"),
  value = c(length(jp), sum(jp$Motif2), sum(jp$GAGA),
            round(unname(ft$estimate),3), signif(ft$p.value,4),
            length(m2), ov, round(100*ov/length(m2),1),
            n_neg, n_m2_neg, round(100*n_m2_neg/max(n_neg,1),1),
            round(100*sum(jp$Motif2 & jp$GAGA)/max(sum(jp$GAGA),1),1)),
  stringsAsFactors = FALSE)
write.table(res, file.path(OUT_DIR, "motif2_vs_gaga_results.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("\nwritten:", file.path(OUT_DIR, "motif2_vs_gaga_results.tsv"), "\n")
