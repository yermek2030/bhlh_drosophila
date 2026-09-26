#!/usr/bin/env Rscript
# ============================================================
# Tests whether the claim holds using only the 9-bp core PWM: reruns or = 1.28, p = 3.5e-3 Fisher test with core PWM in place of 15-bp PWM; TOMTOMs core against the four de novo motifs and redoes 12 pairwise comparisons on core-trimmed PWMs.
# TOMTOM offset: target_index = query_index + Optimal_offset (1-based).
# Depends: fimo, tomtom, GenomicRanges. Prereq: 02_split_scan_register.R (writes m2_core9.meme).
## Repository root: PROJECT_DIR env var, or working directory.
# ============================================================
PROJECT_DIR <- normalizePath(Sys.getenv("PROJECT_DIR", unset = "."), mustWork = TRUE)
if (!file.exists(file.path(PROJECT_DIR, "README.md")))
  stop("PROJECT_DIR does not point at the repository root: ", PROJECT_DIR)

DATA_DIR      <- file.path(PROJECT_DIR, "data")
setwd(PROJECT_DIR)

suppressPackageStartupMessages({
  library(GenomicRanges); library(GenomeInfoDb)   # keepSeqlevels lives here
})

AUX  <- DATA_DIR
OUT  <- file.path(AUX, "motif2_concat_audit")
CH3  <- file.path(AUX, "chapter3")
CH4  <- file.path(AUX, "chapter4")
REFS <- file.path(PROJECT_DIR, "data/reference_gene_lists")

f_core   <- file.path(OUT, "m2_core9.meme")
joint_fa <- file.path(CH3, "joint.fa")                 # 614 peaks, peak_N headers
enh_fa   <- file.path(CH3, "enhancers_genomewide.fa")  # genomic headers
four     <- file.path(CH4, "motif2_4way_masked.meme")
stopifnot(file.exists(f_core), file.exists(joint_fa),
          file.exists(enh_fa),  file.exists(four))

PVAL <- 1e-4   # identical threshold for both sets, this is the like-for-like

# --- MEME Suite runner: prefer `conda run -n phd`, fall back to PATH ----------
meme_run <- function(tool, args) {
  if (nzchar(Sys.which("conda")))
    return(system2("conda", c("run", "-n", "phd", tool, args)))
  if (nzchar(Sys.which(tool)))
    return(system2(tool, args))
  stop("Cannot find ", tool, ": neither `conda` nor `", tool, "` is on PATH.\n",
       "  Run this from a shell where `conda activate phd` works.")
}

# ------------------------------------------------------------ (a) core-only or
run_fimo <- function(meme, fa, dir) {
  if (!file.exists(file.path(dir, "fimo.tsv")))
    meme_run("fimo", c("--oc", shQuote(dir), "--verbosity", "1",
                       "--bgfile", "--nrdb--", "--thresh", "1.0E-4",
                       "--max-stored-scores", "10000000",
                       shQuote(meme), shQuote(fa)))
  f <- read.table(file.path(dir, "fimo.tsv"), header = TRUE, sep = "\t",
                  comment.char = "#", quote = "", stringsAsFactors = FALSE)
  names(f)[names(f) == "p.value"] <- "pval"
  f[f$pval < PVAL, ]
}

fj <- run_fimo(f_core, joint_fa, file.path(OUT, "fimo_core9_joint"))
fe <- run_fimo(f_core, enh_fa,   file.path(OUT, "fimo_core9_enhancers"))

n_joint_pos <- length(unique(fj$sequence_name))   # peak_N headers
n_joint_tot <- 614L
stopifnot(n_joint_pos <= n_joint_tot)

# FIMO collapses chr:start-end headers to bare chromosome names and reports
# ABSOLUTE coordinates, count enhancers by genomic OVERLAP, never by name.
enh_raw <- read.table(file.path(REFS, "merged_embryo_5-20_enhancers_dm6.bed"),
                      header = FALSE, stringsAsFactors = FALSE)[, 1:3]
enh_gr  <- GRanges(enh_raw[[1]], IRanges(enh_raw[[2]] + 1L, enh_raw[[3]]))
enh_gr  <- keepSeqlevels(enh_gr, c("chr2L","chr2R","chr3L","chr3R","chr4","chrX"),
                         pruning.mode = "coarse")
hits_gr <- GRanges(fe$sequence_name, IRanges(fe$start, fe$stop))

or_for_width <- function(min_w, tag) {
  e <- enh_gr[width(enh_gr) >= min_w]
  pos <- sum(countOverlaps(e, hits_gr, ignore.strand = TRUE) > 0)
  ft  <- fisher.test(matrix(c(n_joint_pos, n_joint_tot - n_joint_pos,
                              pos,         length(e) - pos), nrow = 2))
  data.frame(pwm = "core_9bp_STRCACTGA", enh_width_filter = tag,
             joint_pos = n_joint_pos, joint_tot = n_joint_tot,
             joint_pct = round(100 * n_joint_pos / n_joint_tot, 1),
             enh_pos = pos, enh_tot = length(e),
             enh_pct = round(100 * pos / length(e), 1),
             OR = round(unname(ft$estimate), 2),
             CI_lo = round(ft$conf.int[1], 2), CI_hi = round(ft$conf.int[2], 2),
             p = signif(ft$p.value, 3))
}
# width >= 15 reproduces the EXACT denominator of the published or = 1.28 test;
# width >= 9 is the denominator appropriate to a 9-bp PWM.
or_tab <- rbind(or_for_width(15L, "ge15_matches_published"),
                or_for_width( 9L, "ge9_pwm_width"))
cat("\n=== (a) enhancer enrichment, CORE ONLY ===\n"); print(or_tab, row.names = FALSE)
cat("Benchmark -- full 15-bp PWM: OR = 1.28, p = 3.5e-3\n")
ref <- or_tab[or_tab$enh_width_filter == "ge15_matches_published", ]
cat(if (ref$OR >= 1.15 && ref$p < 0.05)
      "SURVIVES: the enrichment is carried by the core; the tail is dispensable.\n"
    else
      "DOES NOT SURVIVE on the core alone: the 15-bp form carries the signal.\n")
write.csv(or_tab, file.path(OUT, "T3_core_only_enhancer_OR.csv"), row.names = FALSE)

# ----------------------------------------------- (b) core-only PWM equivalence
tt_dir <- file.path(OUT, "tomtom_core9_vs_4way")
meme_run("tomtom", c("-no-ssc", "-oc", shQuote(tt_dir), "-verbosity", "1",
                     "-min-overlap", "9", "-dist", "pearson", "-thresh", "1",
                     shQuote(f_core), shQuote(four)))
tt <- read.table(file.path(tt_dir, "tomtom.tsv"), header = TRUE, sep = "\t",
                 comment.char = "#", quote = "", stringsAsFactors = FALSE)
names(tt)[names(tt) == "q.value"] <- "qval"
tt <- tt[!is.na(tt$Target_ID) & tt$Target_ID != "", ]
core_tt <- tt[, c("Query_ID", "Target_ID", "Optimal_offset", "p.value", "qval",
                  "Overlap", "Target_consensus", "Orientation")]
cat("\n=== (b) TOMTOM: 9-bp core vs the four de novo motifs ===\n")
print(core_tt, row.names = FALSE)
cat("Worst q across the four:", signif(max(core_tt$qval), 3), "\n")
cat("CAVEAT: a poor q here can be a LENGTH artefact, not dissimilarity -- a de novo\n",
    " motif that starts inside the core cannot align a full 9-bp query. The\n",
    " length-matched comparison below is the one to report.\n")
write.csv(core_tt, file.path(OUT, "T3_tomtom_core9_vs_4way.csv"), row.names = FALSE)

# ---- 12 pairwise comparisons among the four motifs TRIMMED to the core -------
# Trimming is ANCHORED ON the CACTG BOX, not on the TOMTOM offset. The de novo
# consensuses do not all start at the same point (Sim = TGCACTGA..., i.e. one
# base into the core), so an offset-based window runs off the front of Sim and
# would silently drop it. Anchoring on CACTG is exact and orientation-aware.
read_all_meme <- function(f) {
  x  <- readLines(f, warn = FALSE)
  s  <- grep("^MOTIF", x)
  lp <- grep("^letter-probability matrix", x)
  lapply(seq_along(s), function(k) {
    id <- strsplit(trimws(x[s[k]]), "[[:space:]]+")[[1]][2]
    i  <- lp[lp > s[k]][1]
    w  <- as.integer(sub(".*w= *([0-9]+).*", "\\1", x[i]))
    m  <- do.call(rbind, lapply(strsplit(trimws(x[(i + 1):(i + w)]),
                                         "[[:space:]]+"), as.numeric))
    list(id = id, m = m[, 1:4, drop = FALSE])
  })
}
revcomp_pwm <- function(m) m[nrow(m):1, 4:1, drop = FALSE]
BG <- c(.306, .194, .194, .306)
consensus_of <- function(m) paste(c("A","C","G","T")[max.col(m, "first")], collapse = "")

# canonical 9-bp core, as written by 02_split_scan_register.R
p_core <- read_all_meme(f_core)[[1]]$m
stopifnot(nrow(p_core) == 9L)

# In the canonical PWM the core is STRCACTGA: 3 bases, then CACTG, then A.
# The de novo motifs do not all reach 3 bases upstream of CACTG (Sim starts at
# T, i.e. only 2), so the trim window is the LARGEST WINDOW all FOUR SUPPLY.
# This is reported explicitly rather than padded, and the query core is
# trimmed identically so the comparison stays like-for-like.
CORE_ANCHOR <- "CACTG"; UP_WANT <- 3L; DOWN_WANT <- 1L

mots <- read_all_meme(four)
anchor_of <- function(m0) {
  for (ori in c("+", "-")) {
    m   <- if (ori == "-") revcomp_pwm(m0) else m0
    pos <- regexpr(CORE_ANCHOR, consensus_of(m), fixed = TRUE)
    if (pos > 0) return(list(m = m, ori = ori, anchor = as.integer(pos)))
  }
  NULL
}
anc <- lapply(mots, function(mo) anchor_of(mo$m))
names(anc) <- vapply(mots, `[[`, "", "id")
missing <- names(anc)[vapply(anc, is.null, logical(1))]
if (length(missing)) cat("no CACTG anchor in:", paste(missing, collapse = ", "), "\n")
anc <- anc[!vapply(anc, is.null, logical(1))]
stopifnot(length(anc) == 4L)

up_use   <- min(UP_WANT,   min(vapply(anc, function(a) a$anchor - 1L, 0L)))
down_use <- min(DOWN_WANT, min(vapply(anc, function(a) nrow(a$m) - (a$anchor + 4L), 0L)))
CORE_W   <- up_use + 5L + down_use
cat(sprintf("\nLargest core window all four de novo motifs supply: %d bp (%d up + CACTG + %d down)\n",
            CORE_W, up_use, down_use))
if (CORE_W < UP_WANT + 5L + DOWN_WANT)
  cat("NOTE: shorter than the canonical 9-bp core -- the de novo motifs do not all\n",
      "      extend far enough 5' of CACTG. Reported, not padded.\n")

hdr <- c("MEME version 5", "", "ALPHABET= ACGT", "", "strands: + -", "",
         "Background letter frequencies (from unknown source):",
         sprintf("A %.3f C %.3f G %.3f T %.3f", BG[1], BG[2], BG[3], BG[4]), "")
emit <- function(m, id) c(paste("MOTIF", id),
  sprintf("letter-probability matrix: alength= 4 w= %d nsites= 36 E= 0", nrow(m)),
  apply(m, 1, function(z) paste(sprintf("%.6f", z), collapse = " ")), "")

body <- character(0); kept <- character(0); trim_log <- list()
for (id in names(anc)) {
  a   <- anc[[id]]
  idx <- seq.int(a$anchor - up_use, a$anchor + 4L + down_use)
  sub <- a$m[idx, , drop = FALSE]
  body <- c(body, emit(sub, paste0(id, "_core")))
  kept <- c(kept, id)
  trim_log[[length(trim_log) + 1L]] <- data.frame(
    motif = id, orientation = a$ori, anchor_pos = a$anchor,
    window = sprintf("%d-%d of %d", min(idx), max(idx), nrow(a$m)),
    full_consensus = consensus_of(a$m), core_consensus = consensus_of(sub))
}
f_trim <- file.path(OUT, "motif2_4way_core_trimmed.meme")
writeLines(c(hdr, body), f_trim)

# Trim the CANONICAL core query the same way, so like is compared with like.
core_q  <- p_core[seq.int(4L - up_use, 8L + down_use), , drop = FALSE]
f_core_t <- file.path(OUT, "m2_core_trimmed_query.meme")
writeLines(c(hdr, emit(core_q, "M2_CORE_TRIMMED")), f_core_t)

cat("\n--- core-trimming log ---\n")
trim_df <- do.call(rbind, trim_log); print(trim_df, row.names = FALSE)
write.csv(trim_df, file.path(OUT, "T3_core_trim_log.csv"), row.names = FALSE)
stopifnot(length(kept) == 4L)

# Re-run the query-vs-four comparison with the LENGTH-MATCHED core query.
# The 9-bp query above under-matches Sim (q ~ 0.88) purely because Sim's de
# novo motif starts one base inside the core, so a 9-bp query cannot align
# fully to it. With both sides trimmed to the common window the comparison
# is like-for-like; report this one, and the 9-bp run as the caveat.
ttm_dir <- file.path(OUT, "tomtom_coreTrimmed_vs_4way")
meme_run("tomtom", c("-no-ssc", "-oc", shQuote(ttm_dir), "-verbosity", "1",
                     "-min-overlap", as.character(CORE_W), "-dist", "pearson",
                     "-thresh", "1", shQuote(f_core_t), shQuote(four)))
ttm <- read.table(file.path(ttm_dir, "tomtom.tsv"), header = TRUE, sep = "\t",
                  comment.char = "#", quote = "", stringsAsFactors = FALSE)
names(ttm)[names(ttm) == "q.value"] <- "qval"
ttm <- ttm[!is.na(ttm$Target_ID) & ttm$Target_ID != "", ]
cat(sprintf("\n=== canonical core (%d bp, length-matched) vs the four de novo motifs ===\n",
            CORE_W))
print(ttm[, c("Query_ID", "Target_ID", "p.value", "qval", "Overlap",
              "Target_consensus", "Orientation")], row.names = FALSE)
cat("Worst q:", signif(max(ttm$qval), 3), "\n")
write.csv(ttm, file.path(OUT, "T3_tomtom_coreTrimmed_vs_4way.csv"), row.names = FALSE)

tt2_dir <- file.path(OUT, "tomtom_4way_core_pairwise")
meme_run("tomtom", c("-no-ssc", "-oc", shQuote(tt2_dir), "-verbosity", "1",
                     "-min-overlap", as.character(CORE_W), "-dist", "pearson",
                     "-thresh", "1", shQuote(f_trim), shQuote(f_trim)))
tt2 <- read.table(file.path(tt2_dir, "tomtom.tsv"), header = TRUE, sep = "\t",
                  comment.char = "#", quote = "", stringsAsFactors = FALSE)
names(tt2)[names(tt2) == "q.value"] <- "qval"
tt2 <- tt2[!is.na(tt2$Target_ID) & tt2$Target_ID != "" &
             tt2$Query_ID != tt2$Target_ID, ]
cat("\n=== 12 pairwise comparisons, core-trimmed ===\n")
print(tt2[, c("Query_ID", "Target_ID", "p.value", "qval", "Overlap")], row.names = FALSE)
cat("n comparisons:", nrow(tt2), "  worst q:", signif(max(tt2$qval), 3), "\n")
if (nrow(tt2) != 12L)
  cat("WARNING: expected 12 off-diagonal comparisons, got", nrow(tt2), "\n")
write.csv(tt2, file.path(OUT, "T3_tomtom_4way_core_pairwise.csv"), row.names = FALSE)

cat("\nWrote T3_* to", OUT, "\n")
