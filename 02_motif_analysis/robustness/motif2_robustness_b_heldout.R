# ============================================================
# motif2_robustness_b_heldout.R
#   R3 half vs whole information, A-rich tail discrimination
#   R4 split-half generalisation, discover on A, score held-out B
#   R5 HOT/hyper-ChIPability, promiscuous loci check; R6 positional concentration vs peak centre
# Run: MEME_BIN=... Rscript motif2_robustness_b_heldout.R
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

fa   <- readDNAStringSet(JOINT_FA)
peaks <- import(JOINT_BED, format = "BED")
stopifnot(length(fa) == 614L, length(peaks) == 614L)

resB <- list()
set.seed(20260803)

fimo_txt <- function(meme, fa_path, tag, thresh = "1.0E-4") {
  f <- file.path(OUT, sprintf("fimoB_%s.tsv", tag))
  system(sprintf("fimo --verbosity 1 --bgfile --nrdb-- --thresh %s --text %s %s > %s 2>/dev/null",
                 thresh, shQuote(meme), shQuote(fa_path), shQuote(f)))
  d <- tryCatch(read.delim(f, stringsAsFactors = FALSE), error = function(e) NULL)
  if (is.null(d) || !nrow(d)) return(d)
  d[!startsWith(as.character(d[[1]]), "#") & nzchar(as.character(d[[1]])), , drop = FALSE]
}

# ------------------------------------------------------------
# TEST R3, does the A-rich tail add discrimination over the E-box core alone?
# Enrichment of each sub-motif in real peaks vs 50x dinucleotide-shuffled peaks.
# ------------------------------------------------------------
SHUF <- file.path(OUT, "joint_dishuf2.fa")
if (!file.exists(SHUF))
  system(sprintf("fasta-shuffle-letters -kmer 2 -copies 50 -seed 1 -dna %s > %s 2>/dev/null",
                 shQuote(JOINT_FA), shQuote(SHUF)))

subm <- c(full = PWM_FILE,
          left_ebox_core = file.path(OUT, "m2_left.meme"),
          right_atail    = file.path(OUT, "m2_right.meme"))
stopifnot(all(file.exists(subm)))

r3 <- do.call(rbind, lapply(names(subm), function(nm) {
  ob <- fimo_txt(subm[[nm]], JOINT_FA, paste0("obs_", nm))
  sh <- fimo_txt(subm[[nm]], SHUF,     paste0("shuf_", nm))
  n_obs  <- if (is.null(ob)) 0L else length(unique(ob$sequence_name))
  n_shuf <- if (is.null(sh)) 0L else length(unique(sh$sequence_name)) / 50
  data.frame(motif = nm,
             seqs_with_hit_obs = n_obs, pct_obs = round(100 * n_obs / 614, 1),
             seqs_with_hit_null = round(n_shuf, 1),
             pct_null = round(100 * n_shuf / 614, 1),
             fold = round(n_obs / max(n_shuf, 0.5), 2))
}))
resB$R3 <- r3
cat("\n== TEST R3  half vs whole (vs dinucleotide-shuffled null) ==\n")
print(r3, row.names = FALSE)

# ------------------------------------------------------------
# TEST R4, split-half generalisation.
# Discover de novo on a random half; score the held-out half with the
# discovered PWM; compare with the held-out half's shuffled control.
# A noise/overfit motif does not transfer.
# ------------------------------------------------------------
idx <- sample(seq_along(fa))
A <- fa[idx[1:307]]; B <- fa[idx[308:614]]
faA <- file.path(OUT, "halfA.fa"); faB <- file.path(OUT, "halfB.fa")
writeXStringSet(A, faA); writeXStringSet(B, faB)
system(sprintf("fasta-shuffle-letters -kmer 2 -copies 50 -seed 7 -dna %s > %s 2>/dev/null",
               shQuote(faB), shQuote(file.path(OUT, "halfB_shuf.fa"))))

system(sprintf("streme --p %s --oc %s --dna --nmotifs 5 --minw 6 --maxw 20 --thresh 0.05 >/dev/null 2>&1",
               shQuote(faA), shQuote(file.path(OUT, "streme_halfA"))))
sa <- readLines(file.path(OUT, "streme_halfA", "streme.txt"))
mo <- grep("^MOTIF ", sa, value = TRUE)
cat("\n== TEST R4  split-half generalisation ==\nSTREME on half A recovered:\n")
cat(paste0("  ", mo, collapse = "\n"), "\n")

# TOMTOM: does any half-A motif match the canonical Motif2 PWM?
system(sprintf("tomtom -no-ssc -oc %s -verbosity 1 -min-overlap 5 -dist pearson -thresh 0.5 %s %s >/dev/null 2>&1",
               shQuote(file.path(OUT, "tomtom_halfA_vs_m2")),
               shQuote(file.path(OUT, "streme_halfA", "streme.txt")), shQuote(PWM_FILE)))
tt <- read.delim(file.path(OUT, "tomtom_halfA_vs_m2", "tomtom.tsv"), stringsAsFactors = FALSE)
tt <- tt[!startsWith(as.character(tt[[1]]), "#") & nzchar(as.character(tt[[1]])), , drop = FALSE]
cat("TOMTOM half-A motifs vs canonical Motif2:\n")
if (nrow(tt)) print(tt[, c("Query_ID", "p.value", "q.value", "Query_consensus")], row.names = FALSE) else cat("  none\n")

# score held-out half B with the CANONICAL PWM (independent of half A)
obB <- fimo_txt(PWM_FILE, faB, "heldoutB")
shB <- fimo_txt(PWM_FILE, file.path(OUT, "halfB_shuf.fa"), "heldoutB_shuf")
nB  <- length(unique(obB$sequence_name)); nBs <- length(unique(shB$sequence_name)) / 50
resB$R4 <- data.frame(halfA_top_motif = sub("^MOTIF ", "", mo[1]),
                      tomtom_best_q = if (nrow(tt)) min(as.numeric(tt$q.value)) else NA_real_,
                      heldoutB_seqs_with_motif2 = nB,
                      heldoutB_pct = round(100 * nB / 307, 1),
                      heldoutB_null_seqs = round(nBs, 1),
                      heldoutB_fold = round(nB / max(nBs, 0.5), 2))
print(resB$R4, row.names = FALSE)

# ------------------------------------------------------------
# TEST R5, HOT / hyper-ChIPability.
# For each joint peak count how many of the other embryonic TF peak sets
# overlap it. Then test whether Motif2-positive peaks are the promiscuous ones.
# ------------------------------------------------------------
np_dir <- file.path(PROJECT_DIR, "data/final_peaks")
np <- c(SIM = "SIM_IDR0.05.clean.narrowPeak", DYS = "DYS_IDR0.05.clean.narrowPeak",
        SIMA = "SIMA_IDR0.05.clean.narrowPeak", TWI = "TWI_IDR0.05.final.narrowPeak",
        VVL = "VVL_IDR0.05.final.narrowPeak", RIB = "RIB_IDR0.05.final.narrowPeak")
np <- np[file.exists(file.path(np_dir, np))]
others <- lapply(np, function(f) {
  d <- read.delim(file.path(np_dir, f), header = FALSE, stringsAsFactors = FALSE)
  GRanges(d$V1, IRanges(d$V2 + 1L, d$V3))
})
n_other <- rowSums(vapply(others, function(g) as.integer(overlapsAny(peaks, g)),
                          integer(length(peaks))))
ob <- fimo_txt(PWM_FILE, JOINT_FA, "m2_forHOT")
m2pos <- peaks$name %in% unique(ob$sequence_name)
tt5 <- wilcox.test(n_other[m2pos], n_other[!m2pos])
hot <- n_other >= 4L
ft5 <- fisher.test(table(factor(m2pos, c(FALSE, TRUE)), factor(hot, c(FALSE, TRUE))))
resB$R5 <- data.frame(n_control_sets = length(np),
                      median_other_TFs_m2pos = median(n_other[m2pos]),
                      median_other_TFs_m2neg = median(n_other[!m2pos]),
                      mean_m2pos = round(mean(n_other[m2pos]), 2),
                      mean_m2neg = round(mean(n_other[!m2pos]), 2),
                      wilcox_p = signif(tt5$p.value, 3),
                      hot_OR = round(unname(ft5$estimate), 2),
                      hot_p = signif(ft5$p.value, 3))
cat("\n== TEST R5  HOT / hyper-ChIPability ==\n"); print(resB$R5, row.names = FALSE)
cat("distribution of co-bound TF count (of ", length(np), "):\n", sep = "")
print(table(n_other, Motif2 = m2pos))

# ------------------------------------------------------------
# TEST R6, positional concentration relative to peak centre
# ------------------------------------------------------------
w  <- width(fa); names(w) <- names(fa)
ob$mid <- (as.integer(ob$start) + as.integer(ob$stop)) / 2
ob$rel <- (ob$mid - w[ob$sequence_name] / 2) / (w[ob$sequence_name] / 2)  # -1..1
inner <- mean(abs(ob$rel) <= 0.5)
resB$R6 <- data.frame(n_instances = nrow(ob),
                      median_abs_rel_pos = round(median(abs(ob$rel)), 3),
                      pct_in_central_half = round(100 * inner, 1),
                      binom_p_vs_uniform = signif(
                        binom.test(sum(abs(ob$rel) <= 0.5), nrow(ob), 0.5)$p.value, 3))
cat("\n== TEST R6  positional concentration ==\n"); print(resB$R6, row.names = FALSE)

saveRDS(resB, file.path(OUT, "motif2_robustness_partB.rds"))
cat("\nPart B written to ", OUT, "\n", sep = "")