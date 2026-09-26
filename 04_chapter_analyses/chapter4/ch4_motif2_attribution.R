## ============================================================
## ch4_motif2_attribution.R: tests whether Motif2 marks the four bHLH-PAS-Tgo heterodimers or generic distal enhancers
## Design: one enhancer catalogue split by heterodimer peak overlap, both halves scanned identically
## Conventions: 201 bp windows for all sets; p < 1e-4 not q (q keeps 19% of hits at 35797 enhancers vs 88% at 614 peaks)
## Outputs to data/motif2_attribution/: pooled_attribution.csv (per-stratum carriage), pooled_attribution_tests.csv (Fisher tests + 95% CI)
## ============================================================

## The figure fig4_12_motif2_attribution.png is drawn by
## chapter3_4_figures.Rmd from the tables written here.

## Repository root: the environment variable PROJECT_DIR, or the working directory.
PROJECT_DIR <- normalizePath(Sys.getenv("PROJECT_DIR", unset = "."), mustWork = TRUE)
if (!file.exists(file.path(PROJECT_DIR, "README.md")))
  stop("PROJECT_DIR does not point at the repository root: ", PROJECT_DIR)

suppressPackageStartupMessages({
  library(GenomicRanges)
})

## MEME Suite binaries. Resolution order: $MEME_BIN, then PATH.
## Set MEME_BIN (shell or ~/.Renviron) if the MEME Suite is not on PATH.
MEME_BIN <- Sys.getenv("MEME_BIN", unset = "")
if (!nzchar(MEME_BIN)) MEME_BIN <- dirname(Sys.which("fimo")[[1]])
if (!nzchar(MEME_BIN) || !file.exists(file.path(MEME_BIN, "fimo")))
  stop("MEME Suite not found. Install MEME >= 5.5.5 and set MEME_BIN to its bin/ directory.")
PWM      <- file.path(PROJECT_DIR, "data/chapter3/STRCACTGARARAAA.meme")
ENH_FA   <- file.path(PROJECT_DIR, "data/chapter3/enhancers_genomewide.fa")
OUT_DIR  <- file.path(PROJECT_DIR, "data/motif2_attribution")
WIN      <- 100L      # half-window: 201 bp total
PTHRESH  <- 1e-4

stopifnot(file.exists(PWM), file.exists(ENH_FA))
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

## ---- 1. helpers ----

read_fasta <- function(path) {
  L <- readLines(path, warn = FALSE)
  h <- grep("^>", L)
  nm <- sub("^>", "", L[h])
  starts <- h + 1L
  ends   <- c(h[-1] - 1L, length(L))
  seqs <- vapply(seq_along(h), function(i)
    paste0(L[starts[i]:ends[i]], collapse = ""), character(1))
  setNames(seqs, nm)
}

## centre 201 bp window; FIMO parses "chr:start-end" headers into
## genomic coordinates, so headers are replaced with opaque ids
write_windows <- function(seqs, path, prefix) {
  keep <- nchar(seqs) >= 150L
  seqs <- seqs[keep]
  mid  <- nchar(seqs) %/% 2L
  win  <- substr(seqs, pmax(1L, mid - WIN), mid + WIN)
  ids  <- paste0(prefix, seq_along(win))
  writeLines(as.vector(rbind(paste0(">", ids), win)), path)
  list(ids = ids, kept = keep)
}

## sequences carrying >= 1 hit at p < PTHRESH (q unused)
scan_pwm <- function(fa) {
  out <- file.path(tempdir(), paste0("fimo_", basename(fa)))
  unlink(out, recursive = TRUE)
  system2(file.path(MEME_BIN, "fimo"),
          c("--thresh", format(PTHRESH, scientific = TRUE),
            "--oc", shQuote(out), shQuote(PWM), shQuote(fa)),
          stdout = FALSE, stderr = FALSE)
  tsv <- file.path(out, "fimo.tsv")
  stopifnot(file.exists(tsv))
  d <- read.delim(tsv, comment.char = "#", stringsAsFactors = FALSE)
  d <- d[!is.na(d$p.value) & d$p.value < PTHRESH & nzchar(d$sequence_name), ]
  unique(d$sequence_name)
}

fisher_row <- function(label, a, na, b, nb, note) {
  ft <- fisher.test(matrix(c(a, na - a, b, nb - b), nrow = 2))
  se <- sqrt(1/max(a,1) + 1/max(na-a,1) + 1/max(b,1) + 1/max(nb-b,1))
  or <- as.numeric(ft$estimate)
  data.frame(comparison = label,
             a = a, n_a = na, pct_a = round(100*a/na, 1),
             b = b, n_b = nb, pct_b = round(100*b/nb, 1),
             OR = round(or, 2),
             CI_lo = round(exp(log(or) - 1.96*se), 2),
             CI_hi = round(exp(log(or) + 1.96*se), 2),
             p = signif(ft$p.value, 3),
             note = note, stringsAsFactors = FALSE)
}

## ---- 2. enhancer catalogue, length-normalised ----

enh_seqs <- read_fasta(ENH_FA)
enh_coord <- do.call(rbind, lapply(names(enh_seqs), function(nm) {
  m <- regmatches(nm, regexec("^([^:]+):(\\d+)-(\\d+)", nm))[[1]]
  if (length(m) == 4L) data.frame(chr = m[2], start = as.integer(m[3]),
                                  end = as.integer(m[4]))
  else data.frame(chr = NA_character_, start = NA_integer_, end = NA_integer_)
}))

enh_fa_win <- file.path(OUT_DIR, "enh_201bp.fa")
w <- write_windows(enh_seqs, enh_fa_win, "enh")
enh_coord <- enh_coord[w$kept, , drop = FALSE]
enh_gr <- GRanges(enh_coord$chr, IRanges(enh_coord$start, enh_coord$end))
cat(sprintf("enhancers scanned: %d\n", length(enh_gr)))

enh_hit_ids <- scan_pwm(enh_fa_win)
enh_hit <- as.integer(sub("^enh", "", enh_hit_ids))
cat(sprintf("  carrying Motif2: %d (%.1f%%)\n",
            length(enh_hit), 100*length(enh_hit)/length(enh_gr)))

## ---- 3. pool the four heterodimer peak sets ----

peak_beds <- c(
  Trh  = file.path(PROJECT_DIR, "data/chapter3/Trh_Tgo_joint_peaks_filt.bed"),
  Sim  = file.path(PROJECT_DIR, "data/chapter4/Sim_Tgo_joint_peaks_filt.bed"),
  Dys  = file.path(PROJECT_DIR, "data/chapter4/Dys_Tgo_joint_peaks_filt.bed"),
  Sima = file.path(PROJECT_DIR, "data/chapter4/Sima_Tgo_joint_peaks_filt.bed"))
stopifnot(all(file.exists(peak_beds)))

read_bed <- function(p) {
  d <- read.delim(p, header = FALSE, stringsAsFactors = FALSE)
  GRanges(d$V1, IRanges(d$V2 + 1L, d$V3))
}
peak_gr <- lapply(peak_beds, read_bed)
for (k in names(peak_gr)) cat(sprintf("  %-5s %d peaks\n", k, length(peak_gr[[k]])))

pooled_gr <- reduce(do.call(c, unname(peak_gr)))
cat(sprintf("pooled (merged) heterodimer peaks: %d\n", length(pooled_gr)))

## ---- 4. split the enhancer catalogue by binding ----

hit_flag <- logical(length(enh_gr)); hit_flag[enh_hit] <- TRUE

strata <- list()
for (k in names(peak_gr)) {
  b <- countOverlaps(enh_gr, peak_gr[[k]], ignore.strand = TRUE) > 0
  strata[[k]] <- list(bound = b, n = sum(b), hits = sum(b & hit_flag))
}
b_pool <- countOverlaps(enh_gr, pooled_gr, ignore.strand = TRUE) > 0
strata[["pooled"]] <- list(bound = b_pool, n = sum(b_pool),
                           hits = sum(b_pool & hit_flag))

per_het <- do.call(rbind, lapply(names(strata), function(k)
  data.frame(set = k, bound_enhancers = strata[[k]]$n,
             with_motif2 = strata[[k]]$hits,
             pct = round(100*strata[[k]]$hits/max(strata[[k]]$n,1), 1),
             stringsAsFactors = FALSE)))
print(per_het)

## ---- 5. tests: bound vs unbound, within one catalogue ----

tests <- do.call(rbind, lapply(names(strata), function(k) {
  s <- strata[[k]]
  ub_n <- length(enh_gr) - s$n
  ub_h <- sum(!s$bound & hit_flag)
  fisher_row(sprintf("enhancers bound by %s vs unbound", k),
             s$hits, s$n, ub_h, ub_n,
             "one catalogue, 201bp, p<1e-4, split by binding only")
}))
print(tests[, c("comparison","pct_a","pct_b","OR","CI_lo","CI_hi","p")])

## ---- 6. attribution split on the pooled stratum ----

bg_fa <- file.path(PROJECT_DIR, "data/ctrl_tf/fimo/input_fastas/matched_bg_200bp.fa")
if (file.exists(bg_fa)) {
  bg_n <- length(grep("^>", readLines(bg_fa, warn = FALSE)))
  bg_h <- length(scan_pwm(bg_fa))
  lg <- function(a, n) log((a/n)/(1 - a/n))
  s <- strata[["pooled"]]
  ub_n <- length(enh_gr) - s$n; ub_h <- sum(!s$bound & hit_flag)
  d_enh <- lg(ub_h, ub_n) - lg(bg_h, bg_n)
  d_fac <- lg(s$hits, s$n) - lg(ub_h, ub_n)
  cat(sprintf(
    "\nattribution (log-odds from matched background):\n  enhancer-ness  %.2f (%.0f%%, OR %.1f)\n  factor binding %.2f (%.0f%%, OR %.1f)\n",
    d_enh, 100*d_enh/(d_enh+d_fac), exp(d_enh),
    d_fac, 100*d_fac/(d_enh+d_fac), exp(d_fac)))
  attribution <- data.frame(
    component = c("enhancer-ness", "heterodimer binding"),
    delta_logit = round(c(d_enh, d_fac), 3),
    pct_of_total = round(100*c(d_enh, d_fac)/(d_enh + d_fac)),
    OR = round(exp(c(d_enh, d_fac)), 2))
  write.csv(attribution, file.path(OUT_DIR, "pooled_attribution_split.csv"),
            row.names = FALSE)
}

write.csv(per_het, file.path(OUT_DIR, "pooled_attribution.csv"), row.names = FALSE)
write.csv(tests,   file.path(OUT_DIR, "pooled_attribution_tests.csv"), row.names = FALSE)
cat(sprintf("\nwritten to %s\n", OUT_DIR))
