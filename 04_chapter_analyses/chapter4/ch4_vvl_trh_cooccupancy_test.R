# ============================================================
# ch4_vvl_trh_cooccupancy_test.R
# Tests whether Vvl's Motif2 carriage (36.8%) is inherited from Trh/Tgo
# co-occupancy. Splits existing 258 Vvl 201 bp summit windows (Motif2 FIMO
# scan already run) by Trh/Tgo IDR peak overlap: no peak 38.0%, peak 33.8%
# (or 0.83, p = 0.57). Run: Rscript ch4_vvl_trh_cooccupancy_test.R
# ============================================================

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

suppressPackageStartupMessages({
  library(GenomicRanges); library(GenomeInfoDb)
})


good_chr <- c("chr2L", "chr2R", "chr3L", "chr3R", "chr4", "chrX")
PK_DIR   <- file.path(PROJECT_DIR, "data/final_peaks")
CTRL_DIR <- file.path(PROJECT_DIR, "data/ctrl_tf", "fimo")
MOTIF_ID <- "STRCACTGARARAAA"

## ---- 1. the same peak loader used elsewhere in Chapter 4 -------------------
.load_peaks <- function(tag) {
  f <- file.path(PK_DIR, sprintf("%s_IDR0.05.final.narrowPeak", tag))
  stopifnot(file.exists(f))
  d  <- read.table(f, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
  gr <- GRanges(d[[1]], IRanges(d[[2]] + 1L, d[[3]]))     # BED -> 1-based
  gr <- gr[seqnames(gr) %in% good_chr]
  gr <- keepSeqlevels(gr, good_chr, pruning.mode = "coarse")
  gr[!duplicated(gr)]
}
trh <- .load_peaks("TRH")
tgo <- .load_peaks("TGO")

## ---- 2. the scanned Vvl windows, exactly as FIMO saw them ------------------
bed <- read.table(file.path(CTRL_DIR, "input_fastas", "VVL_enh_summits_200bp.bed"),
                  sep = "\t", header = FALSE, stringsAsFactors = FALSE)
win <- GRanges(bed[[1]], IRanges(bed[[2]] + 1L, bed[[3]]), name = bed[[4]])
stopifnot(length(win) == 258L, !any(duplicated(win$name)))

## ---- 3. Motif2 carriage, from the existing scan ---------------------------
fimo <- read.delim(file.path(CTRL_DIR, "VVL_motif2", "fimo.tsv"),
                   sep = "\t", header = TRUE, comment.char = "#",
                   stringsAsFactors = FALSE)
fimo <- fimo[!is.na(fimo$motif_id) & fimo$motif_id == MOTIF_ID, , drop = FALSE]
pos  <- unique(fimo$sequence_name)
stopifnot(length(pos) == 95L, all(pos %in% win$name))
win$motif2 <- win$name %in% pos

## ---- 4. the split: does a Trh / Tgo peak overlap the scanned window? -------
win$trh_ov <- overlapsAny(win, trh, minoverlap = 1L, ignore.strand = TRUE)
win$tgo_ov <- overlapsAny(win, tgo, minoverlap = 1L, ignore.strand = TRUE)
stopifnot(sum(win$trh_ov) == 74L)   # check: schema/ coordinate-frame drift

## ---- 5. contingency + Fisher ----------------------------------------------
tab2 <- function(keep_a, keep_b, lab_a, lab_b) {
  a <- win[keep_a]; b <- win[keep_b]
  m <- matrix(c(sum(a$motif2), sum(!a$motif2),
                sum(b$motif2), sum(!b$motif2)), nrow = 2,
              dimnames = list(c("Motif2+", "Motif2-"), c(lab_a, lab_b)))
  ft <- fisher.test(m)
  list(m = m, n_a = length(a), k_a = sum(a$motif2),
       pct_a = 100 * sum(a$motif2) / length(a),
       n_b = length(b), k_b = sum(b$motif2),
       pct_b = 100 * sum(b$motif2) / length(b),
       or = unname(ft$estimate), p = ft$p.value,
       ci = unname(ft$conf.int))
}

primary <- tab2( win$trh_ov, !win$trh_ov, "Trh peak", "no Trh peak")
strict  <- tab2( win$trh_ov & win$tgo_ov,
                !win$trh_ov & !win$tgo_ov, "Trh+Tgo", "neither")

cat("\n--- Vvl Motif2 carriage split by Trh co-occupancy ---\n")
print(primary$m)
cat(sprintf("Trh-overlapping  : %d/%d = %.1f%%\n", primary$k_a, primary$n_a, primary$pct_a))
cat(sprintf("no Trh peak      : %d/%d = %.1f%%\n", primary$k_b, primary$n_b, primary$pct_b))
cat(sprintf("OR = %.2f  95%% CI [%.2f, %.2f]  p = %.3g\n",
            primary$or, primary$ci[1], primary$ci[2], primary$p))
cat(sprintf("\nstrict (Trh AND Tgo vs neither): %.1f%% (n=%d) vs %.1f%% (n=%d), OR = %.2f, p = %.3g\n",
            strict$pct_a, strict$n_a, strict$pct_b, strict$n_b, strict$or, strict$p))

## ---- 6. reported values (no thousands separators) ---------------------------
fmt_p <- function(p) if (p >= 0.01) sprintf("%.2f", p) else {
  e <- floor(log10(p)); sprintf("%.2fe%d", p / 10^e, e)
}
vals <- list(
  nVvlWithTrhN    = format(primary$n_a, trim = TRUE),
  nVvlWithTrhPos  = format(primary$k_a, trim = TRUE),
  pctVvlWithTrh   = sprintf("%.1f", primary$pct_a),
  nVvlNoTrhN      = format(primary$n_b, trim = TRUE),
  nVvlNoTrhPos    = format(primary$k_b, trim = TRUE),
  pctVvlNoTrh     = sprintf("%.1f", primary$pct_b),
  vvlTrhSplitOR   = sprintf("%.2f", primary$or),
  vvlTrhSplitCIlo = sprintf("%.2f", primary$ci[1]),
  vvlTrhSplitCIhi = sprintf("%.2f", primary$ci[2]),
  vvlTrhSplitP    = fmt_p(primary$p)
)
save_values("ch4_vvl", vals)

## ---- 7. CSV of the split ---------------------------------------------------
out <- data.frame(
  Set     = c("Vvl windows overlapping a Trh peak", "Vvl windows with no Trh peak"),
  N       = c(primary$n_a, primary$n_b),
  Motif2  = c(primary$k_a, primary$k_b),
  Pct     = round(c(primary$pct_a, primary$pct_b), 1),
  OR      = c(round(primary$or, 2), NA),
  p_value = c(signif(primary$p, 3), NA),
  stringsAsFactors = FALSE
)
dir.create(file.path(PROJECT_DIR, "data", "chapter4"),
           showWarnings = FALSE, recursive = TRUE)
write.csv(out, file.path(PROJECT_DIR, "data", "chapter4",
                         "table4_vvl_trh_cooccupancy.csv"), row.names = FALSE)
cat("csv:", file.path("data", "chapter4",
                      "table4_vvl_trh_cooccupancy.csv"), "\n")
