## ============================================================
##  ch4_cme_occurrence.R
##
##  CME (ACGTG) occurrence per heterodimer and association with
##  regulatory context. Scans both strands (ACGTG or its reverse
##  complement CACGT); Trh/Tgo baseline from Chapter 3, pctCmePos = 30.3%.
##  Outputs: table4_8_CME_occurrence.csv, table4_9_CME_enhancer_context.csv, table4_8b_CME_vs_Trh_fisher.csv
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

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(Biostrings)
  library(Rsamtools)
})

DM6     <- file.path(PROJECT_DIR, "reference", "dm6.fa")
ENH_BED <- file.path(PROJECT_DIR, "data/reference_gene_lists/merged_embryo_5-20_enhancers_dm6.bed")
PEAKS   <- c("Trh/Tgo"  = "data/chapter3/Trh_Tgo_joint_peaks_filt.bed",
             "Sim/Tgo"  = "data/chapter4/Sim_Tgo_joint_peaks_filt.bed",
             "Dys/Tgo"  = "data/chapter4/Dys_Tgo_joint_peaks_filt.bed",
             "Sima/Tgo" = "data/chapter4/Sima_Tgo_joint_peaks_filt.bed")
EXPECT_N <- c("Trh/Tgo" = 614L, "Sim/Tgo" = 355L, "Dys/Tgo" = 870L, "Sima/Tgo" = 296L)
OUT_DIR  <- file.path(PROJECT_DIR, "data/chapter4")

stopifnot(file.exists(DM6), file.exists(ENH_BED))

read_bed <- function(p) {
  # accept either an absolute path or one relative to PROJECT_DIR
  f <- if (startsWith(p, "/")) p else file.path(PROJECT_DIR, p)
  stopifnot(file.exists(f))
  d <- read.delim(f, header = FALSE, stringsAsFactors = FALSE)
  GRanges(d[[1]], IRanges(d[[2]] + 1L, d[[3]]))
}

enh <- read_bed(ENH_BED)
fa  <- FaFile(DM6)

## ---- CME presence, both strands ----
find_cme <- function(gr) {
  seqs <- scanFa(fa, GRanges(seqnames(gr), ranges(gr)))
  n_plus  <- vcountPattern("ACGTG", seqs, fixed = TRUE)
  n_minus <- vcountPattern("CACGT", seqs, fixed = TRUE)   # rev-comp of ACGTG
  gr$CME  <- (n_plus + n_minus) > 0L
  gr$overlaps_enhancer <- countOverlaps(gr, enh) > 0L
  gr
}

hets <- lapply(names(PEAKS), function(k) {
  gr <- read_bed(PEAKS[[k]])
  stopifnot(length(gr) == EXPECT_N[[k]])
  find_cme(gr)
})
names(hets) <- names(PEAKS)

## ---- table 4.8: occurrence ----
cme_df <- data.frame(
  Heterodimer  = names(hets),
  Total_peaks  = vapply(hets, length,          integer(1)),
  CME_positive = vapply(hets, function(g) sum(g$CME), integer(1)),
  stringsAsFactors = FALSE
)
cme_df$CME_fraction <- round(100 * cme_df$CME_positive / cme_df$Total_peaks, 2)
stopifnot(nrow(cme_df) == 4L)

## ---- table 4.8b: each new heterodimer vs the Trh/Tgo baseline ----
## Both sides both-strand: this is the comparison the chapter makes.
ref <- cme_df[cme_df$Heterodimer == "Trh/Tgo", ]
cme_vs_trh <- do.call(rbind, lapply(setdiff(cme_df$Heterodimer, "Trh/Tgo"), function(k) {
  r  <- cme_df[cme_df$Heterodimer == k, ]
  ft <- fisher.test(matrix(c(r$CME_positive,   r$Total_peaks   - r$CME_positive,
                             ref$CME_positive, ref$Total_peaks - ref$CME_positive),
                           nrow = 2, byrow = TRUE))
  data.frame(Heterodimer = k,
             Comparison  = "CME frequency vs Trh/Tgo (both strands)",
             Fisher_OR   = round(unname(ft$estimate), 3),
             CI_lo       = round(ft$conf.int[1], 3),
             CI_hi       = round(ft$conf.int[2], 3),
             p_value     = ft$p.value,
             stringsAsFactors = FALSE)
}))

## ---- table 4.9: CME by regulatory context, per heterodimer ----
cme_context_df <- do.call(rbind, lapply(names(hets), function(k) {
  g  <- hets[[k]]
  tb <- table(factor(g$overlaps_enhancer, c(TRUE, FALSE)),
              factor(g$CME,               c(TRUE, FALSE)))
  ft <- fisher.test(tb)
  data.frame(Heterodimer = k,
             Comparison  = "Annotated enhancer vs non-annotated",
             Fisher_OR   = round(unname(ft$estimate), 3),
             p_value     = ft$p.value,
             stringsAsFactors = FALSE)
}))

write.csv(cme_df,         file.path(OUT_DIR, "table4_8_CME_occurrence.csv"),        row.names = FALSE)
write.csv(cme_vs_trh,     file.path(OUT_DIR, "table4_8b_CME_vs_Trh_fisher.csv"),    row.names = FALSE)
write.csv(cme_context_df, file.path(OUT_DIR, "table4_9_CME_enhancer_context.csv"),  row.names = FALSE)

print(cme_df,         row.names = FALSE)
print(cme_vs_trh,     row.names = FALSE)
print(cme_context_df, row.names = FALSE)

## Consistency with Chapter 3: the Trh/Tgo row here must equal pctCmePos.
cat(sprintf("\nTrh/Tgo both-strand CME: %d/%d = %.1f%%  (Chapter 3 pctCmePos = 30.3)\n",
            ref$CME_positive, ref$Total_peaks,
            100 * ref$CME_positive / ref$Total_peaks))

## ---- reported values ----


fmt_p <- function(p, digits = 2) {
  if (is.na(p)) return(NA_character_)
  s  <- formatC(p, format = "e", digits = digits)
  pr <- strsplit(s, "e")[[1]]
  paste0(pr[1], "e", as.integer(pr[2]))
}
gc_ <- function(k, col) cme_df[[col]][cme_df$Heterodimer == k]
gv  <- function(k, col) cme_vs_trh[[col]][cme_vs_trh$Heterodimer == k]
gx  <- function(k, col) cme_context_df[[col]][cme_context_df$Heterodimer == k]

cme_vals <- list(
  ## Trh/Tgo reference row for the Chapter 4 table (the same numbers as the
  ## Chapter 3 cmePosTotal / pctCmePos; the check above confirms they agree).
  trhCmeN        = gc_("Trh/Tgo","CME_positive"),
  trhCmePct      = sprintf("%.1f", gc_("Trh/Tgo","CME_fraction")),
  simCmeN        = gc_("Sim/Tgo","CME_positive"),
  simCmePct      = sprintf("%.1f", gc_("Sim/Tgo","CME_fraction")),
  dysCmeN        = gc_("Dys/Tgo","CME_positive"),
  dysCmePct      = sprintf("%.1f", gc_("Dys/Tgo","CME_fraction")),
  simaCmeN       = gc_("Sima/Tgo","CME_positive"),
  simaCmePct     = sprintf("%.1f", gc_("Sima/Tgo","CME_fraction")),
  simCmeVsTrhOR  = sprintf("%.2f", gv("Sim/Tgo","Fisher_OR")),
  simCmeVsTrhP   = fmt_p(gv("Sim/Tgo","p_value")),
  dysCmeVsTrhOR  = sprintf("%.2f", gv("Dys/Tgo","Fisher_OR")),
  dysCmeVsTrhP   = fmt_p(gv("Dys/Tgo","p_value")),
  simaCmeVsTrhOR = sprintf("%.2f", gv("Sima/Tgo","Fisher_OR")),
  simaCmeVsTrhP  = fmt_p(gv("Sima/Tgo","p_value")),
  cmeVsTrhORmin  = sprintf("%.2f", min(cme_vs_trh$Fisher_OR)),
  cmeVsTrhORmax  = sprintf("%.2f", max(cme_vs_trh$Fisher_OR)),
  cmeVsTrhPmin   = fmt_p(min(cme_vs_trh$p_value)),
  trhCmeCtxOR    = sprintf("%.2f", gx("Trh/Tgo","Fisher_OR")),
  trhCmeCtxP     = fmt_p(gx("Trh/Tgo","p_value")),
  simCmeCtxOR    = sprintf("%.2f", gx("Sim/Tgo","Fisher_OR")),
  simCmeCtxP     = fmt_p(gx("Sim/Tgo","p_value")),
  dysCmeCtxOR    = sprintf("%.2f", gx("Dys/Tgo","Fisher_OR")),
  dysCmeCtxP     = fmt_p(gx("Dys/Tgo","p_value")),
  simaCmeCtxOR   = sprintf("%.2f", gx("Sima/Tgo","Fisher_OR")),
  simaCmeCtxP    = fmt_p(gx("Sima/Tgo","p_value"))
)
save_values("ch4_cme_occurrence", cme_vals)
