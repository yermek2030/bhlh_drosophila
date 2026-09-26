#!/usr/bin/env Rscript
# ============================================================
# icistarget_trh_tgo_summary.R
# Writes i-cisTarget novelty screen and Motif2-positive/negative
# co-enrichment values to results/reported_values.tsv, parsed from
# data/icistarget/trh_tgo*/icistarget/statistics.tbl.gz and data/motif2_gaga/motif2_vs_gaga_results.tsv.
# Usage: Rscript 02_motif_analysis/icistarget/icistarget_trh_tgo_summary.R
# ============================================================

## Repository root: the environment variable PROJECT_DIR, or the working directory.
PROJECT_DIR <- normalizePath(Sys.getenv("PROJECT_DIR", unset = "."), mustWork = TRUE)
if (!file.exists(file.path(PROJECT_DIR, "README.md")))
  stop("PROJECT_DIR does not point at the repository root: ", PROJECT_DIR)

DATA_DIR      <- file.path(PROJECT_DIR, "data")

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

read_value <- function(name) {
  tab <- read.delim(VALUES_FILE, colClasses = "character", quote = "", comment.char = "")
  tab$value[match(name, tab$name)]
}


AUX <- DATA_DIR
D6  <- file.path(AUX, "icistarget/trh_tgo")
DP  <- file.path(AUX, "icistarget/trh_tgo_motif2_positive")
DN  <- file.path(AUX, "icistarget/trh_tgo_motif2_negative")
DG  <- file.path(AUX, "motif2_gaga")

# ---- formatting helpers ------------------------------------------------------
# p/q to two significant figures in scientific notation; no thousands
# separators.
sci2 <- function(x) {
  if (is.na(x)) return("NA")
  if (x == 0) return("0")
  e <- floor(log10(abs(x))); m <- x / 10^e
  if (round(m, 1) == 10) { m <- 1; e <- e + 1 }
  sprintf("%.1fe%d", m, e)
}
dec2 <- function(x) sprintf("%.2f", x)          # odds ratios, NES
pct1 <- function(x) sprintf("%.1f", x)          # percentages
int0 <- function(x) formatC(as.integer(x), format = "d", big.mark = "")

read_stats <- function(d) {
  f <- file.path(d, "icistarget", "statistics.tbl")
  if (!file.exists(f)) f <- paste0(f, ".gz")
  stopifnot(file.exists(f))
  read.delim(f, quote = "", comment.char = "", stringsAsFactors = FALSE)
}
n_tomtom_hits <- function(d, sub) {
  f <- file.path(d, sub, "tomtom.tsv")
  if (!file.exists(f)) return(NA_integer_)
  x <- readLines(f, warn = FALSE)
  x <- x[nzchar(trimws(x)) & !startsWith(x, "#")]
  max(length(x) - 1L, 0L)                        # drop the header row
}
nearest <- function(d, sub) {
  f <- file.path(d, sub, "tomtom.tsv")
  t <- read.delim(f, comment.char = "#", stringsAsFactors = FALSE)
  t <- t[!is.na(t$p.value), ]
  t[order(t$E.value)[1], ]
}
# Search the same haystack as icistarget_trh_tgo_screen.R: FeatureID
# alone under-reports, because a matrix's factor identity often lives in
# the annotation column rather than the identifier (transfac_pro__M00723 is the
# top Trl feature and carries "Trl" only as an annotation; cisbp__M5928 is
# Grainyhead-family and says so only there). Keep the haystacks and the
# synonym unions identical across the two scripts so the values reported in
# results/reported_values.tsv agree with the tables they describe.
hay_of <- function(st) toupper(paste(st$FeatureID, st$FeatureDescription,
                                     st$FeatureAnnotations))
SYN <- list(Zld   = "ZLD|ZELDA",
            Trl   = "TRL|GAF|GAGA",
            Opa   = "OPA|ODD-PAIRED",
            Grh   = "GRH|GRAINY|GRHL|AACCGGTT|TFCP2",
            CLAMP = "CLAMP")
# Each synonym is anchored with \b, matching icistarget_motif2_analysis's own
# matching exactly. Without the anchor the two disagree: unanchored "GAGA" also
# matches inside longer tokens and lifted Sima/Tgo's Trl score from 3.65 to
# 4.10. Anchored is the stricter and the reported definition; keep both scripts
# on it.
idx_of <- function(st, pat) {
  syns <- strsplit(pat, "\\|", fixed = FALSE)[[1]]
  unique(unlist(lapply(syns, function(p) grep(paste0("\\b", p), hay_of(st)))))
}
best_nes <- function(st, pat) {
  i <- idx_of(st, pat)
  if (!length(i)) return(NA_real_)
  max(st$NES[i], na.rm = TRUE)
}
n_sig <- function(st, pat) { i <- idx_of(st, pat); sum(st$NES[i] >= 3, na.rm = TRUE) }

s6 <- read_stats(D6); sp <- read_stats(DP); sn <- read_stats(DN)
stopifnot(nrow(s6) == nrow(sp), nrow(sp) == nrow(sn))   # same feature space

npos <- length(unique(read.delim(
  file.path(DP, "icistarget", "input_mapped_to_icistarget_regions.bed"),
  header = FALSE, stringsAsFactors = FALSE)$V4))
nneg <- length(unique(read.delim(
  file.path(DN, "icistarget", "input_mapped_to_icistarget_regions.bed"),
  header = FALSE, stringsAsFactors = FALSE)$V4))
stopifnot(npos == 178L, nneg == 436L, npos + nneg == 614L)

np <- nearest(DP, "tt_icis_near")
nn <- nearest(DN, "tt_icis_near")

# ---- information content of the one significant match ----------------------
read_pfm <- function(f, nm) {
  l <- readLines(f, warn = FALSE); s <- grep(paste0("^MOTIF .*", nm), l)
  stopifnot(length(s) == 1L)
  w <- as.integer(sub(".*w=\\s*([0-9]+).*", "\\1", l[s + 1]))
  m <- do.call(rbind, lapply(strsplit(trimws(l[(s + 2):(s + 1 + w)]), "\\s+"), as.numeric))
  m / rowSums(m)
}
fm  <- read_pfm(file.path(DP, "icistarget_289.meme"), "fantom__motif80")
fic <- apply(fm, 1, function(p) { p <- p[p > 0]; 2 + sum(p * log2(p)) })

# ---- GAGA discrimination test ----------------------------------------------
g <- read.delim(file.path(DG, "motif2_vs_gaga_results.tsv"), stringsAsFactors = FALSE)
gv <- function(k) { v <- g$value[g$test == k]; stopifnot(length(v) == 1L); as.numeric(v) }

fantom_tt <- "fantom__motif80_CTCTGAGAGAA"

vals <- list(
  # ---- screen scale ---------------------------------------------------------
  nIcisPwms            = int0(sum(sp$FeatureDatabase == "PWMs")),
  nIcisFeatures        = int0(nrow(sp)),
  nIcisRegionsPos      = int0(npos),
  nIcisRegionsNeg      = int0(nneg),
  nIcisSigPos          = int0(sum(sp$NES >= 3, na.rm = TRUE)),
  nIcisSigNeg          = int0(sum(sn$NES >= 3, na.rm = TRUE)),
  icisNesThresh        = "3",
  # ---- novelty: TOMTOM against the enriched matrices ------------------------
  nIcisTomtomHitsAll   = int0(n_tomtom_hits(D6, "tt_icis_q05")),
  nIcisTomtomHitsPos   = int0(n_tomtom_hits(DP, "tt_icis_q05")),
  nIcisTomtomHitsNeg   = int0(n_tomtom_hits(DN, "tt_icis_q05")),
  nIcisCoreHitsPos     = int0(n_tomtom_hits(DP, "tt_core9_q05")),
  icisNearestPosName   = fantom_tt,
  icisNearestPosE      = sci2(np$E.value),
  icisNearestPosQ      = sci2(np$q.value),
  icisNearestPosOverlap= int0(np$Overlap),
  icisNearestNegName   = "stark__RSWGAGMRHRR",
  icisNearestNegQ      = sci2(nn$q.value),
  # ---- the match is degenerate: IC contrast --------------------------------
  icisFantomWidth      = int0(nrow(fm)),
  icisFantomICtotal    = dec2(sum(fic)),
  icisFantomICmax      = dec2(max(fic)),
  icisFantomPurinePct  = pct1(100 * mean(rowSums(fm[, c(1, 3)]))),
  icisFantomNesPos     = dec2(best_nes(sp, "fantom__motif80")),
  icisFantomNesNeg     = dec2(best_nes(sn, "fantom__motif80")),
  icisFantomNesAll     = dec2(best_nes(s6, "fantom__motif80")),
  # ---- co-enrichment contrast ----------------------------------------------
  # Per-family NES values for all four heterodimers are computed by
  # icistarget_nes_table.R; the GAGA discrimination values by
  # 04_chapter_analyses/chapter3/ch3_tab_motif2_gaga.R.
  nGrhMatricesPos = int0(n_sig(sp, SYN$Grh)),
  nGrhMatricesNeg = int0(n_sig(sn, SYN$Grh)),
  nGrhMatricesAll = int0(n_sig(s6, SYN$Grh))
)

save_values("icistarget_trh_tgo", vals)
for (nm in names(vals)) cat(sprintf("  %-24s %s\n", nm, read_value(nm)))
