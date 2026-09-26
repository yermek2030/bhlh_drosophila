#!/usr/bin/env Rscript
# ============================================================
# icistarget_all_heterodimers_summary.R
# 4 x 2 design: four heterodimers split into Motif2-positive/negative peaks, eight i-cisTarget screens.
# Values parsed from run outputs. Synonym matching anchored with \b, as in icistarget_trh_tgo_screen.R and icistarget_trh_tgo_summary.R.
# Usage: Rscript 02_motif_analysis/icistarget/icistarget_all_heterodimers_summary.R
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


AUX <- DATA_DIR
AH  <- file.path(AUX, "icistarget/all_heterodimers")

# Trh/Tgo was run before the all-heterodimer batch and lives in its own dirs.
runs <- list(
  Trh  = list(pos = file.path(AUX, "icistarget/trh_tgo_motif2_positive"),
              neg = file.path(AUX, "icistarget/trh_tgo_motif2_negative")),
  Sim  = list(pos = file.path(AH, "Sim_Tgo_m2pos"),  neg = file.path(AH, "Sim_Tgo_m2neg")),
  Dys  = list(pos = file.path(AH, "Dys_Tgo_m2pos"),  neg = file.path(AH, "Dys_Tgo_m2neg")),
  Sima = list(pos = file.path(AH, "Sima_Tgo_m2pos"), neg = file.path(AH, "Sima_Tgo_m2neg"))
)
NES_T <- 3

read_stats <- function(d) {
  f <- file.path(d, "icistarget", "statistics.tbl")
  if (!file.exists(f)) f <- paste0(f, ".gz")
  stopifnot(file.exists(f))
  read.delim(f, quote = "", comment.char = "", stringsAsFactors = FALSE)
}
n_regions <- function(d) length(unique(read.delim(
  file.path(d, "icistarget", "input_mapped_to_icistarget_regions.bed"),
  header = FALSE, stringsAsFactors = FALSE)$V4))
n_hits <- function(d) {
  f <- file.path(d, "tt_icis_q05", "tomtom.tsv")
  if (!file.exists(f)) return(NA_integer_)
  x <- readLines(f, warn = FALSE); x <- x[nzchar(trimws(x)) & !startsWith(x, "#")]
  max(length(x) - 1L, 0L)
}
hay <- function(st) toupper(paste(st$FeatureID, st$FeatureDescription, st$FeatureAnnotations))
ix  <- function(st, s) unique(unlist(lapply(s, function(p) grep(paste0("\\b", p), hay(st)))))
bn  <- function(st, s) { i <- ix(st, s); if (!length(i)) NA_real_ else max(st$NES[i], na.rm = TRUE) }
nm  <- function(st, s) { i <- ix(st, s); sum(st$NES[i] >= NES_T, na.rm = TRUE) }

SL <- list(Grh = c("GRH", "GRAINY", "GRHL", "AACCGGTT", "TFCP2"),
           Trl = c("TRL", "GAF", "GAGA"),
           Zld = c("ZLD", "ZELDA"))

dec2 <- function(x) sprintf("%.2f", x)
int0 <- function(x) formatC(as.integer(x), format = "d", big.mark = "")

vals <- list()
tally <- list(GrhPos = 0, GrhNeg = 0, TrlPos = 0, TrlNeg = 0, ZldPos = 0, ZldNeg = 0)
for (h in names(runs)) {
  sp <- read_stats(runs[[h]]$pos); sn <- read_stats(runs[[h]]$neg)
  stopifnot(nrow(sp) == nrow(sn))
  vals[[paste0("nMotifTwoPos", h)]] <- int0(n_regions(runs[[h]]$pos))
  vals[[paste0("nMotifTwoNeg", h)]] <- int0(n_regions(runs[[h]]$neg))
  vals[[paste0("nTomtomPos",   h)]] <- int0(n_hits(runs[[h]]$pos))
  vals[[paste0("nTomtomNeg",   h)]] <- int0(n_hits(runs[[h]]$neg))
  for (fac in names(SL)) {
    bp <- bn(sp, SL[[fac]]); bg <- bn(sn, SL[[fac]])
    vals[[paste0(tolower(fac), "Pos", h)]] <- dec2(bp)
    vals[[paste0(tolower(fac), "Neg", h)]] <- dec2(bg)
    tally[[paste0(fac, "Pos")]] <- tally[[paste0(fac, "Pos")]] + (bp >= NES_T)
    tally[[paste0(fac, "Neg")]] <- tally[[paste0(fac, "Neg")]] + (bg >= NES_T)
  }
  vals[[paste0("grhMatPos", h)]] <- int0(nm(sp, SL$Grh))
  vals[[paste0("grhMatNeg", h)]] <- int0(nm(sn, SL$Grh))
}
word <- c("zero", "one", "two", "three", "four")
vals$nHeterodimersScreened <- word[length(runs) + 1L]
vals$nGrhSigPosSets <- word[tally$GrhPos + 1L]
vals$nGrhSigNegSets <- word[tally$GrhNeg + 1L]
vals$nTrlSigPosSets <- word[tally$TrlPos + 1L]
vals$nTrlSigNegSets <- word[tally$TrlNeg + 1L]
vals$nZldSigSets    <- word[tally$ZldPos + tally$ZldNeg + 1L]
vals$icisNesThreshAll <- as.character(NES_T)

save_values("icistarget_all_heterodimers", vals)
cat(sprintf("Grh significant in %s of %s positive sets, %s of %s negative sets\n",
            vals$nGrhSigPosSets, vals$nHeterodimersScreened,
            vals$nGrhSigNegSets, vals$nHeterodimersScreened))
cat(sprintf("Trl significant in %s positive, %s negative -- peak-set property\n",
            vals$nTrlSigPosSets, vals$nTrlSigNegSets))

