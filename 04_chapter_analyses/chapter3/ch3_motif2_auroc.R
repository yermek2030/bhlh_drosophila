#!/usr/bin/env Rscript
# ============================================================
# ch3_motif2_auroc.R
# Threshold-free AUROC decomposition of Motif2 (full matrix, CACTG core, A-rich tail); parses "R9 threshold-free AUROC" row from Motif2 robustness summary table.
# Writes results/reported_values.tsv (section ch3_motif2_auroc). Usage: Rscript 04_chapter_analyses/chapter3/ch3_motif2_auroc.R
## Repository root: PROJECT_DIR env var, or working directory.
# ============================================================
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

read_value <- function(name) {
  tab <- read.delim(VALUES_FILE, colClasses = "character", quote = "", comment.char = "")
  tab$value[match(name, tab$name)]
}

SUM <- file.path(PROJECT_DIR, "data/motif2_robustness/motif2_robustness_summary.csv")
stopifnot(file.exists(SUM))

## ---- 1. parse the AUROC row ----------------------------------------------
d <- read.csv(SUM, stringsAsFactors = FALSE)
row <- d[grepl("threshold-free AUROC", d$test, fixed = TRUE), , drop = FALSE]
stopifnot(nrow(row) == 1L)
res <- row$result[1]
grab <- function(key) {
  m <- regmatches(res, regexpr(sprintf("%s\\s+([0-9.]+)", key), res))
  stopifnot(length(m) == 1L)
  as.numeric(sub(sprintf("^%s\\s+", key), "", m))
}
vals <- list(
  motifTwoAurocFull = sprintf("%.3f", grab("full")),
  motifTwoAurocCore = sprintf("%.3f", grab("core")),
  motifTwoAurocTail = sprintf("%.3f", grab("tail")))

## ---- 2. stop if a previously recorded value differs ------------------------
for (m in names(vals)) {
  old <- read_value(m)
  if (!is.na(old) && abs(as.numeric(old) - as.numeric(vals[[m]])) > 1e-9)
    stop(sprintf("%s: summary table gives %s, reported_values.tsv has %s.",
                 m, vals[[m]], old))
}
cat("AUROC values parsed from the robustness summary:\n")
print(unlist(vals))

## ---- 3. record ---------------------------------------------------------------
save_values("ch3_motif2_auroc", vals)
