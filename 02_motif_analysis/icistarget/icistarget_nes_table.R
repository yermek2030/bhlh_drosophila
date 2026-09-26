#!/usr/bin/env Rscript
# ============================================================
# icistarget_nes_table.R: best NES per motif family (Grh, Trl, Zld, Odd-paired,
# Clamp) in Motif2-positive/negative peaks per heterodimer, from screen outputs,
# for four-way appendix table. Mismatch with results/reported_values.tsv stops script.
# Usage: Rscript 02_motif_analysis/icistarget/icistarget_nes_table.R
# ============================================================
## Repository root: PROJECT_DIR env var, or working directory.
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

## ---- the nine screens -----------------------------------------------------
PE <- "icistarget_pioneer_enrichment.tsv"
screens <- c(
  Trh_pos  = file.path(AUX, "icistarget/trh_tgo_motif2_positive", PE),
  Trh_neg  = file.path(AUX, "icistarget/trh_tgo_motif2_negative", PE),
  Sim_pos  = file.path(AUX, "icistarget/all_heterodimers/Sim_Tgo_m2pos",  PE),
  Sim_neg  = file.path(AUX, "icistarget/all_heterodimers/Sim_Tgo_m2neg",  PE),
  Dys_pos  = file.path(AUX, "icistarget/all_heterodimers/Dys_Tgo_m2pos",  PE),
  Dys_neg  = file.path(AUX, "icistarget/all_heterodimers/Dys_Tgo_m2neg",  PE),
  Sima_pos = file.path(AUX, "icistarget/all_heterodimers/Sima_Tgo_m2pos", PE),
  Sima_neg = file.path(AUX, "icistarget/all_heterodimers/Sima_Tgo_m2neg", PE),
  All      = file.path(AUX, "icistarget/trh_tgo", PE))
stopifnot(all(file.exists(screens)))

rd <- function(p) {
  d <- read.delim(p, stringsAsFactors = FALSE)
  stopifnot(all(c("factor","best_NES","significant") %in% names(d)))
  d
}
S <- lapply(screens, rd)

## factor label in the TSV -> value key used in results/reported_values.tsv
FAM <- c(Grh = "grh", Trl = "trl", Zld = "zld", Opa = "opa", CLAMP = "clamp")
stopifnot(all(names(FAM) %in% S[[1]]$factor))

nes <- function(sc, fac) {
  d <- S[[sc]]; v <- d$best_NES[d$factor == fac]
  stopifnot(length(v) == 1L); as.numeric(v)
}
sig <- function(sc, fac) {
  d <- S[[sc]]; v <- d$significant[d$factor == fac]
  stopifnot(length(v) == 1L); isTRUE(as.logical(v))
}
dec2 <- function(x) sprintf("%.2f", x)
WORD <- c("zero","one","two","three","four")

## ---- build the value set --------------------------------------------------
vals <- list()
for (fac in names(FAM)) {
  st <- FAM[[fac]]
  for (h in c("Trh","Sim","Dys","Sima")) {
    vals[[paste0(st, "Pos", h)]] <- dec2(nes(paste0(h, "_pos"), fac))
    vals[[paste0(st, "Neg", h)]] <- dec2(nes(paste0(h, "_neg"), fac))
  }
  vals[[paste0(st, "NesAll")]] <- dec2(nes("All", fac))
  np <- sum(vapply(c("Trh","Sim","Dys","Sima"), function(h) sig(paste0(h,"_pos"), fac), logical(1)))
  nn <- sum(vapply(c("Trh","Sim","Dys","Sima"), function(h) sig(paste0(h,"_neg"), fac), logical(1)))
  vals[[paste0("n", toupper(substring(st,1,1)), substring(st,2), "SigPos")]] <- WORD[np + 1L]
  vals[[paste0("n", toupper(substring(st,1,1)), substring(st,2), "SigNeg")]] <- WORD[nn + 1L]
}

## ---- check against previously recorded values ----------------------------
getv <- function(m) read_value(m)
cmp <- do.call(rbind, lapply(names(vals), function(m) {
  old <- getv(m)
  data.frame(name = m, old = old, new = vals[[m]],
             existed = !is.na(old), stringsAsFactors = FALSE)
}))
num <- suppressWarnings(as.numeric(cmp$old)) ; newnum <- suppressWarnings(as.numeric(cmp$new))
bad <- cmp$existed & !is.na(num) & !is.na(newnum) & abs(num - newnum) > 1e-9
if (any(bad)) {
  print(cmp[bad, ], row.names = FALSE)
  stop("Derived NES disagrees with the recorded value.")
}
badw <- cmp$existed & is.na(num) & cmp$old != cmp$new
if (any(badw)) { print(cmp[badw, ], row.names = FALSE); stop("Word-valued entry disagrees.") }

## ---- record -----------------------------------------------------------------
save_values("ch4_icistarget", vals)

cat(sprintf("\nvalues recorded: %d  (reproduced %d previously recorded, %d new)\n",
            nrow(cmp), sum(cmp$existed), sum(!cmp$existed)))
cat("\nfour-way NES matrix (rows = family, Pos/Neg per heterodimer):\n")
M <- t(sapply(names(FAM), function(fac) {
  st <- FAM[[fac]]
  setNames(as.numeric(unlist(vals[paste0(st, rep(c("Pos","Neg"), 4), rep(c("Trh","Sim","Dys","Sima"), each = 2))])),
           paste0(rep(c("Pos","Neg"), 4), rep(c("Trh","Sim","Dys","Sima"), each = 2)))
}))
print(M)
cat("\nsignificant-set tallies (of four heterodimers):\n")
print(unlist(vals[grep("SigPos$|SigNeg$", names(vals))]))
cat("\nnew entries:\n"); print(cmp$name[!cmp$existed])
