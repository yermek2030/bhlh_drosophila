## ============================================================
## dapple_spot_class_qc.R
## Dapple spot-find quality classes from FlyChip vsn+limma output
## (column all:spotfindStatus). One token per array, Cy3-Cy5 pair,
## e.g. "A-S A-A A-R" = 3 arrays x 2 channels = 6 calls (Accept/Show/Reject).
## Writes pooled and per-stage percentages to results/reported_values.tsv and per-stage TSV table.
## ============================================================

## Repository root: the environment variable PROJECT_DIR, or the working directory.
PROJECT_DIR <- normalizePath(Sys.getenv("PROJECT_DIR", unset = "."), mustWork = TRUE)
if (!file.exists(file.path(PROJECT_DIR, "README.md")))
  stop("PROJECT_DIR does not point at the repository root: ", PROJECT_DIR)

DATA_DIR      <- file.path(PROJECT_DIR, "data")
RESULTS_DIR   <- file.path(PROJECT_DIR, "results")
CH5_DIR <- file.path(DATA_DIR, "chapter5")

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


STAGES   <- c(R11 = "Eleven", R12 = "Twelve", R13 = "Thirteen",
              R14 = "Fourteen", R15 = "Fifteen")
STAGENUM <- c(R11 = 11L, R12 = 12L, R13 = 13L, R14 = 14L, R15 = 15L)

## ---- per-stage tabulation ---------------------------------------------------
tab <- do.call(rbind, lapply(names(STAGES), function(s) {
  f <- file.path(CH5_DIR, paste0(s, ".vsn.limma.fdr.txt.gz"))
  stopifnot(file.exists(f))
  d <- read.delim(f, skip = 1L, check.names = FALSE, quote = "",
                  colClasses = "character", stringsAsFactors = FALSE)
  stopifnot("all:spotfindStatus" %in% names(d))

  toks <- strsplit(trimws(d[["all:spotfindStatus"]]), "[ -]+")
  ncal <- vapply(toks, length, 1L)
  stopifnot(length(unique(ncal)) == 1L)          # same n arrays for every probe
  cls  <- unlist(toks, use.names = FALSE)
  stopifnot(all(cls %in% c("A", "S", "R")))      # Dapple classes only

  data.frame(
    stage    = s,
    stagenum = STAGENUM[[s]],
    nProbes  = nrow(d),
    nArrays  = unique(ncal) / 2L,
    nCalls   = length(cls),
    nA = sum(cls == "A"), nS = sum(cls == "S"), nR = sum(cls == "R"),
    pctA = 100 * mean(cls == "A"),
    pctS = 100 * mean(cls == "S"),
    pctR = 100 * mean(cls == "R"),
    pctAllA = 100 * mean(vapply(toks, function(x) all(x == "A"), TRUE)),
    stringsAsFactors = FALSE)
}))
stopifnot(nrow(tab) == 5L, length(unique(tab$nProbes)) == 1L)

print(tab, row.names = FALSE, digits = 4)

## ---- pooled totals ----------------------------------------------------------
tot   <- colSums(tab[, c("nCalls", "nA", "nS", "nR")])
pctA  <- 100 * tot[["nA"]] / tot[["nCalls"]]
pctS  <- 100 * tot[["nS"]] / tot[["nCalls"]]
pctR  <- 100 * tot[["nR"]] / tot[["nCalls"]]
best  <- tab$stage[which.max(tab$pctA)]
worst <- tab$stage[which.min(tab$pctA)]

f1 <- function(x) formatC(x, format = "f", digits = 1)

vals <- list(
  spotClassNProbes      = format(unique(tab$nProbes), big.mark = ""),
  spotClassNArrays      = sum(tab$nArrays),
  spotClassNCalls       = format(tot[["nCalls"]], big.mark = ""),
  spotClassPctAccept    = f1(pctA),
  spotClassPctShow      = f1(pctS),
  spotClassPctReject    = f1(pctR),
  spotClassBestStage    = sub("^R", "", best),
  spotClassBestPct      = f1(max(tab$pctA)),
  spotClassWorstStage   = sub("^R", "", worst),
  spotClassWorstPct     = f1(min(tab$pctA)),
  limmaVersion          = "3.26.9",
  vsnVersion            = "3.38.0",
  biocVersion           = "3.2"
)
for (s in names(STAGES)) {
  r <- tab[tab$stage == s, ]
  vals[[paste0("spotClassPctAcceptSt", STAGES[[s]])]] <- f1(r$pctA)
  vals[[paste0("spotClassPctRejectSt", STAGES[[s]])]] <- f1(r$pctR)
}

save_values("ch2", vals)

## ---- per-stage table for the appendix ---------------------------------------
spot_tab <- data.frame(
  stage                 = c(paste("Stage", tab$stagenum), "All stages"),
  n_calls               = c(format(tab$nCalls, big.mark = ""),
                            format(tot[["nCalls"]], big.mark = "")),
  pct_accept            = f1(c(tab$pctA, pctA)),
  pct_show              = f1(c(tab$pctS, pctS)),
  pct_reject            = f1(c(tab$pctR, pctR)),
  pct_probes_all_accept = c(f1(tab$pctAllA), NA),
  stringsAsFactors = FALSE)
out <- file.path(RESULTS_DIR, "tables", "dapple_spot_class_by_stage.tsv")
dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
write.table(spot_tab, out, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
cat("wrote:", out, "\n")
