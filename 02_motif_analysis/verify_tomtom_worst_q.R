#!/usr/bin/env Rscript
## ============================================================
## verify_tomtom_worst_q.R
## Recomputes worst-of-twelve TOMTOM q-value for four-way Motif2
## cross-comparison from TOMTOM output; checks against results/reported_values.tsv.
## Run: Rscript 02_motif_analysis/verify_tomtom_worst_q.R
## Multiple four-way runs exist in data/; correct run identified by fingerprint, not path.
## ============================================================

## Repository root: the environment variable PROJECT_DIR, or the working directory.
PROJECT_DIR <- normalizePath(Sys.getenv("PROJECT_DIR", unset = "."), mustWork = TRUE)
if (!file.exists(file.path(PROJECT_DIR, "README.md")))
  stop("PROJECT_DIR does not point at the repository root: ", PROJECT_DIR)

RESULTS_DIR   <- file.path(PROJECT_DIR, "results")

VALUES_FILE <- file.path(PROJECT_DIR, "results", "reported_values.tsv")

read_value <- function(name) {
  tab <- read.delim(VALUES_FILE, colClasses = "character", quote = "", comment.char = "")
  tab$value[match(name, tab$name)]
}

## ---- 1. find every TOMTOM result, keep the four-way runs --------------------
tt <- list.files(PROJECT_DIR, pattern = "^tomtom.*\\.tsv$",
                 recursive = TRUE, full.names = TRUE)

read_tt <- function(f) {
  x <- try(read.delim(f, comment.char = "#", stringsAsFactors = FALSE),
           silent = TRUE)
  if (inherits(x, "try-error") || !nrow(x) || !"Query_ID" %in% names(x)) return(NULL)
  x$.file <- f
  x
}

runs <- Filter(Negate(is.null), lapply(tt, read_tt))

## a four-way heterodimer run: 4 queries, 4 targets, IDs naming the heterodimers
is_fourway <- function(x) {
  ids <- unique(c(x$Query_ID, x$Target_ID))
  length(unique(x$Query_ID)) == 4 && length(unique(x$Target_ID)) == 4 &&
    all(sapply(c("Trh", "Sim", "Dys", "Sima"),
               function(p) any(grepl(p, ids, ignore.case = TRUE))))
}
fourway <- Filter(is_fourway, runs)

stopifnot("no four-way TOMTOM run found" = length(fourway) > 0)

## ---- 2. summarise each four-way run ----------------------------------------
short <- function(id) sub("^Motif2_", "", sub("_Tgo$|_core$", "", id))

summarise_run <- function(x) {
  cross <- x[x$Query_ID != x$Target_ID, ]
  w <- cross[which.max(cross$q.value), ]
  b <- cross[which.min(cross$q.value), ]
  data.frame(
    file       = sub(paste0(PROJECT_DIR, "/"), "", x$.file[1], fixed = TRUE),
    n_cross    = nrow(cross),
    worst_q    = w$q.value,
    worst_pair = paste0(short(w$Query_ID), "->", short(w$Target_ID)),
    best_q     = b$q.value,
    best_pair  = paste0(short(b$Query_ID), "->", short(b$Target_ID)),
    min_overlap = min(cross$Overlap),
    max_overlap = max(cross$Overlap),
    stringsAsFactors = FALSE)
}
tab <- do.call(rbind, lapply(fourway, summarise_run))
tab <- tab[!duplicated(tab[, c("worst_q", "best_q")]), ]   # drop mirrored copies

cat("\n=== Four-way TOMTOM runs found (duplicated paths collapsed) ===\n")
print(tab[, c("n_cross", "worst_q", "worst_pair", "best_q", "best_pair",
              "min_overlap", "max_overlap")], row.names = FALSE, digits = 3)
cat("\npaths:\n"); cat(paste0("  ", tab$file, collapse = "\n"), "\n")

## ---- 3. identify the run the thesis describes, by fingerprint --------------
## Chapter 4 states: worst = Sima queried against Sim
##                   best  = Dys  queried against Trh
FP_WORST <- "Sima->Sim"
FP_BEST  <- "Dys->Trh"

hit <- tab[tab$worst_pair == FP_WORST & tab$best_pair == FP_BEST, ]

cat("\n=== Fingerprint match (Chapter 4 text) ===\n")
cat("expected worst pair:", FP_WORST, " expected best pair:", FP_BEST, "\n")
if (nrow(hit) == 1) {
  cat("MATCHED:", hit$file, "\n")
  cat(sprintf("  worst-of-%d q = %.6g   -> %.1e\n", hit$n_cross, hit$worst_q, hit$worst_q))
  cat(sprintf("  best       q = %.6g   -> %.1e\n", hit$best_q, hit$best_q))
} else {
  cat("NO UNIQUE MATCH (", nrow(hit), "candidates ) -- inspect manually.\n")
}

## ---- 4. check against the recorded values ---------------------------------
## chapter4_heterodimer_comparison.Rmd records tomtomWorstQ and tomtomTrhDysQ
## in results/reported_values.tsv and writes the pairwise table to
## results/tables/ch4_motif2_tomtom_pairwise.tsv.
parse_num <- function(s) suppressWarnings(as.numeric(gsub("\\s", "", s)))

checks <- data.frame(
  name     = c("tomtomWorstQ", "tomtomTrhDysQ"),
  observed = c(if (nrow(hit) == 1) hit$worst_q else NA_real_,
               if (nrow(hit) == 1) hit$best_q  else NA_real_),
  stringsAsFactors = FALSE)
checks$recorded       <- vapply(checks$name, read_value, character(1))
checks$recorded_value <- vapply(checks$recorded, parse_num, numeric(1))
## agree if the recorded value equals the observed value rounded to 2 significant figures
checks$agrees <- mapply(function(o, m)
  !is.na(o) && !is.na(m) && isTRUE(all.equal(signif(o, 2), m, tolerance = 1e-6)),
  checks$observed, checks$recorded_value)

cat("\n=== Agreement with results/reported_values.tsv ===\n")
print(checks[, c("name", "recorded", "observed", "agrees")],
      row.names = FALSE, digits = 3)

## order of magnitude below q = 0.05, stated as a word in the text
if (nrow(hit) == 1) {
  oom <- floor(log10(0.05 / hit$worst_q))
  cat(sprintf("\norders of magnitude below q = 0.05 : %.2f -> floor %d (recorded: '%s')\n",
              log10(0.05 / hit$worst_q), oom,
              read_value("tomtomWorstOrderOfMagnitude")))
}

## ---- 5. check every row of the Chapter 4 pairwise table -----------------------
## The table holds all 16 rows (12 cross + 4 self). Each row's q-value and
## overlap is compared against the matched TOMTOM run.
body_rows <- data.frame()
tab_path  <- file.path(RESULTS_DIR, "tables", "ch4_motif2_tomtom_pairwise.tsv")
if (nrow(hit) == 1 && file.exists(tab_path)) {
  tb <- read.delim(tab_path, colClasses = "character")
  nm <- function(s) sub("/Tgo$", "", gsub("\\s", "", s))
  body_rows <- data.frame(Query_ID = nm(tb$query), Target_ID = nm(tb$target),
                          q_tab = tb$q_value, overlap_tab = as.integer(tb$overlap),
                          stringsAsFactors = FALSE)
}

cat("\n=== Chapter 4 pairwise table vs TOMTOM output ===\n")
if (!nrow(body_rows)) {
  cat("results/tables/ch4_motif2_tomtom_pairwise.tsv not found -- SKIPPED\n")
  body_ok <- NA
} else {
  src <- read.delim(file.path(PROJECT_DIR, hit$file),
                    comment.char = "#", stringsAsFactors = FALSE)
  src$Query_ID  <- short(src$Query_ID)
  src$Target_ID <- short(src$Target_ID)
  cmp <- merge(body_rows, src[, c("Query_ID","Target_ID","q.value","Overlap")],
               by = c("Query_ID","Target_ID"), all.x = TRUE)
  cmp$q_recorded <- vapply(cmp$q_tab, parse_num, numeric(1))
  cmp$q_ok      <- mapply(function(a, b) !is.na(a) && !is.na(b) &&
                            isTRUE(all.equal(signif(b, 2), a, tolerance = 1e-6)),
                          cmp$q_recorded, cmp$q.value)
  cmp$ovl_ok    <- !is.na(cmp$Overlap) & cmp$overlap_tab == cmp$Overlap
  cat(sprintf("rows in table: %d (expect 16 = 12 cross + 4 self)\n", nrow(cmp)))
  cat(sprintf("q-values matching to 2 s.f.: %d/%d | overlaps matching: %d/%d\n",
              sum(cmp$q_ok), nrow(cmp), sum(cmp$ovl_ok), nrow(cmp)))
  bad <- cmp[!(cmp$q_ok & cmp$ovl_ok), ]
  if (nrow(bad)) {
    cat("MISMATCHED ROWS:\n")
    print(bad[, c("Query_ID","Target_ID","q_tab","q.value","overlap_tab","Overlap")],
          row.names = FALSE, digits = 3)
  }
  body_ok <- nrow(bad) == 0
}

## ---- 6. verdict ------------------------------------------------------------
cat("\n")
if (nrow(hit) != 1) {
  cat("RESULT: the described TOMTOM run was not identified.\n")
  quit(status = 1)
} else if (all(is.na(checks$recorded)) && is.na(body_ok)) {
  cat(sprintf("RESULT: worst-of-%d q = %.2g (%s); no recorded values to compare yet --\n",
              hit$n_cross, hit$worst_q, hit$worst_pair))
  cat("        run chapter4_heterodimer_comparison.Rmd to record them.\n")
} else if (all(checks$agrees) && isTRUE(body_ok)) {
  cat("RESULT: the recorded values agree with the TOMTOM output they describe,\n")
  cat("        including every row of the Chapter 4 pairwise table.\n")
} else if (all(checks$agrees) && is.na(body_ok)) {
  cat("RESULT: scalar values agree; pairwise table NOT checked (file missing).\n")
  quit(status = 2)
} else {
  cat("RESULT: MISMATCH between the recorded values and the TOMTOM output.\n")
  quit(status = 1)
}
