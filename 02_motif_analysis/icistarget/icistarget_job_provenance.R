#!/usr/bin/env Rscript
# ============================================================
# icistarget_job_provenance.R
# Origin for the eight i-cisTarget screens. Writes summary values to
# results/reported_values.tsv and per-screen job table to results/tables/icistarget_jobs.tsv.
# Job identifier, region mapping and assembly are read from report URLs in the saved result HTML.
# Region mapping is not one-to-one: a peak straddling a boundary maps to two or more regions, so region counts exceed peak counts; both are reported. Read-only. Usage: Rscript 02_motif_analysis/icistarget/icistarget_job_provenance.R
# ============================================================

## Repository root: the environment variable PROJECT_DIR, or the working directory.
PROJECT_DIR <- normalizePath(Sys.getenv("PROJECT_DIR", unset = "."), mustWork = TRUE)
if (!file.exists(file.path(PROJECT_DIR, "README.md")))
  stop("PROJECT_DIR does not point at the repository root: ", PROJECT_DIR)

DATA_DIR      <- file.path(PROJECT_DIR, "data")
RESULTS_DIR   <- file.path(PROJECT_DIR, "results")

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

## ---------------------------------------------------------------- the screens
## Ordered as the text discusses them: Trh/Tgo first (the discovery set),
## then the three heterodimers Chapter 4 tests the PWM against.
screens <- data.frame(
  het     = c("Trh/Tgo", "Trh/Tgo", "Sim/Tgo", "Sim/Tgo",
              "Dys/Tgo", "Dys/Tgo", "Sima/Tgo", "Sima/Tgo"),
  stratum = rep(c("positive", "negative"), 4),
  dir     = c("icistarget/trh_tgo_motif2_positive/icistarget",
              "icistarget/trh_tgo_motif2_negative/icistarget",
              "icistarget/all_heterodimers/Sim_Tgo_m2pos/icistarget",
              "icistarget/all_heterodimers/Sim_Tgo_m2neg/icistarget",
              "icistarget/all_heterodimers/Dys_Tgo_m2pos/icistarget",
              "icistarget/all_heterodimers/Dys_Tgo_m2neg/icistarget",
              "icistarget/all_heterodimers/Sima_Tgo_m2pos/icistarget",
              "icistarget/all_heterodimers/Sima_Tgo_m2neg/icistarget"),
  stringsAsFactors = FALSE
)

## job id: recovered from the report URL the server writes into every result
## page.
job_id_of <- function(d) {
  fs <- list.files(d, pattern = "\\.html$", recursive = TRUE, full.names = TRUE)
  for (f in head(fs, 40)) {
    m <- regmatches(readLines(f, warn = FALSE, n = 400),
                    regexpr("reports/[a-f0-9]{24,}", readLines(f, warn = FALSE, n = 400)))
    m <- m[nzchar(m)]
    if (length(m)) return(sub("reports/", "", m[1]))
  }
  NA_character_
}

n_mapped <- function(d) {
  f <- file.path(d, "input_mapped_to_icistarget_regions.bed")
  if (!file.exists(f)) return(NA_integer_)
  length(readLines(f, warn = FALSE))
}
n_submitted <- function(d) {
  f <- file.path(d, "input_mapped_to_icistarget_regions.bed")
  if (!file.exists(f)) return(NA_integer_)
  ## column 4 is the submitted region's own name; distinct names = peaks sent
  length(unique(read.table(f, sep = "\t", quote = "", comment.char = "")[[4]]))
}

screens$path      <- file.path(AUX, screens$dir)
screens$submitted <- vapply(screens$path, n_submitted, integer(1))
screens$mapped    <- vapply(screens$path, n_mapped,    integer(1))
screens$job       <- vapply(screens$path, job_id_of,   character(1))

stopifnot(all(!is.na(screens$submitted)), all(!is.na(screens$mapped)),
          all(screens$mapped >= screens$submitted))

## assembly the server reports its regions in, read from the query record
asm <- NA_character_
qf <- file.path(screens$path[1], "query.id.txt")
if (file.exists(qf)) {
  s <- readLines(qf, warn = FALSE, n = 1)
  m <- regmatches(s, regexpr("dmel_r[0-9.]+", s))
  if (length(m)) asm <- sub("dmel_r", "", m[1])
}
stopifnot(!is.na(asm))

## ------------------------------------------------------------ summary values
vals <- list(
  nIcisScreens         = as.character(nrow(screens)),
  nIcisSubmittedTotal  = as.character(sum(screens$submitted)),
  nIcisMappedTotal     = as.character(sum(screens$mapped)),
  icisAssembly         = asm,
  icisServerHost       = "gbiomed.kuleuven.be"
)
save_values("icistarget_jobs", vals)

## -------------------------------------------------------------------- table
out <- file.path(RESULTS_DIR, "tables", "icistarget_jobs.tsv")
dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
job_table <- data.frame(heterodimer = screens$het, stratum = screens$stratum,
                        peaks = screens$submitted, regions = screens$mapped,
                        job_identifier = ifelse(is.na(screens$job), "not recorded",
                                                screens$job),
                        stringsAsFactors = FALSE)
write.table(job_table, out, sep = "\t", quote = FALSE, row.names = FALSE)

cat(sprintf("%-9s %-9s %6s %8s  %s\n", "het", "stratum", "peaks", "regions", "job"))
for (i in seq_len(nrow(screens)))
  cat(sprintf("%-9s %-9s %6d %8d  %s\n", screens$het[i], screens$stratum[i],
              screens$submitted[i], screens$mapped[i], substr(screens$job[i], 1, 12)))
cat(sprintf("\ntotals: %d peaks -> %d regions   assembly r%s\n",
            sum(screens$submitted), sum(screens$mapped), asm))
cat("written: ", out, "\n", sep = "")
