#!/usr/bin/env Rscript
# ============================================================
# ch4_vienna_tiles_summary.R
#
# Writes Vienna Tiles (Kvon et al. 2014) validation values into
# results/reported_values.tsv (section ch4_vienna_tiles), read from output CSVs.
# Run after prepare_vienna_tiles_dm6.R and perm_enhancers_vienna_tiles.R.
# Usage: Rscript --vanilla ch4_vienna_tiles_summary.R
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


REF_DIR <- file.path(PROJECT_DIR, "data/reference_gene_lists")
EXP4    <- file.path(PROJECT_DIR, "data/chapter4")

fish <- read.csv(file.path(EXP4, "vt_tile_conditional_fisher.csv"),
                 stringsAsFactors = FALSE)
perm <- read.csv(file.path(EXP4, "perm_enhancers_vt_summary.csv"),
                 stringsAsFactors = FALSE)

## ---- reference-set sizes, recounted from the BEDs ---------------------------
nline <- function(f) length(readLines(file.path(REF_DIR, f), warn = FALSE))
n_tested <- nline("vt_tiles_tested_dm6.bed")
n_active <- nline("vt_enhancers_active_dm6.bed")
n_late   <- nline("vt_enhancers_active_stg9_16_dm6.bed")

## ---- tiled fraction of the six major arms, and median tile width ------------
## Both values are computed here rather than hard-coded, since Chapter 2 quotes
## both figures.
suppressPackageStartupMessages({
  library(GenomicRanges); library(GenomeInfoDb)
  library(BSgenome.Dmelanogaster.UCSC.dm6)
})
good_chr <- c("chr2L", "chr2R", "chr3L", "chr3R", "chr4", "chrX")
tb <- read.table(file.path(REF_DIR, "vt_tiles_tested_dm6.bed"),
                 header = FALSE, stringsAsFactors = FALSE)
tested_gr <- GRanges(tb[[1]], IRanges(tb[[2]] + 1L, tb[[3]]))
tested_gr <- keepSeqlevels(tested_gr, good_chr, pruning.mode = "coarse")
arm_bp    <- sum(seqlengths(BSgenome.Dmelanogaster.UCSC.dm6)[good_chr])
tiled_bp  <- sum(width(reduce(tested_gr)))
stopifnot(length(tested_gr) == n_tested, tiled_bp > 0, tiled_bp < arm_bp)

## ---- Fisher rows: primary activity definition = "any" -----------------------
f <- fish[fish$Activity == "any", ]
stopifnot(nrow(f) == 12L)
HET <- c("Trh/Tgo" = "Trh", "Sim/Tgo" = "Sim",
         "Dys/Tgo" = "Dys", "Sima/Tgo" = "Sima")

get1 <- function(df, het, stratum, col) {
  v <- df[df$Heterodimer == het & df$Stratum == stratum, col]
  stopifnot(length(v) == 1L)
  v
}

vals <- list(
  vtTilesTested       = n_tested,
  vtActiveN           = n_active,
  vtActivePct         = sprintf("%.1f", 100 * n_active / n_tested),
  vtActiveLateN       = n_late,
  vtUnboundActivePct  = sprintf("%.1f", get1(f, "Trh/Tgo", "all", "Unbound_active_pct")),
  ## stage-restricted baseline, quoted in Ch2 as the sensitivity check
  vtUnboundActivePctLate = sprintf("%.1f", get1(fish[fish$Activity == "late", ],
                                                "Trh/Tgo", "all", "Unbound_active_pct")),
  vtCoveragePct       = sprintf("%.1f", 100 * tiled_bp / arm_bp),
  vtMedianWidth       = median(width(tested_gr))
)

## Per-heterodimer distal and promoter arms. Identifiers are letters only.
for (het in names(HET)) {
  tag <- HET[[het]]
  for (st in c("distal", "promoter")) {
    S <- if (st == "distal") "Dist" else "Prom"
    vals[[paste0("vt", tag, S, "OR")]]      <- sprintf("%.2f", get1(f, het, st, "OR"))
    vals[[paste0("vt", tag, S, "ORloCI")]]  <- sprintf("%.2f", get1(f, het, st, "CI_lo"))
    vals[[paste0("vt", tag, S, "ORhiCI")]]  <- sprintf("%.2f", get1(f, het, st, "CI_hi"))
    vals[[paste0("vt", tag, S, "BoundN")]]  <- get1(f, het, st, "Tiles_bound")
    vals[[paste0("vt", tag, S, "ActiveN")]] <- get1(f, het, st, "Bound_active")
    vals[[paste0("vt", tag, S, "ActivePct")]] <- sprintf("%.1f", get1(f, het, st, "Bound_active_pct"))
    ## p in scientific notation below 1e-3, otherwise to three decimals
    p <- get1(f, het, st, "p")
    vals[[paste0("vt", tag, S, "P")]] <- if (p < 1e-3) {
      ex <- floor(log10(p)); sprintf("%.1fe%d", p / 10^ex, ex)
    } else sprintf("%.3f", p)
    ## permutation Z on the VT active-any reference, same stratum
    z <- perm$Z_score[perm$Reference == "VT_active_any" &
                      perm$Heterodimer == het & perm$Stratum == st]
    stopifnot(length(z) == 1L)
    vals[[paste0("vt", tag, S, "Z")]] <- sprintf("%.2f", z)
  }
}

## ---- check: no missing values ------------------------------------------------
stopifnot(!any(vapply(vals, function(x) is.na(x) || !nzchar(as.character(x)), logical(1))))

save_values("ch4_vienna_tiles", vals)
invisible(lapply(names(vals), function(n)
  cat(sprintf("  %-22s %s\n", n, as.character(vals[[n]])))))
