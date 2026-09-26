#!/usr/bin/env Rscript
# ============================================================
# motif2_grh_direct_scan.R
# Direct scan of the 11 Grainyhead-family matrices and confound values for
# the i-cisTarget Grainyhead co-enrichment result: direct scan shows no
# enrichment, Motif2+/- peaks differ 2.4-fold in promoter-proximity, and
# heterodimer Motif2+ sets are ~48% redundant. Reported as non-replication.
# Usage: Rscript 02_motif_analysis/cooccurrence/motif2_grh_direct_scan.R
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


suppressPackageStartupMessages(library(GenomicRanges))
AUX <- DATA_DIR
CF  <- file.path(AUX, "motif2_grh_confounds")
AH  <- file.path(AUX, "icistarget/all_heterodimers")

# ---- confound table written by motif2_grh_confounds.R -------------
cf <- read.delim(file.path(CF, "motif2_grh_confound_results.tsv"),
                 stringsAsFactors = FALSE)
gv <- function(k) { i <- match(k, cf$quantity); stopifnot(!is.na(i)); cf$value[i] }

# ---- direct FIMO contrast, recomputed here so the thresholds are explicit --
jp <- readRDS(file.path(AUX, "chapter3", "joint_peaks_filt.rds"))
stopifnot(length(jp) == 614)
rd <- function(f) { x <- read.table(f, header = TRUE, sep = "\t",
                                    comment.char = "#", stringsAsFactors = FALSE)
                    x[!is.na(x$start), ] }
m2 <- rd(file.path(AUX, "chapter3", "fimo_joint", "fimo.tsv"))
m2 <- m2[m2$motif_id == "STRCACTGARARAAA" & m2$q.value < 0.05, ]
M2 <- seq_along(jp) %in% unique(as.integer(sub("peak_", "", m2$sequence_name)))
gh <- rd(file.path(CF, "fimo_grh", "fimo.tsv"))
gh$peak <- as.integer(sub("peak_", "", gh$sequence_name))
n_grh_matrices <- length(unique(gh$motif_id))

fisher_at <- function(th) {
  g <- seq_along(jp) %in% unique(gh$peak[gh$p.value < th])
  f <- fisher.test(table(M2, g))
  list(pos = 100 * mean(g[M2]), neg = 100 * mean(g[!M2]),
       or = unname(f$estimate), p = f$p.value)
}
a3 <- fisher_at(1e-3); a4 <- fisher_at(1e-4)
nb <- tabulate(gh$peak[gh$p.value < 1e-3], nbins = length(jp)); w <- width(jp)
dens_ratio <- (sum(nb[M2]) / sum(w[M2])) / (sum(nb[!M2]) / sum(w[!M2]))

# ---- peak-set redundancy across the four heterodimers ----------------------
rdbed <- function(f) { b <- read.table(f, sep = "\t", stringsAsFactors = FALSE)
                       GRanges(b$V1, IRanges(b$V2 + 1L, b$V3)) }
P <- lapply(c("Trh_Tgo", "Sim_Tgo", "Dys_Tgo", "Sima_Tgo"), function(h)
  rdbed(file.path(AH, sprintf("%s_motif2_positive_icistarget.bed", h))))
n_sum <- sum(lengths(P)); n_uni <- length(reduce(do.call(c, P)))
pairmax <- max(vapply(seq_along(P), function(i) max(vapply(setdiff(seq_along(P), i),
  function(j) 100 * sum(countOverlaps(P[[i]], P[[j]]) > 0) / length(P[[i]]), 0)), 0))

dec2 <- function(x) sprintf("%.2f", x)
pc1  <- function(x) sprintf("%.1f", x)
sci  <- function(p) { e <- floor(log10(p)); m <- p / 10^e
                      sprintf("%.1fe%d", m, e) }

# The Motif2 instance count used by the GAGA discrimination text. Parsed from
# that run's own output rather than recomputed, so the denominator quoted in
# Chapter 3 is the same number the GAGA test divided by.
gg <- read.delim(file.path(AUX, "motif2_gaga", "motif2_vs_gaga_results.tsv"),
                 stringsAsFactors = FALSE)
gval <- function(k) { i <- match(k, gg$test); stopifnot(!is.na(i)); gg$value[i] }

vals <- list(
  nMotifTwoInstancesGaga = gval("T2 Motif2 instances"),
  nGrhMatricesScanned   = as.character(n_grh_matrices),
  grhFimoPctPos         = pc1(a3$pos),
  grhFimoPctNeg         = pc1(a3$neg),
  grhFimoOR             = dec2(a3$or),
  grhFimoP              = dec2(a3$p),
  grhFimoStrictPctPos   = pc1(a4$pos),
  grhFimoStrictPctNeg   = pc1(a4$neg),
  grhFimoStrictOR       = dec2(a4$or),
  grhFimoStrictP        = dec2(a4$p),
  grhFimoDensityRatio   = dec2(dens_ratio),
  grhInstancesOverlapMotifTwo = gv("Motif2 instances overlapping Grh"),
  motifTwoProxOR        = dec2(as.numeric(gv("distal OR Motif2 x promoter"))),
  motifTwoProxP         = sci(2.59e-9),
  motifTwoDistalPctPos  = pc1(as.numeric(gv("distal % Motif2+"))),
  motifTwoDistalPctNeg  = pc1(as.numeric(gv("distal % Motif2-"))),
  nAllHetMotifTwoSum    = as.character(n_sum),
  nAllHetMotifTwoUnique = as.character(n_uni),
  pctAllHetRedundancy   = pc1(100 * (1 - n_uni / n_sum)),
  pctAllHetPairMax      = pc1(pairmax)
)
save_values("motif2_grh", vals)
for (nm in names(vals)) cat(sprintf("  %-28s %s\n", nm, read_value(nm)))
