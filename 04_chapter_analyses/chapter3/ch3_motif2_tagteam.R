#!/usr/bin/env Rscript
# ============================================================
# ch3_motif2_tagteam.R: scans joint-peak FASTA sequences for the Zelda TAGteam heptamer, checks co-occurrence with Motif2.
# i-cisTarget scored Zelda zero in all eight sets.
# Usage: Rscript 04_chapter_analyses/chapter3/ch3_motif2_tagteam.R
# Repository root: env var PROJECT_DIR, or working directory.
# ============================================================
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

suppressPackageStartupMessages(library(Biostrings))
AUX     <- DATA_DIR
OUT     <- file.path(AUX, "motif2_tagteam")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
N_JOINT <- 614L
TAGTEAM <- "CAGGTAG"            # canonical Zelda heptamer

## ---- 1. joint-peak sequences ---------------------------------------------
fa <- readDNAStringSet(file.path(AUX, "chapter3/Trh_Tgo_joint_peaks_filt_raw.fa"))
stopifnot(length(fa) == N_JOINT)
pk_of <- as.integer(sub("^peak_", "", names(fa)))
stopifnot(!anyNA(pk_of), max(pk_of) == N_JOINT)

## ---- 2. TAGteam occurrences, both strands --------------------------------
hits <- lapply(seq_along(fa), function(i) {
  f <- start(matchPattern(TAGTEAM, fa[[i]], fixed = TRUE))
  r <- start(matchPattern(reverseComplement(DNAString(TAGTEAM)), fa[[i]], fixed = TRUE))
  sort(c(f, r))
})
names(hits) <- pk_of
tt_n  <- vapply(hits, length, integer(1))
tt_pk <- pk_of[tt_n > 0]

## ---- 3. Motif2 instances (conservative q < 0.05 scan) -------------------
m2 <- read.delim(file.path(AUX, "chapter3/fimo_joint/fimo.tsv"),
                 comment.char = "#", stringsAsFactors = FALSE)
m2 <- m2[!is.na(m2$start) & m2$q.value < 0.05 &
           m2$motif_id == "STRCACTGARARAAA", , drop = FALSE]
m2$peak <- as.integer(sub("^peak_", "", m2$sequence_name))
stopifnot(max(m2$peak) <= N_JOINT)
m2_pk <- sort(unique(m2$peak))

## ---- 4. peak-level 2x2 ---------------------------------------------------
both <- length(intersect(m2_pk, tt_pk))
tab  <- matrix(c(both, length(setdiff(m2_pk, tt_pk)),
                 length(setdiff(tt_pk, m2_pk)),
                 N_JOINT - both - length(setdiff(m2_pk, tt_pk)) - length(setdiff(tt_pk, m2_pk))),
               2, 2, byrow = TRUE,
               dimnames = list(c("Motif2+","Motif2-"), c("TAGteam+","TAGteam-")))
stopifnot(sum(tab) == N_JOINT)
ft <- fisher.test(tab)

## ---- 5. instance-level separation ----------------------------------------
ov <- 0L; sep <- numeric(0)
for (i in seq_len(nrow(m2))) {
  h <- hits[[as.character(m2$peak[i])]]
  if (!length(h)) next
  if (any(h <= m2$stop[i] & (h + nchar(TAGTEAM) - 1L) >= m2$start[i])) ov <- ov + 1L
  sep <- c(sep, min(abs(h - m2$start[i])))
}

## ---- 6. outputs ----------------------------------------------------------
write.csv(data.frame(peak = pk_of, tagteam_sites = as.integer(tt_n),
                     motif2 = pk_of %in% m2_pk),
          file.path(OUT, "tagteam_per_peak.csv"), row.names = FALSE)

sci2 <- function(p) { e <- floor(log10(p)); sprintf("%.1fe%d", p/10^e, e) }
vals <- list(
  tagteamWord            = TAGTEAM,
  nTagteamPosPeaks       = sum(tt_n > 0),
  pctTagteamPosPeaks     = sprintf("%.1f", 100 * sum(tt_n > 0) / N_JOINT),
  nTagteamSites          = sum(tt_n),
  nMotifTwoInTagteamPos  = tab[1,1],
  pctMotifTwoInTagteamPos= sprintf("%.1f", 100 * tab[1,1] / sum(tab[,1])),
  pctMotifTwoInTagteamNeg= sprintf("%.1f", 100 * tab[1,2] / sum(tab[,2])),
  tagteamMotifTwoOR      = sprintf("%.2f", unname(ft$estimate)),
  tagteamMotifTwoORloCI  = sprintf("%.2f", ft$conf.int[1]),
  tagteamMotifTwoORhiCI  = sprintf("%.2f", ft$conf.int[2]),
  tagteamMotifTwoP       = sprintf("%.2f", ft$p.value),
  nMotifTwoOverlapTagteam= ov)
if (length(sep)) vals$motifTwoTagteamSepMedian <- sprintf("%.0f", median(sep))
save_values("ch3_motif2_tagteam", vals)

cat(sprintf("\nTAGteam (%s, both strands) in the %d joint peaks\n", TAGTEAM, N_JOINT))
cat(sprintf("  peaks with >=1 site: %d (%.1f%%);  total sites: %d\n",
            sum(tt_n > 0), 100*sum(tt_n>0)/N_JOINT, sum(tt_n)))
print(tab)
cat(sprintf("Fisher OR = %.3f (95%% CI %.3f-%.3f), p = %.4g\n",
            ft$estimate, ft$conf.int[1], ft$conf.int[2], ft$p.value))
cat(sprintf("Motif2 carriage: %.1f%% of TAGteam+ vs %.1f%% of TAGteam- peaks\n",
            100*tab[1,1]/sum(tab[,1]), 100*tab[1,2]/sum(tab[,2])))
cat(sprintf("instances: %d Motif2; %d overlap a TAGteam site\n", nrow(m2), ov))
if (length(sep)) cat(sprintf("co-occurring in one peak: n = %d, median separation %.0f bp (range %.0f-%.0f)\n",
                             length(sep), median(sep), min(sep), max(sep)))
cat("written:", OUT, "+ results/reported_values.tsv (section ch3_motif2_tagteam)\n")
