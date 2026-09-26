## ============================================================
##  ch3_peak_annotation_counts.R
##  Computes chapter 3 counts: joint-peak participation per factor, nested F-test of Motif2 regression, gene-set sizes for ChIP-Enrich, Motif2 fold contrasts between regulatory classes, minimum detectable differences for unique-peak comparisons, and composition of MEME-2 (Motif2) site table.
##  Run after chapter3_trh_tgo_binding.Rmd; values recorded in results/reported_values.tsv.
##  Inputs: data/final_peaks/TRH_I

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


AX  <- file.path(PROJECT_DIR, "data", "chapter3")
FIN <- file.path(PROJECT_DIR, "data/final_peaks")

fmt_p <- function(p, digits = 2) {
  s <- formatC(p, format = "e", digits = digits)
  q <- strsplit(s, "e")[[1]]
  paste0(q[1], "e", as.integer(q[2]))
}

vals <- list()

## ------------------------------------------------------------
##  Joint-peak participation. The 614 joint peaks are pintersect intervals;
##  Trh and Tgo need not contribute the same number of peaks to them, because
##  one broad Tgo peak can overlap two adjacent Trh peaks. Count distinct
##  peaks per side.
## ------------------------------------------------------------
rd <- function(p) {
  x <- read.delim(p, header = FALSE, stringsAsFactors = FALSE)
  names(x)[1:3] <- c("chr", "start", "end")
  good <- c("chr2L", "chr2R", "chr3L", "chr3R", "chr4", "chrX", "chrY", "chrM")
  x <- x[x$chr %in% good, ]
  x[!duplicated(x[, 1:3]), ]
}
trh <- rd(file.path(FIN, "TRH_IDR0.05.clean.narrowPeak"))
tgo <- rd(file.path(FIN, "TGO_IDR0.05.clean.narrowPeak"))

ai <- integer(0); bi <- integer(0)
for (cc in intersect(unique(trh$chr), unique(tgo$chr))) {
  a <- which(trh$chr == cc); b <- which(tgo$chr == cc)
  for (i in a) {
    j <- b[tgo$start[b] < trh$end[i] & tgo$end[b] > trh$start[i]]
    if (length(j)) { ai <- c(ai, rep(i, length(j))); bi <- c(bi, j) }
  }
}
pin <- data.frame(chr = trh$chr[ai],
                  start = pmax(trh$start[ai], tgo$start[bi]),
                  end   = pmin(trh$end[ai],   tgo$end[bi]))
pin <- pin[!duplicated(pin), ]

trh_in <- length(unique(ai))     # distinct Trh peaks in the overlap
tgo_in <- length(unique(bi))     # distinct Tgo peaks in the overlap

stopifnot(nrow(pin) == 614L,          # canonical joint-interval count
          nrow(trh) == 1033L, nrow(tgo) == 1403L,
          trh_in == 614L)             # Trh side happens to equal it; Tgo does not

## Every Trh peak has exactly one Tgo partner; the deficit on the Tgo side
## is entirely due to broad Tgo peaks that span two adjacent Trh peaks.
stopifnot(all(table(ai) == 1L))
n_span <- sum(table(bi) > 1L)
stopifnot(nrow(pin) - tgo_in == n_span)

vals$trhJointPeaks   <- trh_in
vals$tgoJointPeaks   <- tgo_in
vals$tgoSpanningPeaks <- n_span
vals$trhPctJoint   <- round(100 * trh_in / nrow(trh), 1)
vals$tgoPctJoint   <- round(100 * tgo_in / nrow(tgo), 1)
vals$tgoUniqueTotalN <- nrow(tgo) - tgo_in                 # Tgo-unique peaks (tgo_unique.fa)

## ------------------------------------------------------------
##  Motif2 dosage regression: coefficient p-value and nested F-test.
## ------------------------------------------------------------
t1 <- readRDS(file.path(AX, "test1_regression_binding.rds"))
cf <- summary(t1$full)$coefficients
av <- anova(t1$reduced, t1$full)

## OLS coefficient p and the nested-F p are the same test on 1 df; assert it
## rather than quoting two numbers that could drift apart.
stopifnot(abs(cf["n_motif2", "Pr(>|t|)"] - av$`Pr(>F)`[2]) < 1e-12)

vals$motifTwoRegP    <- fmt_p(cf["n_motif2", "Pr(>|t|)"])
vals$motifTwoRegFstat <- round(av$F[2], 2)
vals$motifTwoRegDfRes <- av$Res.Df[2]

## ------------------------------------------------------------
##  Gene-set sizes, read from the gene-set file ChIP-Enrich was given, and the
##  numbers of genes ChIP-Enrich modelled after applying its locus definition
##  (the denominators for the peak-gene counts).
## ------------------------------------------------------------
gs  <- read.delim(file.path(AX, "geneset_tracheal_vs_complement.tsv"),
                  header = FALSE, stringsAsFactors = FALSE)
tab <- table(gs$V1)
ce  <- read.csv(file.path(AX, "chipenrich_ch3_custom.csv"), stringsAsFactors = FALSE)
rownames(ce) <- ce$Geneset.ID

vals$tsSgEntrezN      <- unname(tab[["tracheal_custom"]])          # 2049
vals$compEntrezN      <- unname(tab[["All_dm6_minus_tracheal"]])
vals$ceTrachModelledN <- ce["tracheal_custom",        "N.Geneset.Genes"]
vals$ceCompModelledN  <- ce["All_dm6_minus_tracheal", "N.Geneset.Genes"]

## ------------------------------------------------------------
##  Motif2 carriage fold contrasts between regulatory classes, each contrast
##  reported separately (conservative criterion counts per class).
## ------------------------------------------------------------
pct <- c(distEnh = 60 / 114, promEnh = 20 / 55,
         promNon = 86 / 413, distNon = 12 / 32)
vals$motifTwoFoldEnhVsProm         <- round(pct[["distEnh"]] / pct[["promNon"]], 1)
vals$motifTwoFoldEnhVsDist         <- round(pct[["distEnh"]] / pct[["distNon"]], 1)
vals$motifTwoFoldPromEnhVsPromNon  <- round(pct[["promEnh"]] / pct[["promNon"]], 1)
vals$motifTwoFoldPromEnhVsDistNon  <- round(pct[["promEnh"]] / pct[["distNon"]], 2)

## ------------------------------------------------------------
##  Minimum detectable difference at 80% power, two-sided, for the
##  two unique-peak comparisons. The binding constraint is the smaller
##  comparator (Trh-unique, n = 419).
## ------------------------------------------------------------
mde <- function(n1, n2, p1, alpha = 0.05, target = 0.80) {
  pw <- function(dd) {
    p2 <- p1 + dd
    pbar <- (n1 * p1 + n2 * p2) / (n1 + n2)
    pnorm((dd - qnorm(1 - alpha / 2) *
             sqrt(pbar * (1 - pbar) * (1 / n1 + 1 / n2))) /
          sqrt(p1 * (1 - p1) / n1 + p2 * (1 - p2) / n2))
  }
  uniroot(function(dd) pw(dd) - target, c(0.005, 0.40))$root
}
p_joint <- 221 / 614
vals$mdeUniqueTrhPp   <- round(100 * mde(614, 419, p_joint), 1)
vals$mdeUniqueTgoPp   <- round(100 * mde(614, 794, p_joint), 1)
vals$mdeUniqueArmsPp  <- round(100 * max(mde(614, 419, p_joint),
                                         mde(614, 794, p_joint)), 1)

## ------------------------------------------------------------
##  Core composition, parsed from the MEME-2 site table: number of invariant
##  columns in positions 4-9 and the per-position consensus frequencies.
## ------------------------------------------------------------
mt  <- readLines(file.path(AX, "memechip_Trh_Tgo_masked/meme_out/meme.txt"))
i0  <- grep("Motif GTRCACTGARARAAA MEME-2 sites sorted by position p-value", mt)[1]
hdr <- grep("^Sequence name\\s+Strand\\s+Start\\s+P-value", mt)
hdr <- hdr[hdr > i0][1]

k <- hdr + 2L; body <- character(0)
while (k <= length(mt) && !grepl("^-{30,}", mt[k])) { body <- c(body, mt[k]); k <- k + 1L }
body <- body[nzchar(trimws(body))]

## The site column is the only whitespace-delimited token of exactly the
## motif width; flanks are 10 bp and the name/strand/start/p fields are
## never 15 characters of ACGTN. Assert one hit per row rather than
## slicing by column offset, which breaks on N-padded rows.
sites <- vapply(strsplit(trimws(body), "\\s+"), function(v) {
  cand <- v[nchar(v) == 15L & grepl("^[ACGTN]+$", v)]
  stopifnot(length(cand) == 1L)
  cand
}, character(1))

M    <- do.call(rbind, strsplit(sites, ""))
freq <- apply(M, 2, function(x) max(table(x)) / length(x))

vals$motifTwoSites        <- nrow(M)
vals$motifTwoCoreInvN     <- sum(freq[4:9] == 1)              # invariant cols in 4-9
vals$motifTwoCoreFlexLoPct <- round(100 * min(freq[6:8]), 1)  # weakest of 6-8
vals$motifTwoCoreFlexHiPct <- round(100 * max(freq[6:8]), 1)
vals$motifTwoCorePosNinePct <- round(100 * freq[9], 1)
vals$motifTwoCoreExactN   <- sum(substr(sites, 4, 9) == "CACTGA")

stopifnot(vals$motifTwoCoreInvN == 2L)   # positions 4 and 5 only

## ---- Record ----
save_values("ch3", vals)
for (nm in names(vals)) cat(sprintf("  %-28s %s\n", nm, as.character(vals[[nm]])))
