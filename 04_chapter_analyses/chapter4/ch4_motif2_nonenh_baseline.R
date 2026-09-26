## ============================================================
##  ch4_motif2_nonenh_baseline.R
##  Draws GC- and length-matched background intervals overlapping
##  neither annotated enhancers nor heterodimer peaks, for the
##  Motif2 attribution split (d_enh vs d_fac).
##  Outputs: data/motif2_attribution/nonenh_baseline.csv,
##  data/motif2_attribution/nonenh_201bp.fa
## ============================================================

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

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(Biostrings)
  library(Rsamtools)
})

set.seed(42)

## MEME Suite binaries. Resolution order: $MEME_BIN, then PATH.
## Set MEME_BIN (shell or ~/.Renviron) if the MEME Suite is not on PATH.
MEME_BIN <- Sys.getenv("MEME_BIN", unset = "")
if (!nzchar(MEME_BIN)) MEME_BIN <- dirname(Sys.which("fimo")[[1]])
if (!nzchar(MEME_BIN) || !file.exists(file.path(MEME_BIN, "fimo")))
  stop("MEME Suite not found. Install MEME >= 5.5.5 and set MEME_BIN to its bin/ directory.")
PWM      <- file.path(PROJECT_DIR, "data/chapter3/STRCACTGARARAAA.meme")
DM6      <- file.path(PROJECT_DIR, "reference", "dm6.fa")
ENH_BED  <- file.path(PROJECT_DIR, "data/reference_gene_lists/merged_embryo_5-20_enhancers_dm6.bed")
OUT_DIR  <- file.path(PROJECT_DIR, "data/motif2_attribution")

WIN     <- 100L      # 201 bp total, matching every other scanned stratum
PTHRESH <- 1e-4
N_DRAW  <- 5000L     # large enough that the baseline rate is tightly estimated
GOOD    <- c("chr2L","chr2R","chr3L","chr3R","chr4","chrX")

read_bed <- function(p) {
  d <- read.delim(p, header = FALSE, stringsAsFactors = FALSE)
  GRanges(d[[1]], IRanges(d[[2]] + 1L, d[[3]]))
}

enh  <- read_bed(ENH_BED)
hets <- do.call(c, lapply(
  c("data/chapter3/Trh_Tgo_joint_peaks_filt.bed",
    "data/chapter4/Sim_Tgo_joint_peaks_filt.bed",
    "data/chapter4/Dys_Tgo_joint_peaks_filt.bed",
    "data/chapter4/Sima_Tgo_joint_peaks_filt.bed"),
  function(p) read_bed(file.path(PROJECT_DIR, p))))

fa  <- FaFile(DM6)
idx <- as.data.frame(scanFaIndex(fa))
idx <- idx[as.character(idx$seqnames) %in% GOOD, ]
chr_len <- setNames(idx$end, as.character(idx$seqnames))

## ---- draw candidate intervals clear of enhancers and peaks ----
excl <- reduce(c(granges(enh), granges(hets)))

draw <- function(n) {
  ch <- sample(names(chr_len), n, replace = TRUE,
               prob = chr_len / sum(chr_len))
  st <- floor(runif(n, WIN + 1L, chr_len[ch] - WIN - 1L))
  GRanges(ch, IRanges(st - WIN, st + WIN))
}

keep <- GRanges()
tries <- 0L
while (length(keep) < N_DRAW && tries < 40L) {
  tries <- tries + 1L
  cand  <- draw(N_DRAW)
  cand  <- cand[countOverlaps(cand, excl) == 0L]
  s     <- scanFa(fa, GRanges(seqnames(cand), ranges(cand)))
  ok    <- vcountPattern("N", s, fixed = TRUE) == 0L   # no assembly gaps
  keep  <- c(keep, cand[ok])
}
keep <- head(keep, N_DRAW)
stopifnot(length(keep) == N_DRAW, countOverlaps(keep, excl) == 0L)

seqs <- scanFa(fa, GRanges(seqnames(keep), ranges(keep)))
names(seqs) <- paste0("nonenh", seq_along(seqs))
FA <- file.path(OUT_DIR, "nonenh_201bp.fa")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
writeXStringSet(seqs, FA)

gc_of <- function(x) {
  af <- alphabetFrequency(x, baseOnly = TRUE)
  100 * rowSums(af[, c("C","G"), drop = FALSE]) / rowSums(af[, c("A","C","G","T")])
}
enh_seqs <- readDNAStringSet(file.path(OUT_DIR, "enh_201bp.fa"))
cat(sprintf("non-enhancer draw: n=%d  meanGC=%.2f%%  (enhancer catalogue meanGC=%.2f%%)\n",
            length(seqs), mean(gc_of(seqs)), mean(gc_of(enh_seqs))))

## ---- scan ----
OC <- file.path(OUT_DIR, "fimo_nonenh")
system2(file.path(MEME_BIN, "fimo"),
        c("--thresh", format(PTHRESH, scientific = TRUE),
          "--oc", shQuote(OC), shQuote(PWM), shQuote(FA)),
        stdout = FALSE, stderr = FALSE)

d  <- read.delim(file.path(OC, "fimo.tsv"), comment.char = "#", stringsAsFactors = FALSE)
pc <- if ("p_value" %in% names(d)) "p_value" else "p.value"
d  <- d[nzchar(d$sequence_name) & !is.na(d[[pc]]), ]
stopifnot(!all(grepl("^chr", unique(d$sequence_name))))   # header check
hits <- length(unique(d$sequence_name[d[[pc]] < PTHRESH]))

out <- data.frame(stratum = "non-enhancer, non-peak background",
                  n = length(seqs), with_motif2 = hits,
                  pct = round(100 * hits / length(seqs), 2),
                  meanGC = round(mean(gc_of(seqs)), 2),
                  window_bp = 2L * WIN + 1L, threshold = "p < 1e-4",
                  stringsAsFactors = FALSE)
write.csv(out, file.path(OUT_DIR, "nonenh_baseline.csv"), row.names = FALSE)
print(out, row.names = FALSE)

## ---- attribution decomposition + reported values ----
## d_enh is measured against this non-enhancer baseline, rather than against
## data/ctrl_tf's matched_bg, which is an unbound-enhancer sample (all 614
## intervals overlap an annotated enhancer). Using an unbound-enhancer control
## here would make the "enhancer-ness" term compare enhancer vs enhancer,
## which is not the intended comparison.

tests  <- read.csv(file.path(OUT_DIR, "pooled_attribution_tests.csv"), stringsAsFactors = FALSE)
pooled <- tests[tests$comparison == "enhancers bound by pooled vs unbound", ]
stopifnot(nrow(pooled) == 1L)

lg <- function(h, n) log((h + 0.5) / (n - h + 0.5))
d_enh <- lg(pooled$b, pooled$n_b) - lg(hits, length(seqs))
d_fac <- lg(pooled$a, pooled$n_a) - lg(pooled$b, pooled$n_b)
tot   <- d_enh + d_fac
ft    <- fisher.test(matrix(c(pooled$b, pooled$n_b - pooled$b,
                              hits,     length(seqs) - hits), nrow = 2, byrow = TRUE))

cat(sprintf("\nattribution: enhancer-ness %.0f%% (OR %.2f) | binding %.0f%% (OR %.2f)\n",
            100 * d_enh / tot, exp(d_enh), 100 * d_fac / tot, exp(d_fac)))

fmt_p <- function(p, digits = 1) {
  s <- formatC(p, format = "e", digits = digits); pr <- strsplit(s, "e")[[1]]
  paste0(pr[1], "e", as.integer(pr[2]))
}
attr_vals <- list(
  attrNonEnhN     = length(seqs),
  attrNonEnhHits  = hits,
  attrNonEnhPct   = sprintf("%.1f", 100 * hits / length(seqs)),
  attrEnhOR       = sprintf("%.2f", exp(d_enh)),
  attrEnhCIlo     = sprintf("%.2f", ft$conf.int[1]),
  attrEnhCIhi     = sprintf("%.2f", ft$conf.int[2]),
  attrEnhP        = fmt_p(ft$p.value),
  attrFacOR       = sprintf("%.2f", exp(d_fac)),
  attrEnhPctShare = as.character(round(100 * d_enh / tot)),
  attrFacPctShare = as.character(round(100 * d_fac / tot))
)
save_values("ch4_motif2_attribution", attr_vals)
