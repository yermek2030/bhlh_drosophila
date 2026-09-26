## ============================================================
## joint_peak_count.R
## Rebuilds Trh/Tgo joint peak set from
PROJECT_DIR <- normalizePath(Sys.getenv("PROJECT_DIR", unset = "."), mustWork = TRUE)
if (!file.exists(file.path(PROJECT_DIR, "README.md")))
  stop("PROJECT_DIR does not point at the repository root: ", PROJECT_DIR)

DATA_DIR      <- file.path(PROJECT_DIR, "data")
FINAL_DIR       <- file.path(DATA_DIR, "final_peaks")
trh_narrowPeak_file  <- file.path(FINAL_DIR, "TRH_IDR0.05.final.narrowPeak")
tgo_narrowPeak_file  <- file.path(FINAL_DIR, "TGO_IDR0.05.final.narrowPeak")
sim_narrowPeak_file  <- file.path(FINAL_DIR, "SIM_IDR0.05.final.narrowPeak")
dys_narrowPeak_file  <- file.path(FINAL_DIR, "DYS_IDR0.05.final.narrowPeak")
sima_narrowPeak_file <- file.path(FINAL_DIR, "SIMA_IDR0.05.final.narrowPeak")

good_chr <- c("chr2L", "chr2R", "chr3L", "chr3R", "chr4", "chrX")   # dm6 major arms
suppressPackageStartupMessages({
  library(GenomicRanges); library(GenomeInfoDb); library(rtracklayer)
})

extraCols_narrowPeak <- c(signalValue = "numeric", pValue = "numeric",
                          qValue = "numeric", peak = "integer")

clean_peak_set <- function(gr) {
  gr <- gr[seqnames(gr) %in% good_chr]
  gr <- keepSeqlevels(gr, good_chr, pruning.mode = "coarse")
  gr[!duplicated(gr)]
}
read_peaks <- function(f) clean_peak_set(import(f, format = "BED",
                                                 extraCols = extraCols_narrowPeak))
joint_set <- function(a, b) {
  ov <- findOverlaps(a, b, minoverlap = 1, ignore.strand = TRUE)
  j  <- pintersect(a[queryHits(ov)], b[subjectHits(ov)])
  j  <- keepSeqlevels(j, good_chr, pruning.mode = "coarse")
  j[!duplicated(j)]
}

tgo  <- read_peaks(tgo_narrowPeak_file)
sets <- c(Trh = trh_narrowPeak_file, Sim = sim_narrowPeak_file,
          Dys = dys_narrowPeak_file, Sima = sima_narrowPeak_file)
expected <- c(Trh = 614L, Sim = 355L, Dys = 870L, Sima = 296L)

n <- vapply(names(sets), function(h) length(joint_set(read_peaks(sets[[h]]), tgo)),
            integer(1))
print(data.frame(heterodimer = paste0(names(n), "/Tgo"), joint_peaks = n,
                 reported = expected[names(n)], row.names = NULL))
stopifnot(identical(n[["Trh"]], 614L))
if (!all(n == expected[names(n)])) warning("a joint-set size differs from the text")
cat("Trh/Tgo joint peak set: 614 peaks -- OK\n")
