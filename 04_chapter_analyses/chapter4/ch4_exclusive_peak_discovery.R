## ============================================================
## ch4_exclusive_peak_discovery.R
## Finds peaks unique to one heterodimer joint set, runs repeat-masked STREME
## (--dna --minw 6 --maxw 20 --thresh 0.05) on exclusive sets of >=50 peaks,
## compares motifs to Motif2 PWM by TOMTOM (-dist pearson -min-overlap 5).
## Inputs: data/chapter4/*_Tgo_joint_annotated.csv; repeat-masked FASTAs in
## data/chapter3/ and data/chapter4/. Outputs: data/chapter4/exclusive_peak_discovery/,
## section ch4_exclusive_peaks of results/reported_values.tsv. Root: PROJECT_DIR or cwd.
## ============================================================
PROJECT_DIR <- normalizePath(Sys.getenv("PROJECT_DIR", unset = "."), mustWork = TRUE)
if (!file.exists(file.path(PROJECT_DIR, "README.md")))
  stop("PROJECT_DIR does not point at the repository root: ", PROJECT_DIR)

DATA_DIR      <- file.path(PROJECT_DIR, "data")
CH3_DIR <- file.path(DATA_DIR, "chapter3")
CH4_DIR <- file.path(DATA_DIR, "chapter4")

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
## MEME Suite binaries. Resolution order: $MEME_BIN, then PATH.
MEME_BIN <- Sys.getenv("MEME_BIN", unset = "")
if (!nzchar(MEME_BIN)) MEME_BIN <- dirname(Sys.which("fimo")[[1]])
if (!nzchar(MEME_BIN) || !file.exists(file.path(MEME_BIN, "fimo")))
  stop("MEME Suite not found. Install MEME >= 5.5.5 and set MEME_BIN to its bin/ directory.")
suppressPackageStartupMessages({
  library(GenomicRanges); library(Biostrings); library(BSgenome.Dmelanogaster.UCSC.dm6)
})
OUT <- file.path(CH4_DIR, "exclusive_peak_discovery"); dir.create(OUT, showWarnings = FALSE)
PWM <- file.path(CH3_DIR, "STRCACTGARARAAA.meme")
H <- c("Trh", "Sim", "Dys", "Sima")
gr <- lapply(setNames(H, H), function(h) {
  d <- read.csv(file.path(CH4_DIR, paste0(h, "_Tgo_joint_annotated.csv")))
  GRanges(d$seqnames, IRanges(d$start, d$end)) })
mfa <- c(Trh  = file.path(CH3_DIR, "Trh_Tgo_joint_peaks_filt_repeatmasked.fa"),
         Sim  = file.path(CH4_DIR, "Sim_Tgo_joint_peaks_filt_repeatmasked.fa"),
         Dys  = file.path(CH4_DIR, "Dys_Tgo_joint_peaks_filt_repeatmasked.fa"),
         Sima = file.path(CH4_DIR, "Sima_Tgo_joint_peaks_filt_repeatmasked.fa"))
res <- list(); vals <- list()
for (h in H) {
  others <- reduce(do.call(c, unname(gr[setdiff(H, h)])))
  ex <- countOverlaps(gr[[h]], others, ignore.strand = TRUE) == 0L
  vals[[paste0("excl", h, "N")]] <- sum(ex)
  vals[[paste0("excl", h, "Pct")]] <- formatC(100 * mean(ex), format = "f", digits = 1)
  # the masked FASTA is in the same order as the peak table: verify at unmasked positions
  m <- readDNAStringSet(mfa[[h]])
  stopifnot(length(m) == length(gr[[h]]))
  u <- getSeq(BSgenome.Dmelanogaster.UCSC.dm6, gr[[h]])
  same <- mapply(function(a, b) { a <- as.character(a); b <- as.character(b)
    if (nchar(a) != nchar(b)) return(FALSE); ai <- strsplit(a, "")[[1]]; bi <- strsplit(b, "")[[1]]
    all(ai[ai != "N"] == toupper(bi[ai != "N"])) }, m, u)
  stopifnot(all(same))
  if (sum(ex) < 50L) { res[[h]] <- data.frame(set = h, n_excl = sum(ex), streme = "not run (< 50 peaks)"); next }
  fa <- file.path(OUT, paste0(h, "_exclusive_masked.fa")); writeXStringSet(m[ex], fa)
  so <- file.path(OUT, paste0(h, "_exclusive_streme"))
  system2(file.path(MEME_BIN, "streme"), c("--p", fa, "--oc", so, "--dna", "--minw", "6", "--maxw", "20", "--thresh", "0.05"),
          stdout = FALSE, stderr = FALSE)
  tt <- system2(file.path(MEME_BIN, "tomtom"), c("-no-ssc", "-min-overlap", "5", "-dist", "pearson", "-thresh", "1", "-text",
                file.path(so, "streme.txt"), PWM), stdout = TRUE, stderr = FALSE)
  tt <- read.delim(text = paste(tt, collapse = "\n"), stringsAsFactors = FALSE)
  tt <- tt[!grepl("^#", tt$Query_ID) & nzchar(tt$Query_ID), ]
  st <- readLines(file.path(so, "streme.txt")); hd <- grep("^MOTIF ", st)
  ids <- sub("^MOTIF (\\S+).*$", "\\1", st[hd]); E <- sub(".*E= *([0-9.eE+-]+).*$", "\\1", st[hd + 1])
  best <- if (nrow(tt)) tt[which.min(tt$q.value), ] else NULL
  rk <- if (!is.null(best)) match(best$Query_ID, ids) else NA
  res[[h]] <- data.frame(set = h, n_excl = sum(ex), streme = sprintf("%d motifs", length(ids)),
                         best_match = if (is.null(best)) NA else best$Query_ID, rank = rk,
                         streme_E = if (is.na(rk)) NA else E[rk], tomtom_q = if (is.null(best)) NA else best$q.value)
  vals[[paste0("excl", h, "Rank")]] <- if (is.na(rk)) "NA" else rk
  vals[[paste0("excl", h, "E")]] <- if (is.na(rk)) "NA" else formatC(as.numeric(E[rk]), format = "fg", digits = 2)
  vals[[paste0("excl", h, "Q")]] <- if (is.null(best)) "NA" else formatC(best$q.value, format = "e", digits = 1)
  vals[[paste0("excl", h, "Consensus")]] <- if (is.null(best)) "NA" else sub("^[0-9]+-", "", best$Query_ID)
}
R <- do.call(rbind, lapply(res, function(d) { for (k in c("best_match", "rank", "streme_E", "tomtom_q")) if (!k %in% names(d)) d[[k]] <- NA; d }))
print(R, row.names = FALSE)
write.csv(R, file.path(OUT, "exclusive_peak_discovery.csv"), row.names = FALSE)
save_values("ch4_exclusive_peaks", vals)
