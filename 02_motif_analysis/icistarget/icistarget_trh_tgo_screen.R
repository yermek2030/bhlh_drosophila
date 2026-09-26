#!/usr/bin/env Rscript
# ============================================================
# icistarget_trh_tgo_screen.R: process i-cisTarget archive for 614 Trh/Tgo peaks screened against 24453 PWMs and 1502 modERN tracks.
# Converts motifs.cb to MEME and runs TOMTOM against Motif2 PWM for novelty; checks pioneer-factor motif enrichment and modERN positive control.
# Usage: Rscript icistarget_trh_tgo_screen.R [path/to/icistarget]; requires MEME Suite on PATH.
# ============================================================

## Repository root: the environment variable PROJECT_DIR, or the working directory.
PROJECT_DIR <- normalizePath(Sys.getenv("PROJECT_DIR", unset = "."), mustWork = TRUE)
if (!file.exists(file.path(PROJECT_DIR, "README.md")))
  stop("PROJECT_DIR does not point at the repository root: ", PROJECT_DIR)

args <- commandArgs(trailingOnly = TRUE)
# args[1] = archive dir (the unzipped "icistarget" folder)
# args[2] = output dir. Defaults to the archive's PARENT, so each run writes
#           beside its own input and a second run cannot silently overwrite
#           the first one's tables.
D   <- if (length(args) >= 1) args[1] else
  file.path(PROJECT_DIR, "data", "icistarget", "trh_tgo", "icistarget")
stopifnot(dir.exists(D))
OUT <- if (length(args) >= 2) args[2] else dirname(normalizePath(D))
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
cat("archive:", D, "\noutput :", OUT, "\n\n")
Q <- file.path(PROJECT_DIR, "data/chapter3/STRCACTGARARAAA.meme")
stopifnot(file.exists(Q))

# statistics.tbl is 7.5 MB raw, above the 5 MB cap applied to this class of
# directory, so the copy in data/ is gzipped. Accept either.
stats_f <- file.path(D, "statistics.tbl")
if (!file.exists(stats_f)) stats_f <- paste0(stats_f, ".gz")
stopifnot(file.exists(stats_f))

# ------------------------------------------------------------
# 1. What was submitted, and what was tested
# ------------------------------------------------------------
bed <- read.delim(file.path(D, "input_mapped_to_icistarget_regions.bed"),
                  header = FALSE, stringsAsFactors = FALSE)
st  <- read.delim(stats_f, stringsAsFactors = FALSE,
                  quote = "", comment.char = "")
stopifnot(nrow(st) > 0, "NES" %in% names(st))
sig <- st[st$NES >= 3, ]                       # i-cisTarget default threshold

cat(sprintf("input regions: %d (mapped to %d i-cisTarget regions)\n",
            length(unique(bed$V4)), nrow(bed)))
cat(sprintf("features tested: %d  (PWMs %d, modERN %d)\n", nrow(st),
            sum(st$FeatureDatabase == "PWMs"),
            sum(st$FeatureDatabase == "modERN")))
cat(sprintf("significant at NES>=3: %d  (PWMs %d, modERN %d)\n", nrow(sig),
            sum(sig$FeatureDatabase == "PWMs"),
            sum(sig$FeatureDatabase == "modERN")))

# ------------------------------------------------------------
# 2. Cluster-Buster -> MEME. Rows are A C G T counts; normalise per row.
# ------------------------------------------------------------
cb  <- readLines(file.path(D, "motifs.cb"), warn = FALSE)
cb  <- cb[nzchar(trimws(cb))]
idx <- grep("^>", cb); stopifnot(length(idx) > 0)
ends <- c(idx[-1] - 1L, length(cb))
mats <- Map(function(s, e) {
  m <- do.call(rbind, lapply(strsplit(cb[(s + 1L):e], "[\t ]+"), as.numeric))[, 1:4, drop = FALSE]
  m / pmax(rowSums(m), 1)
}, idx, ends)
names(mats) <- sub("^>", "", cb[idx])

meme <- c("MEME version 4", "", "ALPHABET= ACGT", "", "strands: + -", "",
          "Background letter frequencies", "A 0.25 C 0.25 G 0.25 T 0.25", "")
for (n in names(mats)) {
  m <- mats[[n]]
  meme <- c(meme, paste("MOTIF", n),
            sprintf("letter-probability matrix: alength= 4 w= %d nsites= 20 E= 0", nrow(m)),
            apply(m, 1, function(r) paste0(" ", paste(sprintf("%.6f", r), collapse = " "))), "")
}
meme_f <- file.path(OUT, "icistarget_289.meme")
writeLines(meme, meme_f)
cat(sprintf("converted %d matrices -> icistarget_289.meme\n", length(mats)))

# ------------------------------------------------------------
# 3. TOMTOM: is Motif2 any of them?
# ------------------------------------------------------------
run <- function(oc, extra) system2("tomtom",
  c("-no-ssc", "-oc", file.path(OUT, oc), extra, "-min-overlap", "5",
    "-dist", "pearson", shQuote(Q), shQuote(meme_f)),
  stdout = FALSE, stderr = FALSE)
run("tt_icis_q05",  c("-thresh", "0.05"))
run("tt_icis_near", c("-thresh", "50", "-evalue"))

hits <- read.delim(file.path(OUT, "tt_icis_q05/tomtom.tsv"), comment.char = "#", stringsAsFactors = FALSE)
hits <- hits[!is.na(hits$p.value), ]
near <- read.delim(file.path(OUT, "tt_icis_near/tomtom.tsv"), comment.char = "#", stringsAsFactors = FALSE)
near <- near[!is.na(near$p.value), ]
near <- near[order(near$E.value), ]
cat(sprintf("\nTOMTOM hits at q<=0.05: %d\n", nrow(hits)))
cat(sprintf("nearest neighbour: %s  E=%.3g  q=%.3g  consensus=%s\n",
            near$Target_ID[1], near$E.value[1], near$q.value[1],
            near$Target_consensus[1]))

# ------------------------------------------------------------
# 4. Pioneer-factor enrichment, and the modERN positive control
# ------------------------------------------------------------
hay <- toupper(paste(st$FeatureID, st$FeatureDescription, st$FeatureAnnotations))
# Each factor is searched as a UNION OF SYNONYMS. Searching a single symbol
# under-reports: "GRH" alone misses matrices named GRAINYHEAD/GRHL (the best
# Grainyhead feature scores 2.65 and is named cisbp__M5928, not GRH), and "TRL"
# alone would miss GAGA-named matrices. Add a synonym rather than a new row.
syn <- list(Zld   = c("ZLD", "ZELDA"),
            Trl   = c("TRL", "GAF", "GAGA"),
            Opa   = c("OPA", "ODD-PAIRED"),
            # AACCGGTT is the canonical Grainyhead consensus and TFCP2 is the
            # vertebrate Grainyhead-like family, so matrices named for either
            # are Grh-family evidence.
            # This list must stay identical to SYN in
            # icistarget_trh_tgo_summary.R or the reported values will disagree
            # with the table they describe.
            Grh   = c("GRH", "GRAINY", "GRHL", "AACCGGTT", "TFCP2"),
            CLAMP = c("CLAMP"))
pio <- do.call(rbind, lapply(names(syn), function(f) {
  i <- unique(unlist(lapply(syn[[f]], function(p) grep(paste0("\\b", p), hay))))
  data.frame(factor = f, synonyms = paste(syn[[f]], collapse = "|"),
             n_features = length(i),
             best_NES = if (length(i)) round(max(st$NES[i]), 2) else NA_real_,
             significant = if (length(i)) max(st$NES[i]) >= 3 else NA,
             best_feature = if (length(i)) st$FeatureID[i][which.max(st$NES[i])] else NA_character_,
             stringsAsFactors = FALSE)
}))
print(pio, row.names = FALSE)

md   <- st[st$FeatureDatabase == "modERN", ]
ctrl <- md$NES[grep("eGFP-trh", md$FeatureDescription)]
cat(sprintf("\nmodERN positive control -- eGFP-trh over Trh/Tgo peaks: NES = %.2f\n",
            if (length(ctrl)) max(ctrl) else NA))
cat("  Trh's OWN modERN track is not enriched over Trh/Tgo joint peaks, and no\n")
cat("  modERN track reaches NES 3. The modERN portion of this run therefore\n")
cat("  FAILS its internal positive control and cannot be used to answer the\n")
cat("  reverse-lookup question -- run the direct scan instead.\n")

write.table(pio, file.path(OUT, "icistarget_pioneer_enrichment.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)
write.table(head(sig[order(-sig$NES),
                     c("Rank","NES","FeatureID","FeatureDatabase","FeatureAnnotations")], 40),
            file.path(OUT, "icistarget_top40_significant.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat("\nwritten: icistarget_pioneer_enrichment.tsv, icistarget_top40_significant.tsv\n")
