#!/usr/bin/env Rscript
# ============================================================
# prepare_vienna_tiles_dm6.R: builds VT in vivo reporter enhancer set (7793 tiles), dm3 to dm6 liftover, as comparator to EnhancerAtlas set
# input: data/vienna_tiles/vt_tiles_all_downloaded.csv, reference/liftover/dm3ToDm6.over.chain
# output (data/reference_gene_lists/): vt_tiles_tested_dm6.bed, vt_enhancers_active_dm6.bed, vt_enhancers_active_stg9_16_dm6.bed, vt_

suppressPackageStartupMessages({
  library(GenomicRanges); library(GenomeInfoDb); library(rtracklayer)
})

## Repository root: the environment variable PROJECT_DIR, or the working directory.
PROJECT_DIR <- normalizePath(Sys.getenv("PROJECT_DIR", unset = "."), mustWork = TRUE)
if (!file.exists(file.path(PROJECT_DIR, "README.md")))
  stop("PROJECT_DIR does not point at the repository root: ", PROJECT_DIR)

VT_DIR  <- file.path(PROJECT_DIR, "data", "vienna_tiles")
REF_DIR <- file.path(PROJECT_DIR, "data/reference_gene_lists")
CHAIN   <- file.path(PROJECT_DIR, "reference", "liftover", "dm3ToDm6.over.chain")
stopifnot(dir.exists(VT_DIR), dir.exists(REF_DIR), file.exists(CHAIN))

good_chr    <- c("chr2L", "chr2R", "chr3L", "chr3R", "chr4", "chrX")
STAGE_COLS  <- c("stg4_6", "stg7_8", "stg9_10", "stg11_12", "stg13_14", "stg15_16")
LATE_COLS   <- c("stg9_10", "stg11_12", "stg13_14", "stg15_16")

## ---- read -------------------------------------------------------------------
vt <- read.csv(file.path(VT_DIR, "vt_tiles_all_downloaded.csv"),
               stringsAsFactors = FALSE, check.names = TRUE,
               na.strings = c("NA", ""))
## server header is: VTID, Chrosome [sic], Start, End, Length,
## Verification Status, Positive, stg4_6 ... stg15_16
names(vt)[names(vt) == "Chrosome"]            <- "chrom"
names(vt)[names(vt) == "Verification.Status"] <- "verification"
stopifnot(all(c("VTID", "chrom", "Start", "End", "Length",
                "verification", "Positive", STAGE_COLS) %in% names(vt)))
stopifnot(nrow(vt) == 7793L, !anyDuplicated(vt$VTID))
stopifnot(all(vt$End - vt$Start == vt$Length))   # 0-based half-open confirmed

cat(sprintf("VT download        : %d tiles\n", nrow(vt)))
cat(sprintf("  verification=correct : %d\n", sum(vt$verification == "correct")))
cat(sprintf("  Positive == 1        : %d\n", sum(vt$Positive == 1L)))

## Stage cell holds "tissue;strength" tokens joined by "|"; NA = no activity; strength scale 1-5.
## Positive requires max strength >= 2 across stages 4-16 (n=3
has_activity <- function(x) !is.na(x) & nzchar(x)
stage_any  <- Reduce(`|`, lapply(vt[STAGE_COLS], has_activity))
late_ann   <- Reduce(`|`, lapply(vt[LATE_COLS],  has_activity))

max_strength <- function(row) {
  tk <- unlist(strsplit(na.omit(unlist(row)), "\\|"))
  if (!length(tk)) return(NA_integer_)
  suppressWarnings(max(as.integer(sub("^.*;", "", tk)), na.rm = TRUE))
}
ms <- apply(vt[STAGE_COLS], 1, max_strength)
## Check the partition above; if a future release changes the convention this
## fails loudly instead of quietly redefining the reference set.
stopifnot(all(is.na(ms[vt$Positive == 1L]) == FALSE),
          min(ms[vt$Positive == 1L]) >= 2L,
          all(ms[vt$Positive == 0L & stage_any] == 1L))

vt$active_any  <- vt$Positive == 1L
vt$active_late <- vt$Positive == 1L & late_ann

cat(sprintf("  strength-1 tiles scored Positive==0 : %d\n",
            sum(vt$Positive == 0L & stage_any)))
cat(sprintf("  Positive==1 & annotated stg9_10+     : %d\n", sum(vt$active_late)))

## ---- verified tiles only ----------------------------------------------------
vt_ok <- vt[vt$verification == "correct", ]
cat(sprintf("after verification filter: %d tiles (%d dropped)\n",
            nrow(vt_ok), nrow(vt) - nrow(vt_ok)))

## ---- dm3 GRanges ------------------------------------------------------------
gr3 <- GRanges(seqnames = vt_ok$chrom,
               ranges   = IRanges(vt_ok$Start + 1L, vt_ok$End),
               VTID       = vt_ok$VTID,
               active_any = vt_ok$active_any,
               active_late= vt_ok$active_late)
cat("dm3 chromosome tally:\n"); print(table(as.character(seqnames(gr3))))

## ---- liftOver dm3 -> dm6 ----------------------------------------------------
ch  <- import.chain(CHAIN)
lst <- liftOver(gr3, ch)
names(lst) <- gr3$VTID

n_seg    <- lengths(lst)
n_failed <- sum(n_seg == 0L)
n_split  <- sum(n_seg > 1L)
cat(sprintf("liftOver: %d tiles unmapped, %d split into >1 segment\n",
            n_failed, n_split))

gr6 <- unlist(lst, use.names = FALSE)
seg_vtid <- rep(gr3$VTID,        times = n_seg)
seg_any  <- rep(gr3$active_any,  times = n_seg)
seg_late <- rep(gr3$active_late, times = n_seg)
stopifnot(length(gr6) == length(seg_vtid))

## widest segment per tile -> one dm6 interval per tested fragment
ord  <- order(seg_vtid, -width(gr6))
keep <- ord[!duplicated(seg_vtid[ord])]
gr6  <- gr6[keep]
gr6$VTID        <- seg_vtid[keep]
gr6$active_any  <- seg_any[keep]
gr6$active_late <- seg_late[keep]
stopifnot(!anyDuplicated(gr6$VTID))

## ---- restrict to the six major arms ----------------------------------------
gr6 <- gr6[as.character(seqnames(gr6)) %in% good_chr]
gr6 <- keepSeqlevels(gr6, good_chr, pruning.mode = "coarse")
gr6 <- sort(gr6)

n_tested <- length(gr6)
n_any    <- sum(gr6$active_any)
n_late   <- sum(gr6$active_late)
cat(sprintf("\ndm6, six major arms:\n  tiles tested ......... %d\n", n_tested))
cat(sprintf("  active any stage ..... %d (%.1f%%)\n", n_any,  100 * n_any  / n_tested))
cat(sprintf("  active stg9_10-15_16 . %d (%.1f%%)\n", n_late, 100 * n_late / n_tested))
cat(sprintf("  median width ......... %d bp\n", median(width(gr6))))

## width sanity: tiles are ~2 kb by design; a lift that mangled them would show
stopifnot(median(width(gr6)) > 1500, median(width(gr6)) < 2600)

## ---- write ------------------------------------------------------------------
write_bed <- function(gr, path) {
  df <- data.frame(chrom = as.character(seqnames(gr)),
                   start = start(gr) - 1L,      # back to 0-based half-open
                   end   = end(gr),
                   name  = gr$VTID,
                   stringsAsFactors = FALSE)
  write.table(df, path, sep = "\t", quote = FALSE,
              row.names = FALSE, col.names = FALSE)
  cat(sprintf("wrote %s (%d intervals)\n", basename(path), nrow(df)))
}
write_bed(gr6,                    file.path(REF_DIR, "vt_tiles_tested_dm6.bed"))
write_bed(gr6[gr6$active_any],    file.path(REF_DIR, "vt_enhancers_active_dm6.bed"))
write_bed(gr6[gr6$active_late],   file.path(REF_DIR, "vt_enhancers_active_stg9_16_dm6.bed"))

ann <- data.frame(VTID = gr6$VTID,
                  chrom = as.character(seqnames(gr6)),
                  start = start(gr6) - 1L, end = end(gr6),
                  active_any = as.integer(gr6$active_any),
                  active_late = as.integer(gr6$active_late),
                  stringsAsFactors = FALSE)
ann <- merge(ann, vt_ok[, c("VTID", STAGE_COLS)], by = "VTID", all.x = TRUE)
write.table(ann, file.path(REF_DIR, "vt_tile_annotation.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)
cat(sprintf("wrote vt_tile_annotation.tsv (%d rows)\n", nrow(ann)))

cat("\ndone.\n")
