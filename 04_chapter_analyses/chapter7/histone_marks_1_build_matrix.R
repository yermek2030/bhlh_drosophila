# ============================================================
# histone_marks_1_build_matrix.R
# Stage A: build dm6-anchored probe x (stage,mark) matrix from archived
#          ChIP-chip histone-mark bedGraphs; assign probes to 614 Trh/Tgo joint peaks.
# Archived bedGraphs are on dm2 coordinates (chr2h/3h/4h/Xh/Yh, GPL RANGE_GB AE014134.4).
# Lift dm2 -> dm3 -> dm6 before overlap. Writes data/chapter3/histone_marks_probe_peak.rds
# ============================================================

## Repository root: the environment variable PROJECT_DIR, or the working directory.
PROJECT_DIR <- normalizePath(Sys.getenv("PROJECT_DIR", unset = "."), mustWork = TRUE)
if (!file.exists(file.path(PROJECT_DIR, "README.md")))
  stop("PROJECT_DIR does not point at the repository root: ", PROJECT_DIR)
setwd(PROJECT_DIR)

suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
  library(rtracklayer)
  library(ChIPseeker)
  library(TxDb.Dmelanogaster.UCSC.dm6.ensGene)
})

MARK_DIR  <- file.path(PROJECT_DIR, "reference/chipchip_histone_marks")
CHAIN_DIR <- file.path(PROJECT_DIR, "reference", "liftover")
OUT_DIR   <- file.path(PROJECT_DIR, "data", "chapter3")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

STAGE_LEVELS <- c("E-4-8h", "E-8-12", "E-12-16h", "E-16-20h", "E-20-24h")
MARK_LEVELS  <- c("H3K4Me1", "H3K27Ac", "H3K4Me3", "H3K27Me3")

# ---------------------------------------------------------------- 1. bedGraphs
bg_files <- list.files(MARK_DIR, pattern = "\\.bedGraph$",
                       recursive = TRUE, full.names = TRUE)
bg_files <- bg_files[basename(dirname(bg_files)) != basename(MARK_DIR)]  # drop combined.bedGraph
stopifnot(length(bg_files) == 20L)

parse_tag <- function(p) {
  tag <- sub("\\.bedGraph$", "", basename(p))
  m   <- regmatches(tag, regexec("^(E-[0-9]+-[0-9]+h?)_(H3K[0-9]+[A-Za-z]+[0-9]*)$", tag))[[1]]
  stopifnot(length(m) == 3L)
  list(stage = m[2], mark = m[3])
}
tags <- lapply(bg_files, parse_tag)
stopifnot(all(vapply(tags, `[[`, "", "stage") %in% STAGE_LEVELS),
          all(vapply(tags, `[[`, "", "mark")  %in% MARK_LEVELS))

read_bg <- function(p) {
  dt <- fread(p, header = FALSE, colClasses = c("character","character","character","character"),
              showProgress = FALSE)
  setnames(dt, c("chrom", "start", "end", "score"))
  dt[, `:=`(start = suppressWarnings(as.numeric(start)),
            end   = suppressWarnings(as.numeric(end)),
            score = suppressWarnings(as.numeric(score)))]
  dt <- dt[!is.na(start) & !is.na(end) & !is.na(score) & end > start]
  dt[, key := paste0(chrom, ":", start, "-", end)]
  dt[!duplicated(key)]
}

message("reading 20 bedGraphs ...")
bg_list <- lapply(bg_files, read_bg)
names(bg_list) <- vapply(seq_along(tags), function(i)
  paste(tags[[i]]$stage, tags[[i]]$mark, sep = "|"), "")

# probes present in every file (identical probe design across stages/marks)
common <- Reduce(intersect, lapply(bg_list, `[[`, "key"))
message("common probes across all 20 files: ", length(common))
stopifnot(length(common) > 7e5)

ref <- bg_list[[1]][match(common, key)]
probe_dm2 <- GRanges(ref$chrom, IRanges(ref$start + 1L, ref$end))  # bedGraph is 0-based half-open
mcols(probe_dm2)$key <- ref$key

score_mat <- do.call(cbind, lapply(bg_list, function(dt) dt$score[match(common, dt$key)]))
colnames(score_mat) <- names(bg_list)
stopifnot(!anyNA(score_mat), nrow(score_mat) == length(probe_dm2))

# --------------------------------------------------- 2. liftOver dm2 -> dm6
ch_23 <- import.chain(file.path(CHAIN_DIR, "dm2ToDm3.over.chain"))
ch_36 <- import.chain(file.path(CHAIN_DIR, "dm3ToDm6.over.chain"))

lift_1to1 <- function(gr, ch) {
  L  <- liftOver(gr, ch)
  n  <- lengths(L)
  ok <- n == 1L
  out <- unlist(L[ok], use.names = FALSE)
  mcols(out)$key <- mcols(gr)$key[ok]
  out
}
message("lifting dm2 -> dm3 ...")
probe_dm3 <- lift_1to1(probe_dm2, ch_23)
message("lifting dm3 -> dm6 ...")
probe_dm6 <- lift_1to1(probe_dm3, ch_36)

# unique 1:1 keys only, sane widths, canonical arms
probe_dm6 <- probe_dm6[!(mcols(probe_dm6)$key %in%
                           mcols(probe_dm6)$key[duplicated(mcols(probe_dm6)$key)])]
probe_dm6 <- probe_dm6[width(probe_dm6) >= 30 & width(probe_dm6) <= 200]
good_chr  <- c("chr2L","chr2R","chr3L","chr3R","chr4","chrX")
probe_dm6 <- probe_dm6[as.character(seqnames(probe_dm6)) %in% good_chr]
seqlevels(probe_dm6) <- good_chr

lift_stats <- c(dm2 = length(probe_dm2), dm3 = length(probe_dm3), dm6 = length(probe_dm6))
message("liftOver retention: ", paste(names(lift_stats), lift_stats, sep = "=", collapse = "  "))

score_mat <- score_mat[match(mcols(probe_dm6)$key, common), , drop = FALSE]
stopifnot(nrow(score_mat) == length(probe_dm6))

# ---------------------------------- 3. within-file genome-wide percentile rank
# Three Agilent platforms were concatenated per file and the archived pipeline
# applied no cross-platform scaling, so raw M-values are not comparable across
# files. Percentile rank within each (stage,mark) file makes cross-stage
# differences interpretable; raw M is retained for sensitivity.
pct_mat <- apply(score_mat, 2, function(v) rank(v, ties.method = "average") / length(v))

# ------------------------------------------------------------ 4. the 614 peaks
peaks <- import(file.path(PROJECT_DIR, "data/chapter3/Trh_Tgo_joint_peaks_filt.bed"))
stopifnot(length(peaks) == 614L)
names(peaks) <- peaks$name

txdb <- TxDb.Dmelanogaster.UCSC.dm6.ensGene
# Distance convention: UNSIGNED distance from the peak interval to the
# nearest 1-bp transcript TSS (dm6 ensGene), ignore.strand = TRUE.  Reproduces
# joint_peaks_filt.rds$distanceToTSS exactly (614/614) -> 468 proximal / 146
# distal at 1 kb, matching Chapter 3 Fig 3.12, the threshold sweep and ChIP-Enrich.
good_chr <- c("chr2L", "chr2R", "chr3L", "chr3R", "chr4", "chrX")
tss_pt   <- resize(transcripts(txdb), width = 1, fix = "start")
tss_pt   <- keepSeqlevels(tss_pt, good_chr, pruning.mode = "coarse")
peaks    <- keepSeqlevels(peaks, good_chr, pruning.mode = "coarse")
stopifnot(length(peaks) == 614L)

ov_tss <- distanceToNearest(peaks, tss_pt, ignore.strand = TRUE)
mcols(peaks)$distanceToTSS <- NA_integer_
mcols(peaks)$distanceToTSS[queryHits(ov_tss)] <- mcols(ov_tss)$distance
stopifnot(!anyNA(mcols(peaks)$distanceToTSS),
          all(mcols(peaks)$distanceToTSS >= 0))

fimo <- fread(file.path(PROJECT_DIR, "data/chapter3/fimo_joint/fimo.tsv"),
              header = TRUE, showProgress = FALSE)
fimo <- fimo[motif_id == "STRCACTGARARAAA"]          # MANDATORY: fimo.tsv is multi-motif
m2_peaks <- unique(fimo$sequence_name)
message("Motif2-positive peaks: ", length(m2_peaks), " / 614")
mcols(peaks)$motif2 <- peaks$name %in% m2_peaks
mcols(peaks)$distal <- abs(mcols(peaks)$distanceToTSS) > 1000
stopifnot(sum(mcols(peaks)$distal) == 146L)   # Ch3/Ch4 canonical: 468 / 146
# On the signed annotatePeak distance, `> 1000` yields 90 distal peaks because
# 56 far-upstream peaks carry negative distances and fall into the promoter
# class. The unsigned distance yields the 146 distal peaks used by Ch3 Fig
# 3.12, the threshold sweep and ChIP-Enrich.

# --------------------------------------------- 5. probe -> peak assignment
hits <- findOverlaps(peaks, probe_dm6)
message("peaks with >=1 probe: ", length(unique(queryHits(hits))), " / 614")

saveRDS(list(
  probe_dm6   = probe_dm6,
  score_mat   = score_mat,
  pct_mat     = pct_mat,
  peaks       = peaks,
  hits        = hits,
  lift_stats  = lift_stats,
  n_common    = length(common),
  stage_levels = STAGE_LEVELS,
  mark_levels  = MARK_LEVELS,
  built_at    = Sys.time()
), file.path(OUT_DIR, "histone_marks_probe_peak.rds"))

message("wrote ", file.path(OUT_DIR, "histone_marks_probe_peak.rds"))
