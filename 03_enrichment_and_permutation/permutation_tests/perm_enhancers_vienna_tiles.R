#!/usr/bin/env Rscript
# ============================================================
# perm_enhancers_vienna_tiles.R
# Tests heterodimer peak overlap against Gao-Qian EnhancerAtlas (35798 intervals, genome-wide) and Vienna Tiles (Kvon 2014, 7477 tiles, 3457 active). Runs tile-conditional Fisher test (A) and genome-wide permutation (B).
# Input: data/reference_gene_lists/vt_enhancers_active_dm6.bed, vt_enhancers_active_stg9_16_dm6.bed, vt_tiles_tested_dm6.bed, merged_embryo_5-20_enhancers_dm6.bed (run prepare_vienna_tiles_dm6.R first)
# Output: data/chapter4/perm_enhancers_vt_summary.csv, perm_enhancers_vt.rds, vt_tile_conditional_fisher.csv; data/chapter3/perm_enhancers_vt_trh_tgo.csv
# Usage: Rscript --vanilla perm_enhancers_vienna_tiles.R [ntimes=5000]
# ============================================================

suppressPackageStartupMessages({
  library(GenomicRanges); library(GenomeInfoDb); library(rtracklayer)
  library(regioneR); library(TxDb.Dmelanogaster.UCSC.dm6.ensGene)
  library(BSgenome.Dmelanogaster.UCSC.dm6)
})

## Repository root: the environment variable PROJECT_DIR, or the working directory.
PROJECT_DIR <- normalizePath(Sys.getenv("PROJECT_DIR", unset = "."), mustWork = TRUE)
if (!file.exists(file.path(PROJECT_DIR, "README.md")))
  stop("PROJECT_DIR does not point at the repository root: ", PROJECT_DIR)

REF_DIR   <- file.path(PROJECT_DIR, "data/reference_gene_lists")
FINAL_DIR <- file.path(PROJECT_DIR, "data/final_peaks")
EXP3      <- file.path(PROJECT_DIR, "data/chapter3")
EXP4      <- file.path(PROJECT_DIR, "data/chapter4")
CKPT_DIR  <- file.path(PROJECT_DIR, "results/checkpoints/perm_vienna_tiles")
dir.create(CKPT_DIR, showWarnings = FALSE, recursive = TRUE)

args   <- commandArgs(trailingOnly = TRUE)
NTIMES <- if (length(args) >= 1 && !is.na(as.integer(args[1]))) as.integer(args[1]) else 5000L
NCORES <- suppressWarnings(as.integer(Sys.getenv("SLURM_CPUS_PER_TASK",
             unset = Sys.getenv("PERM_CORES", unset = NA))))
if (is.na(NCORES)) NCORES <- max(1L, parallel::detectCores() - 2L)
if (is.na(NCORES)) NCORES <- 6L
cat(sprintf("ntimes=%d | cores=%d\n", NTIMES, NCORES))

good_chr  <- c("chr2L", "chr2R", "chr3L", "chr3R", "chr4", "chrX")
txdb      <- TxDb.Dmelanogaster.UCSC.dm6.ensGene
extraCols <- c(signalValue = "numeric", pValue = "numeric",
               qValue = "numeric", peak = "integer")

## ---- canonical TSS points ---------------------------------------------------
tss_pt <- resize(transcripts(txdb), width = 1, fix = "start")
tss_pt <- keepSeqlevels(tss_pt, good_chr, pruning.mode = "coarse")

annotate_tss <- function(gr) {
  gr <- keepSeqlevels(gr, good_chr, pruning.mode = "coarse")
  ov <- distanceToNearest(gr, tss_pt, ignore.strand = TRUE)
  gr$distanceToTSS <- NA_integer_
  gr$distanceToTSS[queryHits(ov)] <- mcols(ov)$distance
  stopifnot(!anyNA(gr$distanceToTSS), all(gr$distanceToTSS >= 0))
  gr
}

## ---- peak sets (identical to perm_enhancer_overlap.R) --------------------
load_peaks <- function(path) {
  gr <- import(path, format = "BED", extraCols = extraCols)
  gr <- gr[seqnames(gr) %in% good_chr]
  gr <- keepSeqlevels(gr, good_chr, pruning.mode = "coarse")
  gr[!duplicated(gr)]
}
Tgo_clean  <- load_peaks(file.path(FINAL_DIR, "TGO_IDR0.05.final.narrowPeak"))
Sim_clean  <- load_peaks(file.path(FINAL_DIR, "SIM_IDR0.05.final.narrowPeak"))
Dys_clean  <- load_peaks(file.path(FINAL_DIR, "DYS_IDR0.05.final.narrowPeak"))
Sima_clean <- load_peaks(file.path(FINAL_DIR, "SIMA_IDR0.05.final.narrowPeak"))

trh_df <- read.table(file.path(EXP3, "Trh_Tgo_joint_peaks_native_width_for_MEME_ChIP.bed"),
                     header = FALSE, stringsAsFactors = FALSE,
                     col.names = c("chr","start","end","name","score","strand"))
Trh_Tgo_joint_filt <- GRanges(trh_df$chr,
                              IRanges(trh_df$start + 1L, trh_df$end),
                              name = trh_df$name)
Trh_Tgo_joint_filt <- keepSeqlevels(Trh_Tgo_joint_filt, good_chr, pruning.mode = "coarse")
stopifnot(length(Trh_Tgo_joint_filt) == 614L)

make_joint <- function(tf, tgo) {
  ov <- findOverlaps(tf, tgo, minoverlap = 1, ignore.strand = TRUE)
  gr <- pintersect(tf[queryHits(ov)], tgo[subjectHits(ov)])
  mcols(gr) <- NULL
  gr <- keepSeqlevels(gr, good_chr, pruning.mode = "coarse")
  gr[!duplicated(gr)]
}
Sim_Tgo_joint_filt  <- make_joint(Sim_clean,  Tgo_clean)
Dys_Tgo_joint_filt  <- make_joint(Dys_clean,  Tgo_clean)
Sima_Tgo_joint_filt <- make_joint(Sima_clean, Tgo_clean)

Trh_Tgo_joint_filt  <- annotate_tss(Trh_Tgo_joint_filt)
Sim_Tgo_joint_filt  <- annotate_tss(Sim_Tgo_joint_filt)
Dys_Tgo_joint_filt  <- annotate_tss(Dys_Tgo_joint_filt)
Sima_Tgo_joint_filt <- annotate_tss(Sima_Tgo_joint_filt)

## Check: canonical splits must match the values reported in Ch3/Ch4 text.
stopifnot(sum(abs(Trh_Tgo_joint_filt$distanceToTSS) <= 1000) == 468L,
          sum(abs(Trh_Tgo_joint_filt$distanceToTSS)  > 1000) == 146L,
          length(Sim_Tgo_joint_filt)  == 355L,
          length(Dys_Tgo_joint_filt)  == 870L,
          length(Sima_Tgo_joint_filt) == 296L)

HETS <- list("Trh/Tgo"  = Trh_Tgo_joint_filt,
             "Sim/Tgo"  = Sim_Tgo_joint_filt,
             "Dys/Tgo"  = Dys_Tgo_joint_filt,
             "Sima/Tgo" = Sima_Tgo_joint_filt)
HKEY <- c("Trh/Tgo" = "trh", "Sim/Tgo" = "sim",
          "Dys/Tgo" = "dys", "Sima/Tgo" = "sima")

## ---- reference interval sets ------------------------------------------------
read_ref <- function(f) {
  df <- read.table(file.path(REF_DIR, f), header = FALSE, stringsAsFactors = FALSE)
  colnames(df)[1:3] <- c("seqnames", "start", "end")
  gr <- GRanges(df$seqnames, IRanges(df$start + 1L, df$end))
  gr <- keepSeqlevels(gr, good_chr, pruning.mode = "coarse")
  gr[!duplicated(gr)]
}
REFS <- list(
  "GaoQian_EnhancerAtlas" = read_ref("merged_embryo_5-20_enhancers_dm6.bed"),
  "VT_active_any"         = read_ref("vt_enhancers_active_dm6.bed"),
  "VT_active_stg9_16"     = read_ref("vt_enhancers_active_stg9_16_dm6.bed"))
vt_tested <- read_ref("vt_tiles_tested_dm6.bed")

for (nm in names(REFS))
  cat(sprintf("reference %-24s : %6d intervals, median width %d bp\n",
              nm, length(REFS[[nm]]), median(width(REFS[[nm]]))))
cat(sprintf("VT tested universe       : %6d tiles,  %.1f%% of dm6 major arms\n",
            length(vt_tested),
            100 * sum(width(reduce(vt_tested))) /
              sum(seqlengths(keepSeqlevels(
                GRanges(seqinfo(BSgenome.Dmelanogaster.UCSC.dm6)),
                good_chr, pruning.mode = "coarse")))))
stopifnot(length(REFS$GaoQian_EnhancerAtlas) == 35798L)

## ============================================================
## (A) Tile-conditional Fisher test (primary VT statistic)
## ============================================================
## Tile activity flags, recovered by exact-coordinate match to the active sets.
key   <- function(gr) paste0(seqnames(gr), ":", start(gr), "-", end(gr))
k_any <- key(REFS$VT_active_any)
k_lat <- key(REFS$VT_active_stg9_16)
vt_tested$active_any  <- key(vt_tested) %in% k_any
vt_tested$active_late <- key(vt_tested) %in% k_lat
stopifnot(sum(vt_tested$active_any)  == length(REFS$VT_active_any),
          sum(vt_tested$active_late) == length(REFS$VT_active_stg9_16))

fisher_tiles <- function(joint, het_label, activity = c("any", "late")) {
  activity <- match.arg(activity)
  act <- if (activity == "any") vt_tested$active_any else vt_tested$active_late
  d   <- abs(joint$distanceToTSS)
  arms <- list(all      = joint,
               promoter = joint[d <= 1000],
               distal   = joint[d  > 1000])
  ## reference group held constant: tiles overlapping no peak of the full set
  unbound <- countOverlaps(vt_tested, joint, ignore.strand = TRUE) == 0L
  do.call(rbind, lapply(names(arms), function(a) {
    bound <- countOverlaps(vt_tested, arms[[a]], ignore.strand = TRUE) > 0L
    if (!any(bound)) return(NULL)
    tab <- matrix(c(sum(bound   & act), sum(bound   & !act),
                    sum(unbound & act), sum(unbound & !act)),
                  nrow = 2, byrow = TRUE,
                  dimnames = list(c("bound", "unbound"), c("active", "inactive")))
    ft <- fisher.test(tab)
    data.frame(Heterodimer = het_label, Activity = activity, Stratum = a,
               N_peaks = length(arms[[a]]),
               Tiles_bound = sum(bound), Bound_active = tab[1, 1],
               Bound_active_pct = round(100 * tab[1, 1] / sum(tab[1, ]), 1),
               Tiles_unbound = sum(unbound), Unbound_active = tab[2, 1],
               Unbound_active_pct = round(100 * tab[2, 1] / sum(tab[2, ]), 1),
               OR = round(unname(ft$estimate), 3),
               CI_lo = round(ft$conf.int[1], 3), CI_hi = round(ft$conf.int[2], 3),
               p = signif(ft$p.value, 3), stringsAsFactors = FALSE)
  }))
}

fisher_res <- do.call(rbind, c(
  lapply(names(HETS), function(h) fisher_tiles(HETS[[h]], h, "any")),
  lapply(names(HETS), function(h) fisher_tiles(HETS[[h]], h, "late"))))
write.csv(fisher_res, file.path(EXP4, "vt_tile_conditional_fisher.csv"),
          row.names = FALSE)
cat("\n=========== (A) VT TILE-CONDITIONAL FISHER ===========\n")
print(fisher_res, row.names = FALSE)

## ============================================================
## (B) GENOME-WIDE PERMUTATION, continuity with the Ch4 table
## ============================================================
genome_gr <- getGenomeAndMask("dm6", mask = NA)$genome
genome_gr <- keepSeqlevels(genome_gr, good_chr, pruning.mode = "coarse")

run_one_perm <- function(A, B, label, n = NTIMES, seed = 777) {
  set.seed(seed)
  t0 <- Sys.time()
  pt <- permTest(A = A, B = B, genome = genome_gr,
                 randomize.function = randomizeRegions,
                 evaluate.function  = numOverlaps,
                 ntimes = n, per.chromosome = TRUE, count.once = TRUE,
                 allow.overlaps = FALSE, mc.set.seed = FALSE,
                 mc.cores = NCORES, verbose = FALSE)
  cat(sprintf("  %-42s n=%4d obs=%4d  [%.1f min]\n", label, length(A),
              pt$numOverlaps$observed,
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  list(label = label, observed = pt$numOverlaps$observed,
       null = pt$numOverlaps$permuted, pval = pt$numOverlaps$pval, full = pt)
}

run_three_strata <- function(joint, het_label, ref, ref_label) {
  d <- abs(joint$distanceToTSS)
  arms <- list(all = joint, promoter = joint[d <= 1000], distal = joint[d > 1000])
  cat(sprintf("\n=== %s x %s ===  all=%d promoter=%d distal=%d\n",
              ref_label, het_label, length(arms$all),
              length(arms$promoter), length(arms$distal)))
  out <- lapply(names(arms), function(a)
    run_one_perm(arms[[a]], ref, sprintf("%s | %s %s", ref_label, het_label, a)))
  names(out) <- names(arms)
  out$n <- vapply(arms, length, integer(1))
  out
}

extract_perm <- function(pt, n, ref_label, het_label, stratum) {
  obs <- pt$observed; nv <- pt$null
  data.frame(Reference = ref_label, Heterodimer = het_label, Stratum = stratum,
             Context = pt$label, N_peaks = n, Observed = obs,
             Null_mean = round(mean(nv), 2), Null_SD = round(sd(nv), 2),
             Z_score = round((obs - mean(nv)) / sd(nv), 2),
             p_emp = pt$pval, N_perm = length(nv), stringsAsFactors = FALSE)
}

t_start <- Sys.time()
perm_store <- list(); rows <- list()
for (ref_label in names(REFS)) {
  for (het_label in names(HETS)) {
    ## ntimes is part of the key: a checkpoint from a low-ntimes smoke run
    ## must never be silently reused for the full run.
    f <- file.path(CKPT_DIR, sprintf("%s__%s__n%d.rds", ref_label,
                                     HKEY[[het_label]], NTIMES))
    if (file.exists(f)) {
      res <- readRDS(f)
      cat(sprintf("\n=== %s x %s === (cached)\n", ref_label, het_label))
    } else {
      res <- run_three_strata(HETS[[het_label]], het_label, REFS[[ref_label]], ref_label)
      saveRDS(res, f)
    }
    perm_store[[paste(ref_label, het_label)]] <- res
    for (a in c("all", "promoter", "distal"))
      rows[[length(rows) + 1L]] <- extract_perm(res[[a]], res$n[[a]],
                                                ref_label, het_label, a)
  }
}
perm_vt_results <- do.call(rbind, rows)
stopifnot(nrow(perm_vt_results) == 3L * length(REFS) * length(HETS))

saveRDS(list(perm = perm_store, summary = perm_vt_results, fisher = fisher_res),
        file.path(EXP4, "perm_enhancers_vt.rds"))
write.csv(perm_vt_results, file.path(EXP4, "perm_enhancers_vt_summary.csv"),
          row.names = FALSE)
write.csv(perm_vt_results[perm_vt_results$Heterodimer == "Trh/Tgo", ],
          file.path(EXP3, "perm_enhancers_vt_trh_tgo.csv"), row.names = FALSE)

cat("\n=========== (B) GENOME-WIDE PERMUTATION ===========\n")
print(perm_vt_results[, c("Reference","Heterodimer","Stratum","N_peaks",
                          "Observed","Null_mean","Null_SD","Z_score","p_emp")],
      row.names = FALSE)
cat(sprintf("\nTotal wall time: %.1f min\n",
            as.numeric(difftime(Sys.time(), t_start, units = "mins"))))
