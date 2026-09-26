#!/usr/bin/env Rscript
# ------------------------------------------------------------
# perm_enhancer_overlap.R, HPC driver; reads PROJECT_DIR and SLURM_CPUS_PER_TASK; checkpoints per heterodimer under results/checkpoints/perm_enhancers/, resumable.
# Runs Chapter 3 (3 strata) and Chapter 4 (12 strata) enhancer-overlap permutation tests using the proximal/distal convention

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

NTIMES <- 5000L
NCORES <- suppressWarnings(as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = Sys.getenv("PERM_CORES", unset = NA))))
if (is.na(NCORES)) NCORES <- max(1L, parallel::detectCores() - 2L)
if (is.na(NCORES)) NCORES <- 6L
cat(sprintf("ntimes=%d | cores=%d\n", NTIMES, NCORES))

good_chr  <- c("chr2L", "chr2R", "chr3L", "chr3R", "chr4", "chrX")
txdb      <- TxDb.Dmelanogaster.UCSC.dm6.ensGene
extraCols <- c(signalValue = "numeric", pValue = "numeric",
               qValue = "numeric", peak = "integer")

## ---- canonical TSS points -------------------------------------------------
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

## ---- peak sets ------------------------------------------------------------
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

## ---- enhancers + genome ---------------------------------------------------
enh_df <- read.table(file.path(REF_DIR, "merged_embryo_5-20_enhancers_dm6.bed"),
                     header = FALSE, stringsAsFactors = FALSE)
colnames(enh_df)[1:3] <- c("seqnames","start","end")
enhancers_filt <- GRanges(enh_df$seqnames, IRanges(enh_df$start + 1L, enh_df$end))
enhancers_filt <- keepSeqlevels(enhancers_filt, good_chr, pruning.mode = "coarse")
cat(sprintf("enhancers_filt: %d\n", length(enhancers_filt)))

genome_gr <- getGenomeAndMask("dm6", mask = NA)$genome
genome_gr <- keepSeqlevels(genome_gr, good_chr, pruning.mode = "coarse")

## ---- permutation helpers --------------------------------------------------
run_one_perm <- function(A, B, label, n = NTIMES, seed = 777) {
  set.seed(seed)
  t0 <- Sys.time()
  pt <- permTest(A = A, B = B, genome = genome_gr,
                 randomize.function = randomizeRegions,
                 evaluate.function  = numOverlaps,
                 ntimes = n, per.chromosome = TRUE, count.once = TRUE,
                 allow.overlaps = FALSE, mc.set.seed = FALSE,
                 mc.cores = NCORES, verbose = FALSE)
  cat(sprintf("  %-24s n=%4d obs=%4d  [%.1f min]\n", label, length(A),
              pt$numOverlaps$observed,
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  list(label = label, observed = pt$numOverlaps$observed,
       null = pt$numOverlaps$permuted, pval = pt$numOverlaps$pval,
       full = pt)
}

run_three_strata <- function(joint_peaks, het_label) {
  d <- abs(joint_peaks$distanceToTSS)
  promoter_peaks <- joint_peaks[d <= 1000]
  distal_peaks   <- joint_peaks[d  > 1000]
  cat(sprintf("\n=== %s ===  all=%d promoter=%d distal=%d\n", het_label,
              length(joint_peaks), length(promoter_peaks), length(distal_peaks)))
  list(n_all = length(joint_peaks),
       n_promoter = length(promoter_peaks), n_distal = length(distal_peaks),
       all      = run_one_perm(joint_peaks,    enhancers_filt, paste(het_label, "all")),
       promoter = run_one_perm(promoter_peaks, enhancers_filt, paste(het_label, "promoter")),
       distal   = run_one_perm(distal_peaks,   enhancers_filt, paste(het_label, "distal")))
}

extract_perm <- function(pt, n) {
  obs <- pt$observed; nv <- pt$null
  data.frame(Context = pt$label, N_peaks = n, Observed = obs,
             Null_mean = round(mean(nv), 2), Null_SD = round(sd(nv), 2),
             Z_score = round((obs - mean(nv)) / sd(nv), 2),
             p_emp = pt$pval, N_perm = length(nv), stringsAsFactors = FALSE)
}

## ---- run ------------------------------------------------------------------
## Each heterodimer is checkpointed to CKPT_DIR immediately after its three
## strata finish, so an interrupted run resumes instead of restarting.
CKPT_DIR <- file.path(PROJECT_DIR, "results/checkpoints/perm_enhancers")
dir.create(CKPT_DIR, showWarnings = FALSE, recursive = TRUE)

cached_strata <- function(gr, het_label, key) {
  f <- file.path(CKPT_DIR, paste0(key, ".rds"))
  if (file.exists(f)) { cat(sprintf("\n=== %s === (cached)\n", het_label)); return(readRDS(f)) }
  res <- run_three_strata(gr, het_label)
  saveRDS(res, f)
  res
}

t_start <- Sys.time()
pt_enh_trh  <- cached_strata(Trh_Tgo_joint_filt,  "Trh/Tgo",  "trh")
pt_enh_sim  <- cached_strata(Sim_Tgo_joint_filt,  "Sim/Tgo",  "sim")
pt_enh_dys  <- cached_strata(Dys_Tgo_joint_filt,  "Dys/Tgo",  "dys")
pt_enh_sima <- cached_strata(Sima_Tgo_joint_filt, "Sima/Tgo", "sima")

perm_enh_results <- do.call(rbind, list(
  extract_perm(pt_enh_trh$all,       pt_enh_trh$n_all),
  extract_perm(pt_enh_trh$promoter,  pt_enh_trh$n_promoter),
  extract_perm(pt_enh_trh$distal,    pt_enh_trh$n_distal),
  extract_perm(pt_enh_sim$all,       pt_enh_sim$n_all),
  extract_perm(pt_enh_sim$promoter,  pt_enh_sim$n_promoter),
  extract_perm(pt_enh_sim$distal,    pt_enh_sim$n_distal),
  extract_perm(pt_enh_dys$all,       pt_enh_dys$n_all),
  extract_perm(pt_enh_dys$promoter,  pt_enh_dys$n_promoter),
  extract_perm(pt_enh_dys$distal,    pt_enh_dys$n_distal),
  extract_perm(pt_enh_sima$all,      pt_enh_sima$n_all),
  extract_perm(pt_enh_sima$promoter, pt_enh_sima$n_promoter),
  extract_perm(pt_enh_sima$distal,   pt_enh_sima$n_distal)))

## ---- Ch4 output (12 rows) -------------------------------------------------
stopifnot(nrow(perm_enh_results) == 12L)
saveRDS(list(pt_enh_trh = pt_enh_trh, pt_enh_sim = pt_enh_sim,
             pt_enh_dys = pt_enh_dys, pt_enh_sima = pt_enh_sima,
             summary = perm_enh_results),
        file.path(EXP4, "perm_enhancers.rds"))
write.csv(perm_enh_results, file.path(EXP4, "perm_enhancers_summary.csv"),
          row.names = FALSE)

## ---- Ch3 output (3 rows, Trh/Tgo only) ------------------------------------
ch3_summary <- data.frame(
  Context   = c("All Trh/Tgo peaks", "Promoter-proximal (<=1 kb)", "Distal (>1 kb)"),
  Observed  = c(pt_enh_trh$all$observed, pt_enh_trh$promoter$observed,
                pt_enh_trh$distal$observed),
  Null_mean = round(sapply(list(pt_enh_trh$all, pt_enh_trh$promoter,
                                pt_enh_trh$distal), function(z) mean(z$null)), 2),
  Null_SD   = round(sapply(list(pt_enh_trh$all, pt_enh_trh$promoter,
                                pt_enh_trh$distal), function(z) sd(z$null)), 2),
  Z_score   = round(sapply(list(pt_enh_trh$all, pt_enh_trh$promoter,
                                pt_enh_trh$distal),
                           function(z) (z$observed - mean(z$null)) / sd(z$null)), 2),
  p_emp     = sapply(list(pt_enh_trh$all, pt_enh_trh$promoter,
                          pt_enh_trh$distal), function(z) z$pval),
  N_perm    = NTIMES, stringsAsFactors = FALSE)

saveRDS(list(pt_all = pt_enh_trh$all$full, pt_promoter = pt_enh_trh$promoter$full,
             pt_distal = pt_enh_trh$distal$full, summary = ch3_summary),
        file.path(EXP3, "perm_enhancers_trh_tgo.rds"))
write.csv(ch3_summary, file.path(EXP3, "perm_enhancers_trh_tgo_summary.csv"),
          row.names = FALSE)

cat("\n================ CHAPTER 3 ================\n"); print(ch3_summary, row.names = FALSE)
cat("\n================ CHAPTER 4 ================\n"); print(perm_enh_results, row.names = FALSE)
cat(sprintf("\nTotal wall time: %.1f min\n",
            as.numeric(difftime(Sys.time(), t_start, units = "mins"))))
