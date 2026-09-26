#!/usr/bin/env Rscript
# chipenrich_4x2.R, ChIP-Enrich gene-set enrichment.
# 4 heterodimers (Trh/Tgo, Sim/Tgo, Dys/Tgo, Sima/Tgo) x 2 custom tissue
# genesets (tracheal+SG vs complement; CNS+midline vs complement) = 8 runs.
#
# Usage:  Rscript --no-restore --no-save chipenrich_4x2.R
# SLURM:  see chipenrich_4x2.slurm

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(GenomeInfoDb)
  library(rtracklayer)
  library(chipenrich)
  library(AnnotationDbi)
  library(org.Dm.eg.db)
  library(dplyr)
  library(readr)
})

# ---- Paths -------------------------------------------------------------------
## Repository root: the environment variable PROJECT_DIR, or the working directory.
PROJECT_DIR <- normalizePath(Sys.getenv("PROJECT_DIR", unset = "."), mustWork = TRUE)
if (!file.exists(file.path(PROJECT_DIR, "README.md")))
  stop("PROJECT_DIR does not point at the repository root: ", PROJECT_DIR)
REF_DIR    <- file.path(PROJECT_DIR, "data/reference_gene_lists")
FINAL_DIR  <- file.path(PROJECT_DIR, "data/final_peaks")
EXPORT_DIR    <- file.path(PROJECT_DIR, "data/chapter4")
dir.create(EXPORT_DIR, showWarnings = FALSE, recursive = TRUE)

NCORES <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "1"))
cat(sprintf("cores=%d\n", NCORES))

# ---- Constants ---------------------------------------------------------------
good_chr  <- c("chr2L", "chr2R", "chr3L", "chr3R", "chr4", "chrX")
extraCols <- c(signalValue = "numeric", pValue = "numeric",
               qValue = "numeric", peak = "integer")

# ---- Load peak sets ----------------------------------------------------------
load_peaks <- function(path) {
  gr <- import(path, format = "BED", extraCols = extraCols)
  gr <- gr[seqnames(gr) %in% good_chr]
  gr <- keepSeqlevels(gr, good_chr, pruning.mode = "coarse")
  gr[!duplicated(gr)]
}

cat("Loading peaks...\n")
Tgo_clean  <- load_peaks(file.path(FINAL_DIR, "TGO_IDR0.05.final.narrowPeak"))
Sim_clean  <- load_peaks(file.path(FINAL_DIR, "SIM_IDR0.05.final.narrowPeak"))
Dys_clean  <- load_peaks(file.path(FINAL_DIR, "DYS_IDR0.05.final.narrowPeak"))
Sima_clean <- load_peaks(file.path(FINAL_DIR, "SIMA_IDR0.05.final.narrowPeak"))

# Trh/Tgo joint from Chapter 3
trh_bed <- file.path(PROJECT_DIR, "data/chapter3",
                     "Trh_Tgo_joint_peaks_native_width_for_MEME_ChIP.bed")
trh_df <- read.table(trh_bed, header = FALSE, stringsAsFactors = FALSE,
                     col.names = c("chr","start","end","name","score","strand"))
Trh_Tgo_joint_filt <- GRanges(
  seqnames = trh_df$chr,
  ranges   = IRanges(start = trh_df$start + 1L, end = trh_df$end)
)
Trh_Tgo_joint_filt <- keepSeqlevels(Trh_Tgo_joint_filt, good_chr,
                                    pruning.mode = "coarse")

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

# ---- Load gene lists ---------------------------------------------------------
read_fbgn <- function(path) {
  df <- read.csv(path, stringsAsFactors = FALSE)
  v  <- if ("genes" %in% names(df)) df[["genes"]] else df[[1]]
  v  <- unique(trimws(as.character(v)))
  v[grepl("^FBgn\\d+$", v)]
}
tracheal_sg_genes <- read_fbgn(file.path(REF_DIR,
                                "trh_ts_sg_tissues_reference.csv"))
cns_midline_genes <- read_fbgn(file.path(REF_DIR,
                                "sim_cns_tissues_DEGs.csv"))

# ---- Build 2-way geneset files (target vs complement) ------------------------
build_geneset_file <- function(fbgn_vec, label_target, label_comp, out_file) {
  m <- AnnotationDbi::select(org.Dm.eg.db,
                             keys    = fbgn_vec,
                             keytype = "FLYBASE",
                             columns = "ENTREZID")
  m <- dplyr::filter(m, !is.na(ENTREZID)) |>
       dplyr::distinct(ENTREZID)
  target_entrez <- as.character(m$ENTREZID)

  all_entrez <- unique(as.character(keys(org.Dm.eg.db, keytype = "ENTREZID")))
  target_entrez <- intersect(target_entrez, all_entrez)
  comp_entrez   <- setdiff(all_entrez, target_entrez)

  gs <- bind_rows(
    tibble(geneset = label_target, gene_id = target_entrez),
    tibble(geneset = label_comp,   gene_id = comp_entrez)
  )
  readr::write_delim(gs, out_file, delim = "\t", col_names = FALSE)
  list(n_target = length(target_entrez), n_comp = length(comp_entrez))
}

tracheal_gs_file <- file.path(EXPORT_DIR, "geneset_tracheal_vs_complement.tsv")
cns_gs_file      <- file.path(EXPORT_DIR, "geneset_cns_vs_complement.tsv")
stats_tracheal <- build_geneset_file(
  tracheal_sg_genes, "tracheal_sg", "all_dm6_minus_tracheal_sg",
  tracheal_gs_file)
stats_cns <- build_geneset_file(
  cns_midline_genes, "cns_midline", "all_dm6_minus_cns_midline",
  cns_gs_file)
cat(sprintf("Tracheal+SG: %d target / %d complement\n",
            stats_tracheal$n_target, stats_tracheal$n_comp))
cat(sprintf("CNS+midline: %d target / %d complement\n",
            stats_cns$n_target, stats_cns$n_comp))

# ---- chipenrich --------------------------------------------------------------
build_peaks_df <- function(gr) {
  data.frame(chr   = as.character(seqnames(gr)),
             start = start(gr),
             end   = end(gr),
             stringsAsFactors = FALSE)
}

run_chipenrich_custom <- function(joint_peaks, geneset_file, label) {
  res <- chipenrich(
    peaks            = build_peaks_df(joint_peaks),
    genome           = "dm6",
    locusdef         = "1kb",
    genesets         = geneset_file,
    method           = "chipenrich",
    min_geneset_size = 10,
    max_geneset_size = 15000,
    qc_plots         = FALSE,
    n_cores          = NCORES,
    out_name         = NULL
  )$results
  res$Pair <- label
  res
}

cat("Running chipenrich (8 calls)...\n")
ce_trh_ts   <- run_chipenrich_custom(Trh_Tgo_joint_filt,  tracheal_gs_file, "Trh/Tgo")
ce_trh_cns  <- run_chipenrich_custom(Trh_Tgo_joint_filt,  cns_gs_file,      "Trh/Tgo")
ce_sim_ts   <- run_chipenrich_custom(Sim_Tgo_joint_filt,  tracheal_gs_file, "Sim/Tgo")
ce_sim_cns  <- run_chipenrich_custom(Sim_Tgo_joint_filt,  cns_gs_file,      "Sim/Tgo")
ce_dys_ts   <- run_chipenrich_custom(Dys_Tgo_joint_filt,  tracheal_gs_file, "Dys/Tgo")
ce_dys_cns  <- run_chipenrich_custom(Dys_Tgo_joint_filt,  cns_gs_file,      "Dys/Tgo")
ce_sima_ts  <- run_chipenrich_custom(Sima_Tgo_joint_filt, tracheal_gs_file, "Sima/Tgo")
ce_sima_cns <- run_chipenrich_custom(Sima_Tgo_joint_filt, cns_gs_file,      "Sima/Tgo")

# ---- Compile -----------------------------------------------------------------
extract_target_row <- function(res, tissue_label) {
  tgt <- res[grepl("tracheal_sg|cns_midline", res$Description) &
             !grepl("minus", res$Description), , drop = FALSE]
  tgt$Tissue <- tissue_label
  tgt
}

ce_all <- rbind(
  extract_target_row(ce_trh_ts,   "Tracheal + SG"),
  extract_target_row(ce_trh_cns,  "CNS + Midline"),
  extract_target_row(ce_sim_ts,   "Tracheal + SG"),
  extract_target_row(ce_sim_cns,  "CNS + Midline"),
  extract_target_row(ce_dys_ts,   "Tracheal + SG"),
  extract_target_row(ce_dys_cns,  "CNS + Midline"),
  extract_target_row(ce_sima_ts,  "Tracheal + SG"),
  extract_target_row(ce_sima_cns, "CNS + Midline")
)

chipenrich_4x2 <- data.frame(
  Pair                 = ce_all$Pair,
  Tissue               = ce_all$Tissue,
  Odds_Ratio           = round(ce_all$Odds.Ratio, 3),
  P_value              = ce_all$P.value,
  FDR                  = ce_all$FDR,
  N_Geneset_Peak_Genes = ce_all$N.Geneset.Peak.Genes,
  stringsAsFactors = FALSE
)
print(chipenrich_4x2, row.names = FALSE)

# ---- Save --------------------------------------------------------------------
saveRDS(list(
  ce_raw = list(ce_trh_ts = ce_trh_ts, ce_trh_cns = ce_trh_cns,
                ce_sim_ts = ce_sim_ts, ce_sim_cns = ce_sim_cns,
                ce_dys_ts = ce_dys_ts, ce_dys_cns = ce_dys_cns,
                ce_sima_ts = ce_sima_ts, ce_sima_cns = ce_sima_cns),
  summary = chipenrich_4x2
), file.path(EXPORT_DIR, "chipenrich_4x2.rds"))
write.csv(chipenrich_4x2,
          file.path(EXPORT_DIR, "chipenrich_4x2_summary.csv"), row.names = FALSE)
cat(sprintf("Saved: %s/chipenrich_4x2.rds\n", EXPORT_DIR))
cat("Done.\n")
