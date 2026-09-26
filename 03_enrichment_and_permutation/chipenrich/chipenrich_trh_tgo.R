#!/usr/bin/env Rscript
# chipenrich_trh_tgo.R, HPC ChIP-Enrich, consolidated:
#   1. res_custom       (custom tracheal+SG vs complement, method=chipenrich)
#   2. res_custom_fet   (same gene-set, method=fet)
#   3. res_BP           (GOBP, method=chipenrich)
#   4. res_MF           (GOMF, method=chipenrich)
#   5. res_CC           (GOCC, method=chipenrich)
#   6. res_reactome     (Reactome, method=chipenrich)

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(GenomeInfoDb)
  library(rtracklayer)
  library(chipenrich)
  library(AnnotationDbi)
  library(org.Dm.eg.db)
  library(dplyr)
  library(readr)
  library(tibble)
})

# ---- Paths -------------------------------------------------------------------
## Repository root: the environment variable PROJECT_DIR, or the working directory.
PROJECT_DIR <- normalizePath(Sys.getenv("PROJECT_DIR", unset = "."), mustWork = TRUE)
if (!file.exists(file.path(PROJECT_DIR, "README.md")))
  stop("PROJECT_DIR does not point at the repository root: ", PROJECT_DIR)
REF_DIR   <- file.path(PROJECT_DIR, "data/reference_gene_lists")
EXPORT_DIR   <- file.path(PROJECT_DIR, "data/chapter3")
dir.create(EXPORT_DIR, showWarnings = FALSE, recursive = TRUE)

NCORES <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "1"))
cat(sprintf("cores=%d\n", NCORES))

good_chr <- c("chr2L", "chr2R", "chr3L", "chr3R", "chr4", "chrX")

# ---- joint_peaks_filt -> peaks_df -------------------------------------------
trh_bed <- file.path(PROJECT_DIR, "data/chapter3",
                     "Trh_Tgo_joint_peaks_native_width_for_MEME_ChIP.bed")
stopifnot(file.exists(trh_bed))
trh_df <- read.table(trh_bed, header = FALSE, stringsAsFactors = FALSE,
                     col.names = c("chr","start","end","name","score","strand"))
joint_peaks_filt <- GRanges(
  seqnames = trh_df$chr,
  ranges   = IRanges(start = trh_df$start + 1L, end = trh_df$end)
)
joint_peaks_filt <- keepSeqlevels(joint_peaks_filt, good_chr,
                                  pruning.mode = "coarse")
stopifnot(length(joint_peaks_filt) == 614)

peaks_df <- data.frame(
  chr   = as.character(seqnames(joint_peaks_filt)),
  start = start(joint_peaks_filt),
  end   = end(joint_peaks_filt),
  stringsAsFactors = FALSE
)
cat(sprintf("peaks_df: %d rows\n", nrow(peaks_df)))

# ---- Build the custom tracheal-vs-complement gene-set file ------------------
read_fbgn <- function(path) {
  df <- read.csv(path, stringsAsFactors = FALSE)
  v  <- if ("genes" %in% names(df)) df[["genes"]] else df[[1]]
  v  <- unique(trimws(as.character(v)))
  v[grepl("^FBgn\\d+$", v)]
}
map_fbgn_to_entrez <- function(fbgn_vec) {
  if (length(fbgn_vec) == 0) return(character(0))
  m <- AnnotationDbi::select(org.Dm.eg.db, keys = fbgn_vec,
                             keytype = "FLYBASE", columns = "ENTREZID")
  m <- dplyr::filter(m, !is.na(ENTREZID)) |> dplyr::distinct(ENTREZID)
  as.character(m$ENTREZID)
}

all_dm6      <- unique(as.character(keys(org.Dm.eg.db, keytype = "ENTREZID")))
trach_fbgn   <- read_fbgn(file.path(REF_DIR, "trh_ts_sg_tissues_reference.csv"))
trach_entrez <- intersect(map_fbgn_to_entrez(trach_fbgn), all_dm6)
comp_entrez  <- setdiff(all_dm6, trach_entrez)
stopifnot(length(intersect(trach_entrez, comp_entrez)) == 0)

cat(sprintf("Tracheal geneset:   %d Entrez IDs\n", length(trach_entrez)))
cat(sprintf("Complement geneset: %d Entrez IDs\n", length(comp_entrez)))

genesets_file <- file.path(EXPORT_DIR, "genesets_tracheal_vs_complement.tsv")
custom_gs <- bind_rows(
  tibble(geneset = "tracheal_custom",        gene_id = trach_entrez),
  tibble(geneset = "All_dm6_minus_tracheal", gene_id = comp_entrez)
)
readr::write_delim(custom_gs, genesets_file, delim = "\t", col_names = FALSE)

# ---- (1) PRIMARY custom, method = chipenrich (logistic regression) --------
cat("\n[1/6] res_custom (custom, method=chipenrich)...\n")
res_custom <- chipenrich(
  peaks            = peaks_df,
  genome           = "dm6",
  locusdef         = "1kb",
  genesets         = genesets_file,
  method           = "chipenrich",
  min_geneset_size = 10,
  max_geneset_size = 15000,
  qc_plots         = FALSE,
  n_cores          = NCORES,
  out_name         = NULL
)$results

# ---- (2) Sensitivity custom, method = fet ---------------------------------
cat("[2/6] res_custom_fet (custom, method=fet)...\n")
res_custom_fet <- chipenrich(
  peaks            = peaks_df,
  genome           = "dm6",
  locusdef         = "1kb",
  genesets         = genesets_file,
  method           = "fet",
  min_geneset_size = 10,
  max_geneset_size = 15000,
  qc_plots         = FALSE,
  n_cores          = NCORES,
  out_name         = NULL
)$results

# ---- (3)-(6) Built-in ontologies (chipenrich method) ------------------------
run_builtin <- function(gs_name, label, idx) {
  cat(sprintf("[%d/6] res_%s (%s, method=chipenrich)...\n",
              idx, gs_name, gs_name))
  chipenrich(
    peaks            = peaks_df,
    genome           = "dm6",
    locusdef         = "1kb",
    genesets         = gs_name,
    method           = "chipenrich",
    qc_plots         = FALSE,
    n_cores          = NCORES,
    out_name         = NULL
  )$results
}

res_BP       <- run_builtin("GOBP",     "BP",       3)
res_MF       <- run_builtin("GOMF",     "MF",       4)
res_CC       <- run_builtin("GOCC",     "CC",       5)
res_reactome <- run_builtin("reactome", "reactome", 6)

# Significant subsets (FDR < 0.05), used by downstream chunks
res_sig_BP       <- res_BP[res_BP$FDR             < 0.05, ]
res_sig_MF       <- res_MF[res_MF$FDR             < 0.05, ]
res_sig_CC       <- res_CC[res_CC$FDR             < 0.05, ]
res_sig_reactome <- res_reactome[res_reactome$FDR < 0.05, ]

cat(sprintf("\nGOBP: %d total, %d FDR<0.05\n",
            nrow(res_BP), nrow(res_sig_BP)))
cat(sprintf("GOMF: %d total, %d FDR<0.05\n",
            nrow(res_MF), nrow(res_sig_MF)))
cat(sprintf("GOCC: %d total, %d FDR<0.05\n",
            nrow(res_CC), nrow(res_sig_CC)))
cat(sprintf("Reactome: %d total, %d FDR<0.05\n",
            nrow(res_reactome), nrow(res_sig_reactome)))

# ---- Save -------------------------------------------------------------------
saveRDS(list(
  res_custom       = res_custom,
  res_custom_fet   = res_custom_fet,
  res_BP           = res_BP,
  res_MF           = res_MF,
  res_CC           = res_CC,
  res_reactome     = res_reactome,
  res_sig_BP       = res_sig_BP,
  res_sig_MF       = res_sig_MF,
  res_sig_CC       = res_sig_CC,
  res_sig_reactome = res_sig_reactome,
  trach_entrez_n   = length(trach_entrez),
  comp_entrez_n    = length(comp_entrez)
), file.path(EXPORT_DIR, "chipenrich_ch3.rds"))

# Concise CSVs of the primary tables
write.csv(res_custom,
          file.path(EXPORT_DIR, "chipenrich_ch3_custom.csv"),     row.names = FALSE)
write.csv(res_custom_fet,
          file.path(EXPORT_DIR, "chipenrich_ch3_custom_fet.csv"), row.names = FALSE)
write.csv(res_sig_BP,
          file.path(EXPORT_DIR, "chipenrich_ch3_GOBP_sig.csv"),   row.names = FALSE)

cat(sprintf("\nSaved: %s/chipenrich_ch3.rds\n", EXPORT_DIR))
cat("Done.\n")
