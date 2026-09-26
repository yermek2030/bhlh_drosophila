# ============================================================
# ch4_specificity_test.R: interaction statistic for the 4x2 design
#
# tests whether TSG-vs-CNS preference differs between heterodimers, using
# disjoint gene lists and a within-row contrast. Adds peak-number-matched
# ranking and control-TF row. Run after tss_gr, all_dm6_fbgn, the four
# *_joint_filt objects and tracheal_sg_genes/cns_midline_genes; otherwise rebuilds peaks from data/final_peaks/.
# ============================================================

## The figure fig4_9_S_specificity_reanalysis.png is drawn by
## chapter3_4_figures.Rmd from the tables written here.

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
  library(GenomicRanges); library(GenomeInfoDb); library(rtracklayer)
})

good_chr <- c("chr2L","chr2R","chr3L","chr3R","chr4","chrX")
PK_DIR   <- file.path(PROJECT_DIR, "data/final_peaks")

# ---- 0. peak sets (use in-session objects if present, else rebuild) ---------
.load_peaks <- function(tag) {
  f <- file.path(PK_DIR, sprintf("%s_IDR0.05.final.narrowPeak", tag))
  d <- read.table(f, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
  gr <- GRanges(d[[1]], IRanges(d[[2]] + 1L, d[[3]]))
  gr <- gr[seqnames(gr) %in% good_chr]
  gr <- keepSeqlevels(gr, good_chr, pruning.mode = "coarse")
  gr[!duplicated(gr)]
}
.joint <- function(a_tag) {
  a <- .load_peaks(a_tag); b <- .load_peaks("TGO")
  o <- findOverlaps(a, b, minoverlap = 1L, ignore.strand = TRUE)
  j <- pintersect(a[queryHits(o)], b[subjectHits(o)]); mcols(j) <- NULL
  j <- keepSeqlevels(j, good_chr, pruning.mode = "coarse")
  j[!duplicated(j)]
}
J <- list(
  `Trh/Tgo`  = if (exists("Trh_Tgo_joint_filt"))  Trh_Tgo_joint_filt  else .joint("TRH"),
  `Sim/Tgo`  = if (exists("Sim_Tgo_joint_filt"))  Sim_Tgo_joint_filt  else .joint("SIM"),
  `Dys/Tgo`  = if (exists("Dys_Tgo_joint_filt"))  Dys_Tgo_joint_filt  else .joint("DYS"),
  `Sima/Tgo` = if (exists("Sima_Tgo_joint_filt")) Sima_Tgo_joint_filt else .joint("SIMA"))
stopifnot(length(J[["Trh/Tgo"]]) == 614L, length(J[["Sim/Tgo"]])  == 355L,
          length(J[["Dys/Tgo"]]) == 870L, length(J[["Sima/Tgo"]]) == 296L)
## List names are title-case (Twi/Vvl/Rib) for y-axis labels; .load_peaks() keys stay uppercase to match filenames.
## Da (ENCODE ENCSR416MJQ, IDR 0.05, 923 peaks) is the control for a different obligate partner.
## Pipeline: peak -> nearest TSS -> tissue list; no FIMO scan, no GC/length-matched background; IDR peak set is the only input.
CTRL <- list(Twi = .load_peaks("TWI"), Vvl = .load_peaks("VVL"),
             Rib = .load_peaks("RIB"), Da  = .load_peaks("DA"))
SETS <- c(J, CTRL)
hets <- names(J)

# ---- 0b. annotation + tissue lists (rebuilt only if not already in session) --
if (!exists("tss_gr")) {
  suppressPackageStartupMessages({
    library(TxDb.Dmelanogaster.UCSC.dm6.ensGene); library(org.Dm.eg.db); library(AnnotationDbi)
  })
  txdb   <- TxDb.Dmelanogaster.UCSC.dm6.ensGene
  tss_gr <- suppressMessages(promoters(genes(txdb), upstream = 0, downstream = 1))
  tss_gr <- keepSeqlevels(tss_gr, good_chr, pruning.mode = "coarse")
  tss_gr <- trim(tss_gr)
  tss_gr$flybase <- unname(suppressMessages(
    mapIds(org.Dm.eg.db, keys = names(tss_gr), keytype = "ENSEMBL",
           column = "FLYBASE", multiVals = "first")))
  tss_gr <- tss_gr[!is.na(tss_gr$flybase)]
}
if (!exists("tracheal_sg_genes") || !exists("cns_midline_genes")) {
  REFERENCE_LISTS <- file.path(PROJECT_DIR, "data/reference_gene_lists")
  .read_fbgn <- function(path) {
    df <- if (grepl("\\.csv$", path, ignore.case = TRUE))
      suppressMessages(readr::read_csv(path, show_col_types = FALSE))
    else suppressMessages(readr::read_tsv(path, show_col_types = FALSE))
    v <- if ("genes" %in% names(df)) df[["genes"]] else df[[1]]
    v <- unique(trimws(as.character(v))); v[grepl("^FBgn\\d+$", v)]
  }
  tracheal_sg_genes <- .read_fbgn(file.path(REFERENCE_LISTS, "trh_ts_sg_tissues_reference.csv"))
  cns_midline_genes <- .read_fbgn(file.path(REFERENCE_LISTS, "sim_cns_tissues_DEGs.csv"))
}

# ---- 1. DISJOINT tissue lists (the 347 shared genes are removed) ------------
stopifnot(exists("tracheal_sg_genes"), exists("cns_midline_genes"))
shared  <- intersect(tracheal_sg_genes, cns_midline_genes)
TSGonly <- setdiff(tracheal_sg_genes, cns_midline_genes)
CNSonly <- setdiff(cns_midline_genes, tracheal_sg_genes)
cat(sprintf("TSG %d | CNS %d | shared %d (%.1f%% of CNS) -> TSGonly %d, CNSonly %d\n",
            length(tracheal_sg_genes), length(cns_midline_genes), length(shared),
            100 * length(shared) / length(cns_midline_genes),
            length(TSGonly), length(CNSonly)))

.nearby <- function(gr, maxgap = 1000L) {
  h <- findOverlaps(gr, tss_gr, maxgap = maxgap, ignore.strand = TRUE)
  unique(na.omit(tss_gr$flybase[subjectHits(h)]))
}
NB <- lapply(SETS, .nearby)

# ---- 2. within-row TSG-vs-CNS preference (correctly oriented specificity) ---
spec <- do.call(rbind, lapply(names(SETS), function(s) {
  pg <- NB[[s]]
  a <- sum(TSGonly %in% pg); na <- length(TSGonly)
  b <- sum(CNSonly %in% pg); nb <- length(CNSonly)
  ft <- fisher.test(matrix(c(a, na - a, b, nb - b), nrow = 2, byrow = TRUE))
  data.frame(Set = s,
             Class   = if (s %in% hets) "bHLH-PAS/Tgo" else "non-bHLH-PAS control",
             N_peaks = length(SETS[[s]]),
             pct_TSGonly = round(100 * a / na, 2),
             pct_CNSonly = round(100 * b / nb, 2),
             OR_TSG_over_CNS = round(unname(ft$estimate), 3),
             CI_lo = round(ft$conf.int[1], 3), CI_hi = round(ft$conf.int[2], 3),
             p = signif(ft$p.value, 3), row.names = NULL)
}))
cat("\n=== Within-row tissue preference, disjoint lists ===\n")
print(spec, row.names = FALSE)

# ---- 3. het x tissue interaction: the claim the design actually makes -------
D <- do.call(rbind, lapply(hets, function(h) {
  pg <- NB[[h]]
  rbind(data.frame(set = h, tissue = "TSGonly", prox = TSGonly %in% pg),
        data.frame(set = h, tissue = "CNSonly", prox = CNSonly %in% pg))
}))
D$set    <- factor(D$set, levels = hets)
D$tissue <- factor(D$tissue, levels = c("CNSonly", "TSGonly"))
lrt <- anova(glm(prox ~ set + tissue, D, family = binomial),
             glm(prox ~ set * tissue, D, family = binomial), test = "LRT")
cat(sprintf("\n=== Omnibus het x tissue interaction: LRT = %.3f on %d df, p = %.4f ===\n",
            lrt$Deviance[2], lrt$Df[2], lrt$`Pr(>Chi)`[2]))

pw_int <- do.call(rbind, lapply(combn(hets, 2, simplify = FALSE), function(p) {
  d <- subset(D, set %in% p); d$set <- factor(d$set, levels = p)
  m  <- glm(prox ~ set * tissue, d, family = binomial)
  co <- summary(m)$coefficients; ix <- grep(":", rownames(co))
  ci <- confint.default(m)[ix, ]
  data.frame(Comparison = paste(p, collapse = " vs "),
             Interaction_OR = round(exp(co[ix, 1]), 3),
             CI_lo = round(exp(ci[1]), 3), CI_hi = round(exp(ci[2]), 3),
             p = signif(co[ix, 4], 3), row.names = NULL)
}))
cat("\n=== Pairwise interaction ORs ===\n"); print(pw_int, row.names = FALSE)

# ---- 4. peak-number-matched ranking (removes the peak-count confound) -------
set.seed(20260801)
match_n <- min(vapply(J, length, integer(1)))
matched <- do.call(rbind, lapply(hets, function(h) {
  gr <- J[[h]]
  m <- replicate(500, {
    pg <- .nearby(gr[sample(length(gr), match_n)])
    c(sum(TSGonly %in% pg) / length(TSGonly), sum(CNSonly %in% pg) / length(CNSonly))
  })
  data.frame(Heterodimer = h, N_peaks = length(gr),
             TSGonly_pct = round(100 * mean(m[1, ]), 2), TSG_sd = round(100 * sd(m[1, ]), 2),
             CNSonly_pct = round(100 * mean(m[2, ]), 2), CNS_sd = round(100 * sd(m[2, ]), 2),
             row.names = NULL)
}))
matched$rank_TSG <- rank(-matched$TSGonly_pct)
matched$rank_CNS <- rank(-matched$CNSonly_pct)
cat(sprintf("\n=== Matched to n = %d peaks (500 draws) ===\n", match_n))
print(matched, row.names = FALSE)

# ---- 5. peak-count confound, all seven sets --------------------------------
npk <- vapply(SETS, length, integer(1))
pTS <- vapply(NB, function(pg) 100 * sum(TSGonly %in% pg) / length(TSGonly), numeric(1))
pCN <- vapply(NB, function(pg) 100 * sum(CNSonly %in% pg) / length(CNSonly), numeric(1))
rho_tsg <- cor(npk, pTS, method = "spearman")
rho_cns <- cor(npk, pCN, method = "spearman")
cat(sprintf("\n=== Peak count vs column occupancy: rho_TSG = %.3f, rho_CNS = %.3f (n = 7 sets) ===\n",
            rho_tsg, rho_cns))

# ---- 6. exports ------------------------------------------------------------
EXPORT_DIR <- file.path(PROJECT_DIR, "data/chapter4")
write.csv(spec,    file.path(EXPORT_DIR, "table4_specificity_disjoint.csv"), row.names = FALSE)
write.csv(pw_int,  file.path(EXPORT_DIR, "table4_interaction_pairwise.csv"), row.names = FALSE)
write.csv(matched, file.path(EXPORT_DIR, "table4_peakmatched_ranking.csv"),  row.names = FALSE)

gk <- function(s, col) spec[spec$Set == s, col]
vals <- c(
  specSharedGenes      = length(shared),
  specSharedPctCns     = sprintf("%.1f", 100 * length(shared) / length(cns_midline_genes)),
  specTsgOnlyN         = length(TSGonly),
  specCnsOnlyN         = length(CNSonly),
  specTrhOR            = sprintf("%.2f", gk("Trh/Tgo","OR_TSG_over_CNS")),
  specTrhCIlo          = sprintf("%.2f", gk("Trh/Tgo","CI_lo")),
  specTrhCIhi          = sprintf("%.2f", gk("Trh/Tgo","CI_hi")),
  specSimOR            = sprintf("%.2f", gk("Sim/Tgo","OR_TSG_over_CNS")),
  specSimCIlo          = sprintf("%.2f", gk("Sim/Tgo","CI_lo")),
  specSimCIhi          = sprintf("%.2f", gk("Sim/Tgo","CI_hi")),
  specDysOR            = sprintf("%.2f", gk("Dys/Tgo","OR_TSG_over_CNS")),
  specDysP             = sprintf("%.2f", gk("Dys/Tgo","p")),
  specSimaOR           = sprintf("%.2f", gk("Sima/Tgo","OR_TSG_over_CNS")),
  specSimaP            = sprintf("%.2f", gk("Sima/Tgo","p")),
  specInteractionP     = sprintf("%.2f", lrt$`Pr(>Chi)`[2]),
  specInteractionDev   = sprintf("%.2f", lrt$Deviance[2]),
  specInteractionDf    = lrt$Df[2],
  specMatchedN         = match_n,
  specMatchedTopTsg    = matched$Heterodimer[which.min(matched$rank_TSG)],
  specMatchedTopCns    = matched$Heterodimer[which.min(matched$rank_CNS)],
  specRhoPeaksTsg      = sprintf("%.2f", rho_tsg),
  specRhoPeaksCns      = sprintf("%.2f", rho_cns)
)
save_values("ch4_specificity", vals)

