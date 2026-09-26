## ch3_occupancy_fold_enrichment.R
## Trh/Tgo occupancy on the input-normalised scale (Chapter 3).
## (1) Promoter-proximal against distal annotated-enhancer peaks: MACS3 fold-enrichment
##     medians and a Wilcoxon rank-sum test (the raw-coverage contrast is Figure 3.15).
## (2) Regression of fold-enrichment on Motif2 dosage with GC content, peak width,
##     log10 TSS distance and CME presence as covariates; HC3 standard errors; binary form;
##     nested F-test. Writes the fitted models that ch3_tab_motif2_regression.R tabulates.
## Inputs:  data/chapter3/trh_joint_macs3_fe.csv (614 joint peaks; signalValue = mean pooled-ChIP
##          coverage, macs3_fold_enrichment = MACS3 fold-enrichment over input at the summit of the
##          IDR-paired replicate peaks, averaged over the two replicates),
##          data/chapter3/Trh_Tgo_joint_annotated.csv, data/chapter3/Trh_Tgo_joint_peaks_filt_raw.fa,
##          data/chapter3/fimo_joint/fimo.tsv, data/reference_gene_lists/merged_embryo_5-20_enhancers_dm6.bed.
## Outputs: data/chapter3/test1_regression_binding.rds; results/reported_values.tsv sections
##          ch3_signal_fold_enrichment and ch3_motif2_regression.
PROJECT_DIR <- normalizePath(Sys.getenv("PROJECT_DIR", unset = "."), mustWork = TRUE)
if (!file.exists(file.path(PROJECT_DIR, "README.md")))
  stop("PROJECT_DIR does not point at the repository root: ", PROJECT_DIR)

DATA_DIR <- file.path(PROJECT_DIR, "data")
CH3_DIR <- file.path(DATA_DIR, "chapter3")
REF_DIR <- file.path(DATA_DIR, "reference_gene_lists")
std_chr <- c("chr2L", "chr2R", "chr3L", "chr3R", "chr4", "chrX")

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
fmt_p <- function(p) {                        # 3 dp at or above 0.001, otherwise 2 s.f. as m x 10^e
  if (p >= 1e-3) return(sprintf("%.3f", p))
  e <- floor(log10(p)); m <- signif(p / 10^e, 2); if (m >= 10) { m <- m / 10; e <- e + 1 }
  sprintf("%.2f \\times 10^{%d}", m, e)
}

suppressPackageStartupMessages({
  library(GenomicRanges); library(GenomeInfoDb); library(Biostrings)
  library(TxDb.Dmelanogaster.UCSC.dm6.ensGene)
  library(sandwich); library(lmtest); library(car)
})

## ---- 1. peaks, covariates, dosage --------------------------------------------
fe_tab <- read.csv(file.path(CH3_DIR, "trh_joint_macs3_fe.csv"))
stopifnot(nrow(fe_tab) == 614L, !anyNA(fe_tab$macs3_fold_enrichment))
joint <- GRanges(fe_tab$seqnames, IRanges(fe_tab$start, fe_tab$end))
ann <- read.csv(file.path(CH3_DIR, "Trh_Tgo_joint_annotated.csv"))
stopifnot(nrow(ann) == 614L, all(ann$start == start(joint)), all(ann$end == end(joint)), sum(ann$CME) == 186L)

tss <- keepSeqlevels(promoters(TxDb.Dmelanogaster.UCSC.dm6.ensGene, upstream = 0, downstream = 1),
                     std_chr, pruning.mode = "coarse")
dh   <- distanceToNearest(joint, tss, ignore.strand = TRUE)
dist <- rep(NA_integer_, 614L); dist[queryHits(dh)] <- mcols(dh)$distance
stopifnot(!anyNA(dist), all(dist >= 0), sum(dist <= 1000) == 468L)

fa <- readDNAStringSet(file.path(CH3_DIR, "Trh_Tgo_joint_peaks_filt_raw.fa"))
stopifnot(length(fa) == 614L, all(width(fa) == width(joint)))

fimo <- read.table(file.path(CH3_DIR, "fimo_joint", "fimo.tsv"), header = TRUE, sep = "\t",
                   comment.char = "#", stringsAsFactors = FALSE)
fimo <- fimo[fimo$motif_id == "STRCACTGARARAAA" & fimo$q.value < 0.05, ]
n_m2 <- tabulate(as.integer(sub("^peak_", "", fimo$sequence_name)), nbins = 614L)
stopifnot(sum(n_m2) == 269L, sum(n_m2 > 0) == 178L)

## ---- 2. promoter against distal annotated enhancer, fold-enrichment scale --------
enh <- read.table(file.path(REF_DIR, "merged_embryo_5-20_enhancers_dm6.bed"), header = FALSE,
                  stringsAsFactors = FALSE)[, 1:3]
enh <- keepSeqlevels(GRanges(enh[[1]], IRanges(enh[[2]], enh[[3]])), std_chr, pruning.mode = "coarse")
P <- dist <= 1000
D <- !P & overlapsAny(joint, enh, ignore.strand = TRUE)
stopifnot(sum(P) == 468L, sum(D) == 114L)
fe <- fe_tab$macs3_fold_enrichment
w  <- wilcox.test(fe[P], fe[D])
cat(sprintf("fold-enrichment median: promoter %.2f (n = %d), distal enhancer %.2f (n = %d); Wilcoxon p = %.3f\n",
            median(fe[P]), sum(P), median(fe[D]), sum(D), w$p.value))
save_values("ch3_signal_fold_enrichment", list(
  sigFEMedPromoter  = sprintf("%.2f", median(fe[P])),
  sigFEMedDistalEnh = sprintf("%.2f", median(fe[D])),
  sigFEPromDistP    = fmt_p(w$p.value)))

## ---- 3. regression of fold-enrichment on Motif2 dosage -------------------------
df <- data.frame(
  peak_score   = fe,
  signal_value = fe_tab$signalValue,
  width        = width(joint),
  dist_tss     = dist,
  has_cme      = as.integer(ann$CME),
  n_motif2     = n_m2,
  gc_content   = as.numeric(letterFrequency(fa, "GC", as.prob = TRUE)))
df$has_motif2   <- as.integer(df$n_motif2 > 0)
df$log_dist_tss <- log10(df$dist_tss + 1)

fit_full    <- lm(peak_score ~ gc_content + width + log_dist_tss + has_cme + n_motif2,   data = df)
fit_reduced <- lm(peak_score ~ gc_content + width + log_dist_tss + has_cme,              data = df)
fit_bin     <- lm(peak_score ~ gc_content + width + log_dist_tss + has_cme + has_motif2, data = df)
fit_log     <- lm(log1p(peak_score) ~ gc_content + width + log_dist_tss + has_cme + n_motif2, data = df)
fit_gc      <- lm(peak_score ~ gc_content + n_motif2, data = df)
## the same model on raw coverage, kept for comparison with Figure 3.15
fit_signal  <- lm(signal_value ~ gc_content + width + log_dist_tss + has_cme + n_motif2, data = df)

av   <- anova(fit_reduced, fit_full)
cf   <- summary(fit_full)$coefficients
cfb  <- summary(fit_bin)$coefficients
p_hc <- coeftest(fit_full, vcov. = vcovHC(fit_full, type = "HC3"))["n_motif2", "Pr(>|t|)"]
stopifnot(abs(cf["n_motif2", "Pr(>|t|)"] - av$`Pr(>F)`[2]) < 1e-12)
print(round(cf, 4))
cat(sprintf("HC3 p = %.3f | binary beta = %.3f, p = %.3f | adjusted R2 %.3f -> %.3f | F(1, %d) = %.2f\n",
            p_hc, cfb["has_motif2", 1], cfb["has_motif2", 4], summary(fit_reduced)$adj.r.squared,
            summary(fit_full)$adj.r.squared, av$Res.Df[2], av$F[2]))
cat(sprintf("raw coverage as outcome: beta = %.2f, p = %.4f\n",
            coef(fit_signal)["n_motif2"], summary(fit_signal)$coefficients["n_motif2", 4]))

saveRDS(list(full = fit_full, reduced = fit_reduced, binary = fit_bin, log = fit_log, gc_only = fit_gc,
             signal = fit_signal, vif = car::vif(fit_full), df = df,
             outcome = "MACS3 fold-enrichment over input; element signal holds the same model on raw coverage"),
        file.path(CH3_DIR, "test1_regression_binding.rds"))

save_values("ch3_motif2_regression", list(
  motifTwoBeta        = sprintf("%.2f", cf["n_motif2", "Estimate"]),
  motifTwoRegP        = fmt_p(cf["n_motif2", "Pr(>|t|)"]),
  motifTwoRegPHC      = fmt_p(p_hc),
  motifTwoBinBeta     = sprintf("%.2f", cfb["has_motif2", "Estimate"]),
  motifTwoBinP        = fmt_p(cfb["has_motif2", "Pr(>|t|)"]),
  adjRReduced         = sprintf("%.3f", summary(fit_reduced)$adj.r.squared),
  adjRFull            = sprintf("%.3f", summary(fit_full)$adj.r.squared),
  motifTwoRegFstat    = sprintf("%.2f", av$F[2]),
  motifTwoRegDfRes    = as.character(av$Res.Df[2]),
  motifTwoRegBetaCIlo = sprintf("%.2f", cf["n_motif2", "Estimate"] - 1.96 * cf["n_motif2", "Std. Error"]),
  motifTwoRegBetaCIhi = sprintf("%.2f", cf["n_motif2", "Estimate"] + 1.96 * cf["n_motif2", "Std. Error"]),
  motifTwoRegOutcomeMedian = sprintf("%.2f", median(df$peak_score))))
