## ============================================================
## ch5_rule_on_control_factors.R: benchmarks Chapter 5 direct-target rule (nearest TSS within +/-2kb of a peak, gene in Fisher et al. 2023 combined DEG list) on the four control factors alone, intersected with Tgo, and the other three Tgo heterodimers; universe is 11120 genes.
## Inputs: data/final_peaks/, data/chapter5/Fisher_combined_DEGs.csv, data/chapter5/merged_microarray_fisher_st13.csv. Outputs: data/chapter5/ctrl_rule_ben
PROJECT_DIR <- normalizePath(Sys.getenv("PROJECT_DIR", unset = "."), mustWork = TRUE)
if (!file.exists(file.path(PROJECT_DIR, "README.md")))
  stop("PROJECT_DIR does not point at the repository root: ", PROJECT_DIR)

DATA_DIR      <- file.path(PROJECT_DIR, "data")
FINAL_DIR       <- file.path(DATA_DIR, "final_peaks")
CH5_DIR <- file.path(DATA_DIR, "chapter5")

good_chr <- c("chr2L", "chr2R", "chr3L", "chr3R", "chr4", "chrX")   # dm6 major arms

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
  library(GenomicRanges); library(GenomeInfoDb); library(rtracklayer); library(ChIPseeker)
  library(TxDb.Dmelanogaster.UCSC.dm6.ensGene); library(org.Dm.eg.db)
})
OUT <- file.path(CH5_DIR, "ctrl_rule_benchmark")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
txdb <- TxDb.Dmelanogaster.UCSC.dm6.ensGene

# ---- peak sets: Chapter 4/5 recipe (chr filter, dedup, pintersect with Tgo) --
clean <- function(f) {
  # read.table, not import(): some control files carry a non-integer 10th column
  df <- read.table(file.path(FINAL_DIR, paste0(f, "_IDR0.05.final.narrowPeak")), sep = "\t",
                   header = FALSE, stringsAsFactors = FALSE)
  gr <- GRanges(df[[1]], IRanges(df[[2]] + 1L, df[[3]]))   # BED 0-based -> 1-based, as import() does
  gr <- keepSeqlevels(gr[seqnames(gr) %in% good_chr], good_chr, pruning.mode = "coarse")
  gr[!duplicated(gr)]
}
P <- lapply(c(TRH = "TRH", TGO = "TGO", TWI = "TWI", VVL = "VVL", RIB = "RIB", DA = "DA",
              SIM = "SIM", DYS = "DYS", SIMA = "SIMA"), clean)
with_tgo <- function(x) {
  ov <- findOverlaps(x, P$TGO, minoverlap = 1, ignore.strand = TRUE)
  j <- pintersect(x[queryHits(ov)], P$TGO[subjectHits(ov)])
  j <- keepSeqlevels(j, good_chr, pruning.mode = "coarse"); j[!duplicated(j)]
}
S <- list(`Trh/Tgo` = with_tgo(P$TRH), `Sim/Tgo` = with_tgo(P$SIM), `Dys/Tgo` = with_tgo(P$DYS),
          `Sima/Tgo` = with_tgo(P$SIMA), Twi = P$TWI, Vvl = P$VVL, Rib = P$RIB, Da = P$DA,
          `Twi/Tgo` = with_tgo(P$TWI), `Vvl/Tgo` = with_tgo(P$VVL), `Rib/Tgo` = with_tgo(P$RIB),
          `Da/Tgo` = with_tgo(P$DA))
stopifnot(length(S$`Trh/Tgo`) == 614L)

# ---- the Chapter 5 gene assignment --------------------------------------------
genes_of <- function(gr) {
  a <- as.data.frame(annotatePeak(gr, TxDb = txdb, tssRegion = c(-2000L, 2000L), level = "gene",
                                  annoDb = "org.Dm.eg.db", verbose = FALSE))
  a <- subset(a, abs(distanceToTSS) <= 2000L)
  fb <- intersect(c("FLYBASE", "flybase", "ENSEMBL", "geneId"), names(a))[1]
  g <- unique(a[[fb]]); g[!is.na(g) & nzchar(g)]
}
G <- lapply(S, genes_of)

# ---- DEG list and testable universe (as in chapter5_direct_targets.Rmd) -------
degs <- read.csv(file.path(CH5_DIR, "Fisher_combined_DEGs.csv"), stringsAsFactors = FALSE)
degs <- unique(as.character(degs[[grep("FBgn|gene|flybase", names(degs), ignore.case = TRUE)[1]]]))
m13  <- read.csv(file.path(CH5_DIR, "merged_microarray_fisher_st13.csv"), stringsAsFactors = FALSE)
U    <- unique(m13$FBgn[!is.na(m13$FBgn)])
stopifnot(length(degs) == 934L, length(U) == 11120L, length(G$`Trh/Tgo`) == 448L,
          length(intersect(G$`Trh/Tgo`, degs)) == 21L)

KU <- length(intersect(degs, U)); N <- length(U)
res <- do.call(rbind, lapply(names(G), function(s) {
  g <- G[[s]]; gu <- intersect(g, U); k <- length(intersect(gu, degs)); n <- length(gu)
  ft <- fisher.test(matrix(c(k, n - k, KU - k, N - KU - n + k), 2, 2))
  data.frame(set = s, peaks = length(S[[s]]), genes = length(g), targets = length(intersect(g, degs)),
             genes_u = n, targets_u = k, expected = n * KU / N, OR = unname(ft$estimate),
             lo = ft$conf.int[1], hi = ft$conf.int[2],
             p_enrich = phyper(k - 1, KU, N - KU, n, lower.tail = FALSE))
}))
print(res, digits = 3, row.names = FALSE)
write.csv(res, file.path(OUT, "ch5_rule_on_control_factors.csv"), row.names = FALSE)

# ---- are the 21 Trh/Tgo candidates also 'targets' of the control factors? ----
dt <- intersect(G$`Trh/Tgo`, degs)
nb <- rowSums(sapply(G[c("Twi", "Vvl", "Rib", "Da")], function(g) dt %in% g))
shared <- data.frame(FBgn = dt, n_control_factors = nb)
write.csv(shared, file.path(OUT, "trh_candidates_bound_by_controls.csv"), row.names = FALSE)

# ---- values quoted in the text -------------------------------------------------
ctrl <- res[!res$set %in% c("Trh/Tgo", "Sim/Tgo", "Dys/Tgo", "Sima/Tgo"), ]
f2 <- function(x) formatC(x, format = "f", digits = 2)
save_values("ch5_rule_on_controls", list(
  ctrlRuleTargetsMin = min(ctrl$targets), ctrlRuleTargetsMax = max(ctrl$targets),
  ctrlRuleORmin = f2(min(ctrl$OR)), ctrlRuleORmax = f2(max(ctrl$OR)),
  ctrlRulePmin = formatC(min(ctrl$p_enrich), format = "fg", digits = 2, flag = "#"),
  ctrlRuleTwiTargets = res$targets[res$set == "Twi"],
  ctrlRuleTwiGenes = res$genes[res$set == "Twi"],
  hetRuleORmin = f2(min(res$OR[res$set %in% c("Sim/Tgo", "Dys/Tgo", "Sima/Tgo")])),
  hetRuleORmax = f2(max(res$OR[res$set %in% c("Sim/Tgo", "Dys/Tgo", "Sima/Tgo")])),
  dtBoundByCtrlAny = sum(nb >= 1), dtBoundByCtrlTwo = sum(nb >= 2)))
cat("candidates bound by >=1 / >=2 control factors:", sum(nb >= 1), sum(nb >= 2), "\n")
