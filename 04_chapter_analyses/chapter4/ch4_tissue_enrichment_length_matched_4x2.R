## ============================================================
## ch4_tissue_enrichment_length_matched_4x2.R
## Gene-length-matched null (20 quantile bins, 5000 draws) for the 8 cells of the Chapter 4 heterodimer x tissue matrix. Statistic: list genes with TSS within 1kb of a joint peak. Observed counts 114/50, 59/34, 136/48, 57/23.
## Outputs: data/chapter4/lenmatch_4x2/lenmatch_4x2_5000.csv; reported_values.tsv section ch4_length_matched_null.
## Run: Rscript 04_chapter_analyses/chapter4/ch4_tissue_enrichment_length_matched_4x2.R
## ============================================================
PROJECT_DIR <- normalizePath(Sys.getenv("PROJECT_DIR", unset = "."), mustWork = TRUE)
if (!file.exists(file.path(PROJECT_DIR, "README.md")))
  stop("PROJECT_DIR does not point at the repository root: ", PROJECT_DIR)

DATA_DIR      <- file.path(PROJECT_DIR, "data")
REFERENCE_LISTS <- file.path(DATA_DIR, "reference_gene_lists")
CH4_DIR <- file.path(DATA_DIR, "chapter4")

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
  library(GenomicRanges); library(GenomeInfoDb); library(GenomicFeatures)
  library(TxDb.Dmelanogaster.UCSC.dm6.ensGene)
})
OUT <- file.path(CH4_DIR, "lenmatch_4x2")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
std <- c("chr2L", "chr2R", "chr3L", "chr3R", "chr4", "chrX")

# ---- peak sets (Chapter 4 joint sets) -----------------------------------------
rd_peaks <- function(het) {
  d <- read.csv(file.path(CH4_DIR, sprintf("%s_Tgo_joint_annotated.csv", het)))
  GRanges(d$seqnames, IRanges(d$start, d$end))
}
hets <- c(Trh = "Trh/Tgo", Sim = "Sim/Tgo", Dys = "Dys/Tgo", Sima = "Sima/Tgo")
P <- lapply(setNames(names(hets), names(hets)), rd_peaks)
stopifnot(identical(unname(lengths(P)), c(614L, 355L, 870L, 296L)))

# ---- tissue lists (as read in chapter4_heterodimer_comparison.Rmd) ------------
read_fbgn <- function(f) {
  x <- read.csv(file.path(REFERENCE_LISTS, f), stringsAsFactors = FALSE)
  v <- if ("genes" %in% names(x)) x[["genes"]] else x[[1]]
  v <- unique(trimws(as.character(v))); v[grepl("^FBgn[0-9]+$", v)]
}
Lraw <- list(Tsg = read_fbgn("trh_ts_sg_tissues_reference.csv"), Cns = read_fbgn("sim_cns_tissues_DEGs.csv"))
stopifnot(length(Lraw$Tsg) == 2053L, length(Lraw$Cns) == 684L)

# ---- genes, TSS +/- 1 kb windows, hit matrix ---------------------------------
g   <- suppressMessages(genes(TxDb.Dmelanogaster.UCSC.dm6.ensGene)); g <- g[as.character(seqnames(g)) %in% std]
win <- suppressWarnings(trim(promoters(resize(g, 1, fix = "start"), upstream = 1000, downstream = 1001)))
H   <- sapply(P, function(p) countOverlaps(win, p, ignore.strand = TRUE) > 0); rownames(H) <- names(g)
L   <- lapply(Lraw, function(x) intersect(x, names(g)))
obs <- sapply(L, function(x) colSums(H[x, , drop = FALSE]))          # rows = het, cols = list
stopifnot(identical(unname(as.integer(obs[, "Tsg"])), c(114L, 59L, 136L, 57L)),
          identical(unname(as.integer(obs[, "Cns"])), c(50L, 34L, 48L, 23L)))

# ---- length-matched draws (TS/SG first, as in the Chapter 3 script) ----------
len  <- width(g); names(len) <- names(g)
bins <- cut(len, unique(quantile(len, seq(0, 1, 0.05))), include.lowest = TRUE); names(bins) <- names(g)
pools <- split(names(g), bins)
B <- 5000L; set.seed(20260925)
draw_null <- function(ref) {
  need <- table(bins[ref])
  t(vapply(seq_len(B), function(i) {
    s <- unlist(lapply(names(need)[need > 0], function(b) sample(pools[[b]], need[[b]])))
    colSums(H[s, , drop = FALSE]) }, numeric(ncol(H))))
}
N <- list(Tsg = draw_null(L$Tsg), Cns = draw_null(L$Cns))

res <- do.call(rbind, lapply(names(L), function(l) do.call(rbind, lapply(names(hets), function(h) {
  o <- obs[h, l]; nl <- N[[l]][, h]
  data.frame(het = hets[[h]], key = h, list = l, n_genes = length(L[[l]]), peaks = length(P[[h]]),
             observed = o, null_mean = mean(nl), sd = sd(nl), fold = o / mean(nl),
             Z = (o - mean(nl)) / sd(nl), p_emp = (1 + sum(nl >= o)) / (1 + B), n_perm = B)
}))))
res <- res[order(match(res$key, names(hets)), match(res$list, c("Tsg", "Cns"))), ]
print(res[, c("het", "list", "observed", "null_mean", "fold", "Z", "p_emp")], digits = 3, row.names = FALSE)
write.csv(res, file.path(OUT, "lenmatch_4x2_5000.csv"), row.names = FALSE)

# identical to the Chapter 3 null (same seed, loop and list)
trhT <- res[res$key == "Trh" & res$list == "Tsg", ]
stopifnot(abs(trhT$null_mean - 93.20) < 0.005, abs(trhT$sd - 8.35) < 0.005)

# ---- values quoted in the text -------------------------------------------------
f2 <- function(x) formatC(x, format = "f", digits = 2)
fp <- function(p) formatC(p, format = "fg", digits = 2, flag = "#")
cell <- function(h, l) res[res$key == h & res$list == l, ]
save_values("ch4_length_matched_null", list(
  lenMatchFourNperm = B,
  lenMatchFourTsgN  = length(L$Tsg), lenMatchFourCnsN = length(L$Cns),
  lenMatchFourCnsDroppedN = length(Lraw$Cns) - length(L$Cns),
  lenMatchFourFoldMin = f2(min(res$fold)), lenMatchFourFoldMax = f2(max(res$fold)),
  lenMatchFourNsig = sum(res$p_emp < 0.05), lenMatchFourPmax = fp(max(res$p_emp)),
  lenMatchFourSimCnsFold = f2(cell("Sim", "Cns")$fold), lenMatchFourSimCnsP = fp(cell("Sim", "Cns")$p_emp),
  lenMatchFourTrhCnsFold = f2(cell("Trh", "Cns")$fold), lenMatchFourTrhCnsP = fp(cell("Trh", "Cns")$p_emp),
  lenMatchFourSimTsgFold = f2(cell("Sim", "Tsg")$fold), lenMatchFourSimTsgP = fp(cell("Sim", "Tsg")$p_emp),
  lenMatchFourDysTsgFold = f2(cell("Dys", "Tsg")$fold), lenMatchFourDysTsgP = fp(cell("Dys", "Tsg")$p_emp),
  lenMatchFourSimaTsgFold = f2(cell("Sima", "Tsg")$fold), lenMatchFourSimaTsgP = fp(cell("Sima", "Tsg")$p_emp),
  lenMatchFourSimaCnsFold = f2(cell("Sima", "Cns")$fold), lenMatchFourSimaCnsP = fp(cell("Sima", "Cns")$p_emp),
  lenMatchFourDysCnsFold = f2(cell("Dys", "Cns")$fold), lenMatchFourDysCnsP = fp(cell("Dys", "Cns")$p_emp)))
cat("cells with p < 0.05:", sum(res$p_emp < 0.05), "of", nrow(res), "\n")
