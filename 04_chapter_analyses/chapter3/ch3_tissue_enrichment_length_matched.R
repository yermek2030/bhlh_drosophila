## ============================================================
## ch3_tissue_enrichment_length_matched.R
## Gene-length-matched null for TS/SG reference gene set (Chapter 3).
## Reference: 114 genes with Trh/Tgo joint peak within 1kb of TSS. Null:
## peaks fixed, gene sets drawn matching reference on gene length (20
## quantile bins), 5000 times, for full set and each data/reference_sets/ref_*.txt.
## Outputs: data/reference_sets/refset_lenmatched_5000.csv, section ch3_length_matched_null
##          of results/reported_values.tsv
## ============================================================
PROJECT_DIR <- normalizePath(Sys.getenv("PROJECT_DIR", unset = "."), mustWork = TRUE)
if (!file.exists(file.path(PROJECT_DIR, "README.md")))
  stop("PROJECT_DIR does not point at the repository root: ", PROJECT_DIR)

DATA_DIR      <- file.path(PROJECT_DIR, "data")
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
  library(GenomicRanges); library(GenomeInfoDb); library(GenomicFeatures); library(TxDb.Dmelanogaster.UCSC.dm6.ensGene)
})
A <- file.path(DATA_DIR, "reference_sets")
std <- c("chr2L", "chr2R", "chr3L", "chr3R", "chr4", "chrX")
d  <- read.csv(file.path(CH4_DIR, "Trh_Tgo_joint_annotated.csv"))
jp <- GRanges(d$seqnames, IRanges(d$start, d$end)); stopifnot(length(jp) == 614L)
g  <- suppressMessages(genes(TxDb.Dmelanogaster.UCSC.dm6.ensGene)); g <- g[as.character(seqnames(g)) %in% std]
win <- suppressWarnings(trim(promoters(resize(g, 1, fix = "start"), upstream = 1000, downstream = 1001)))
hit <- names(g)[countOverlaps(win, jp, ignore.strand = TRUE) > 0]
rd  <- function(f) { x <- trimws(readLines(file.path(A, f))); intersect(x[nzchar(x)], names(g)) }
parts <- list(full = rd("ref2053.txt"), trach_any = rd("ref_trach.txt"), sg_any = rd("ref_sg.txt"),
              trach_only = rd("ref_trach_only.txt"), sg_only = rd("ref_sg_only.txt"),
              both = rd("ref_both.txt"), screen_only = rd("ref_screen_only.txt"))
stopifnot(sum(parts$full %in% hit) == 114L)
len <- width(g); names(len) <- names(g)
bins <- cut(len, unique(quantile(len, seq(0, 1, 0.05))), include.lowest = TRUE); names(bins) <- names(g)
is_hit <- names(g) %in% hit; names(is_hit) <- names(g)
B <- 5000L; set.seed(20260925)
res <- do.call(rbind, lapply(names(parts), function(p) {
  ref <- parts[[p]]; need <- table(bins[ref]); pools <- split(names(g), bins)
  obs <- sum(is_hit[ref])
  null <- vapply(seq_len(B), function(i) {
    s <- unlist(lapply(names(need)[need > 0], function(b) sample(pools[[b]], need[[b]])))
    sum(is_hit[s]) }, numeric(1))
  data.frame(partition = p, n_genes = length(ref), observed = obs, null_mean = mean(null),
             fold = obs / mean(null), sd = sd(null), Z = (obs - mean(null)) / sd(null),
             p_emp = (1 + sum(null >= obs)) / (1 + B), n_perm = B)
}))
print(res, digits = 3, row.names = FALSE)
write.csv(res, file.path(A, "refset_lenmatched_5000.csv"), row.names = FALSE)
key <- c(full = "Full", trach_any = "TrachAny", sg_any = "SgAny", trach_only = "TrachOnly",
         sg_only = "SgOnly", both = "Both", screen_only = "ScreenOnly")
f2 <- function(x) formatC(x, format = "f", digits = 2)
vals <- list()
for (i in seq_len(nrow(res))) {
  k <- key[[res$partition[i]]]
  vals[[paste0("lenMatch", k, "N")]]    <- res$n_genes[i]
  vals[[paste0("lenMatch", k, "Obs")]]  <- res$observed[i]
  vals[[paste0("lenMatch", k, "Null")]] <- f2(res$null_mean[i])
  vals[[paste0("lenMatch", k, "Fold")]] <- f2(res$fold[i])
  vals[[paste0("lenMatch", k, "Z")]]    <- f2(res$Z[i])
  vals[[paste0("lenMatch", k, "Sd")]]   <- f2(res$sd[i])
  vals[[paste0("lenMatch", k, "P")]]    <- formatC(res$p_emp[i], format = "fg", digits = 2, flag = "#")
}
vals$lenMatchNperm <- B
ref_raw <- trimws(readLines(file.path(A, "ref2053.txt"))); ref_raw <- unique(ref_raw[nzchar(ref_raw)])
vals$lenMatchDroppedN <- length(ref_raw) - res$n_genes[res$partition == "full"]
save_values("ch3_length_matched_null", vals)
