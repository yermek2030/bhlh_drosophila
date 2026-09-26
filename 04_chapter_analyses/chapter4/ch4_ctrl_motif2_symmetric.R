## ============================================================
## ch4_ctrl_motif2_symmetric.R
## Motif2 carriage for control factors, processed as the Trh/Tgo reference (Chapter 4): single, withTgo (pintersect with Tgo), noTgo, each as enhancer-overlapping 201-bp windows, FIMO --bgfile --nrdb-- --thresh 1e-4, >=1 hit. Scans enhancer-matched and window-GC-matched backgrounds; logistic regression on window GC per 10 points.
## Outputs: data/ctrl_tf/symmetric/ (FASTA per set, ctrl_symmetric_motif2.csv/.rds), section ch4_control_panel_symmetric of results/reported_values.tsv.
## Run: Rscript 04_chapter_analyses/chapter4/ch4_ctrl_motif2_symmetric.R
## Repository root: PROJECT_DIR env var, or working directory.
## ============================================================
PROJECT_DIR <- normalizePath(Sys.getenv("PROJECT_DIR", unset = "."), mustWork = TRUE)
if (!file.exists(file.path(PROJECT_DIR, "README.md")))
  stop("PROJECT_DIR does not point at the repository root: ", PROJECT_DIR)

DATA_DIR      <- file.path(PROJECT_DIR, "data")
FINAL_DIR       <- file.path(DATA_DIR, "final_peaks")
REFERENCE_LISTS <- file.path(DATA_DIR, "reference_gene_lists")
CH3_DIR <- file.path(DATA_DIR, "chapter3")
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
## MEME Suite binaries. Resolution order: $MEME_BIN, then PATH.
MEME_BIN <- Sys.getenv("MEME_BIN", unset = "")
if (!nzchar(MEME_BIN)) MEME_BIN <- dirname(Sys.which("fimo")[[1]])
if (!nzchar(MEME_BIN) || !file.exists(file.path(MEME_BIN, "fimo")))
  stop("MEME Suite not found. Install MEME >= 5.5.5 and set MEME_BIN to its bin/ directory.")
suppressPackageStartupMessages({
  library(GenomicRanges); library(GenomeInfoDb); library(Biostrings); library(BSgenome.Dmelanogaster.UCSC.dm6)
})
OUT <- file.path(DATA_DIR, "ctrl_tf", "symmetric"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
g <- BSgenome.Dmelanogaster.UCSC.dm6
std <- c("chr2L", "chr2R", "chr3L", "chr3R", "chr4", "chrX")
PWM <- file.path(CH3_DIR, "STRCACTGARARAAA.meme")

read_np <- function(f) {
  df <- read.table(file.path(FINAL_DIR, f), sep = "\t", header = FALSE, stringsAsFactors = FALSE)
  gr <- GRanges(df[[1]], IRanges(df[[2]] + 1L, df[[3]]))
  sort(keepSeqlevels(gr[as.character(seqnames(gr)) %in% std], std, pruning.mode = "coarse"))
}
P <- lapply(c(TRH = "TRH", TGO = "TGO", TWI = "TWI", VVL = "VVL", RIB = "RIB", DA = "DA"),
            function(x) read_np(paste0(x, "_IDR0.05.final.narrowPeak")))
enh <- read.table(file.path(REFERENCE_LISTS, "merged_embryo_5-20_enhancers_dm6.bed"),
                  header = FALSE, stringsAsFactors = FALSE)[, 1:3]
enh <- GRanges(enh[[1]], IRanges(enh[[2]], enh[[3]]))
enh <- keepSeqlevels(enh[as.character(seqnames(enh)) %in% std], std, pruning.mode = "coarse")

with_tgo <- function(x) {
  ov <- findOverlaps(x, P$TGO, ignore.strand = TRUE)
  j <- pintersect(x[queryHits(ov)], P$TGO[subjectHits(ov)]); j[!duplicated(j)]
}
no_tgo <- function(x) x[countOverlaps(x, P$TGO, ignore.strand = TRUE) == 0L]
win <- function(gr) {
  gr <- gr[countOverlaps(gr, enh, ignore.strand = TRUE) > 0L]
  c0 <- (start(gr) + end(gr)) %/% 2L
  GRanges(seqnames(gr), IRanges(c0 - 100L, c0 + 100L))
}
S <- list()
for (f in c("TRH", "TWI", "VVL", "RIB", "DA")) {
  S[[paste0(f, "_single")]]  <- win(P[[f]])
  S[[paste0(f, "_withTgo")]] <- win(with_tgo(P[[f]]))
  S[[paste0(f, "_noTgo")]]   <- win(no_tgo(P[[f]]))
}
# enhancer-matched background of the control table (whole-enhancer GC/width matched), as scanned
bgb <- read.table(file.path(DATA_DIR, "ctrl_tf/fimo/input_fastas/matched_bg_200bp.bed"), stringsAsFactors = FALSE)
S$BG_enhancer_matched <- GRanges(bgb[[1]], IRanges(bgb[[2]] + 1L, bgb[[3]]))
# window-GC-matched background from enhancers bound by none of the four heterodimers
het <- do.call(c, lapply(c("Trh", "Sim", "Dys", "Sima"), function(h) {
  d <- read.csv(file.path(CH4_DIR, paste0(h, "_Tgo_joint_annotated.csv")))
  GRanges(d$seqnames, IRanges(d$start, d$end)) }))
non_tgo_enh <- enh[countOverlaps(enh, reduce(het), ignore.strand = TRUE) == 0L]
cw <- { c0 <- (start(non_tgo_enh) + end(non_tgo_enh)) %/% 2L; GRanges(seqnames(non_tgo_enh), IRanges(c0 - 100L, c0 + 100L)) }
gcw <- function(gr) as.numeric(letterFrequency(getSeq(g, gr), "GC", as.prob = TRUE))
ref_gc <- gcw(S$TRH_withTgo); cand_gc <- gcw(cw)
br <- unique(quantile(ref_gc, seq(0, 1, 0.1)))
rb <- findInterval(ref_gc, br, rightmost.closed = TRUE); cb <- findInterval(cand_gc, br, rightmost.closed = TRUE)
set.seed(20260925)
idx <- unlist(lapply(unique(rb), function(b) { pool <- which(cb == b); k <- 5L * sum(rb == b)
  if (length(pool) < k) stop("GC bin ", b, " too small"); sample(pool, k) }))
S$BG_gcmatched <- cw[idx]

# the window counts of the Chapter 4 control table are reproduced
print(sapply(S, length))
stopifnot(length(S$TRH_withTgo) == 169L, length(S$TWI_single) == 1892L, length(S$VVL_single) == 258L,
          length(S$RIB_single) == 138L, length(S$BG_enhancer_matched) == 614L)

scan_hits <- function(gr, tag) {
  fa <- file.path(OUT, paste0(tag, ".fa")); s <- getSeq(g, gr); names(s) <- paste0(tag, "_", seq_along(s))
  writeXStringSet(s, fa)
  x <- system2(file.path(MEME_BIN, "fimo"), c("--bgfile", "--nrdb--", "--thresh", "1e-4", "--text", PWM, fa),
               stdout = TRUE, stderr = FALSE)
  x <- read.delim(text = paste(x, collapse = "\n"), stringsAsFactors = FALSE)
  names(s) %in% unique(x$sequence_name)
}
D <- do.call(rbind, lapply(names(S), function(nm) {
  gr <- S[[nm]]; data.frame(set = nm, hit = scan_hits(gr, nm), gc = gcw(gr), stringsAsFactors = FALSE) }))
tab <- do.call(rbind, lapply(split(D, D$set), function(d) data.frame(set = d$set[1], n = nrow(d),
            hits = sum(d$hit), pct = 100 * mean(d$hit), gc = 100 * mean(d$gc))))
ref <- tab[tab$set == "TRH_withTgo", ]
tab$or_vs_trhtgo <- NA; tab$p_vs_trhtgo <- NA
for (i in seq_len(nrow(tab))) if (tab$set[i] != "TRH_withTgo") {
  ft <- fisher.test(matrix(c(tab$hits[i], tab$n[i] - tab$hits[i], ref$hits, ref$n - ref$hits), 2))
  tab$or_vs_trhtgo[i] <- unname(ft$estimate); tab$p_vs_trhtgo[i] <- ft$p.value }
# GC-adjusted logistic regression (window GC per 10 percentage points)
D$gc10 <- 10 * D$gc
D$set <- relevel(factor(D$set), ref = "TRH_withTgo")
fit <- glm(hit ~ set + gc10, family = binomial, data = D)
co <- summary(fit)$coefficients; ci <- confint.default(fit)
adj <- data.frame(set = sub("^set", "", rownames(co)), adj_or = exp(co[, 1]), lo = exp(ci[, 1]), hi = exp(ci[, 2]), p = co[, 4])
adj <- adj[grepl("^set", rownames(co)) | rownames(co) == "gc10", ]
tab <- merge(tab, adj[adj$set != "gc10", ], by = "set", all.x = TRUE)
write.csv(tab, file.path(OUT, "ctrl_symmetric_motif2.csv"), row.names = FALSE)
gc_or <- adj[adj$set == "gc10", ]
print(tab[order(tab$set), c("set", "n", "hits", "pct", "gc", "or_vs_trhtgo", "p_vs_trhtgo", "adj_or", "lo", "hi")], digits = 3, row.names = FALSE)
cat(sprintf("GC per +10 points: OR %.2f (%.2f-%.2f)\n", gc_or$adj_or, gc_or$lo, gc_or$hi))
saveRDS(list(tab = tab, fit = fit, gc_or = gc_or), file.path(OUT, "ctrl_symmetric_motif2.rds"))

# ---- values quoted in the text -------------------------------------------------
f1 <- function(x) formatC(x, format = "f", digits = 1); f2 <- function(x) formatC(x, format = "f", digits = 2)
w  <- tab[tab$set %in% c("TWI_withTgo", "VVL_withTgo", "RIB_withTgo", "DA_withTgo"), ]
nt <- tab[tab$set %in% c("TWI_noTgo", "VVL_noTgo", "RIB_noTgo", "DA_noTgo"), ]
g1 <- function(s, col) tab[tab$set == s, col]
save_values("ch4_control_panel_symmetric", list(
  symTgoPctMin = f1(min(w$pct)), symTgoPctMax = f1(max(w$pct)),
  symTgoPMin   = f2(min(w$p_vs_trhtgo)),
  symTgoAdjMin = f2(min(w$adj_or)), symTgoAdjMax = f2(max(w$adj_or)),
  trhTgoWinGC  = f1(g1("TRH_withTgo", "gc")), bgEnhMatchedGC = f1(g1("BG_enhancer_matched", "gc")),
  gcbgPct = f1(g1("BG_gcmatched", "pct")), gcbgN = g1("BG_gcmatched", "n"), gcbgGC = f1(g1("BG_gcmatched", "gc")),
  gcbgAdjOR = f2(g1("BG_gcmatched", "adj_or")),
  gcOrPerTen = f2(gc_or$adj_or), gcOrLo = f2(gc_or$lo), gcOrHi = f2(gc_or$hi),
  trhNoTgoPct = f1(g1("TRH_noTgo", "pct")), vvlNoTgoPct = f1(g1("VVL_noTgo", "pct")),
  ctrlNoTgoPctMin = f1(min(nt$pct)), ctrlNoTgoPctMax = f1(max(nt$pct)),
  vvlSingleAdjOR = f2(g1("VVL_single", "adj_or")),
  daSingleRebuiltN = g1("DA_single", "n")
))
