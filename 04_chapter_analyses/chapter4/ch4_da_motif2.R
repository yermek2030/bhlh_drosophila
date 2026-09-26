## ch4_da_motif2.R, score Daughterless in the identical frame as the nine
## published single-factor rows (ch4_singlefactor_motif2.R). Read-only: prints
## and writes a CSV, writes no values to results/reported_values.tsv.
## Repository root: the environment variable PROJECT_DIR, or the working directory.
PROJECT_DIR <- normalizePath(Sys.getenv("PROJECT_DIR", unset = "."), mustWork = TRUE)
if (!file.exists(file.path(PROJECT_DIR, "README.md")))
  stop("PROJECT_DIR does not point at the repository root: ", PROJECT_DIR)
suppressPackageStartupMessages({
  library(GenomicRanges); library(GenomeInfoDb); library(Biostrings)
  library(BSgenome.Dmelanogaster.UCSC.dm6)
})
good_chr <- c("chr2L","chr2R","chr3L","chr3R","chr4","chrX")
MOTIF_ID <- "STRCACTGARARAAA"; HALF <- 100L
GENOME   <- BSgenome.Dmelanogaster.UCSC.dm6

enh <- read.table(file.path(PROJECT_DIR,"data/reference_gene_lists",
                            "merged_embryo_5-20_enhancers_dm6.bed"),
                  sep="\t", header=FALSE, stringsAsFactors=FALSE)
enh <- GRanges(enh[[1]], IRanges(enh[[2]]+1L, enh[[3]]))
enh <- keepSeqlevels(enh, good_chr, pruning.mode="coarse")
stopifnot(length(enh) == 35798L)

fimo <- read.delim(file.path(PROJECT_DIR,"data", "chapter3",
                             "fimo_enhancers_genomewide","fimo.tsv"),
                   sep="\t", header=TRUE, comment.char="#", stringsAsFactors=FALSE)
fimo <- fimo[!is.na(fimo$motif_id) & fimo$motif_id == MOTIF_ID, , drop=FALSE]
stopifnot(nrow(fimo) == 15044L)
mid  <- (fimo$start + fimo$stop) %/% 2L
hits <- keepSeqlevels(GRanges(fimo$sequence_name, IRanges(mid,mid)),
                      good_chr, pruning.mode="coarse")

load_peaks <- function(tag) {
  f <- file.path(PROJECT_DIR,"data/final_peaks", sprintf("%s_IDR0.05.final.narrowPeak", tag))
  stopifnot(file.exists(f))
  d  <- read.table(f, sep="\t", header=FALSE, stringsAsFactors=FALSE)
  gr <- GRanges(d[[1]], IRanges(d[[2]]+1L, d[[3]]))
  gr <- gr[seqnames(gr) %in% good_chr]
  gr <- keepSeqlevels(gr, good_chr, pruning.mode="coarse")
  gr[!duplicated(gr)]
}
PUB  <- c("TRH","SIM","DYS","SIMA","TGO","TWI","VVL","RIB")
pk   <- lapply(PUB, load_peaks); names(pk) <- PUB
da   <- load_peaks("DA")
cat(sprintf("DA peaks loaded: %d (of 923 in file; %d dropped off good_chr/dupes)\n",
            length(da), 923L - length(da)))

score_set <- function(gr) {
  ctr  <- (start(gr) + end(gr)) %/% 2L
  win  <- GRanges(seqnames(gr), IRanges(ctr-HALF, ctr+HALF))
  win  <- win[win %within% enh]
  if (!length(win)) return(list(n=0L,k=0L,gc=NA_real_))
  k  <- sum(overlapsAny(win, hits, ignore.strand=TRUE))
  gc <- 100*mean(letterFrequency(getSeq(GENOME,win),"GC",as.prob=TRUE))
  list(n=length(win), k=k, gc=gc)
}
ref <- score_set(pk$TRH)
sda <- score_set(da)

fish <- function(s) {
  m <- matrix(c(s$k, s$n-s$k, ref$k, ref$n-ref$k), 2, 2)
  f <- fisher.test(m); c(or=unname(f$estimate), p=f$p.value,
                         lo=f$conf.int[1], hi=f$conf.int[2])
}
fd <- fish(sda)

## unbound, both definitions
bound8 <- reduce(unlist(GRangesList(unname(pk))), ignore.strand=TRUE)
bound9 <- reduce(c(bound8, da), ignore.strand=TRUE)
u8 <- score_set(enh[!overlapsAny(enh, bound8, ignore.strand=TRUE)])
u9 <- score_set(enh[!overlapsAny(enh, bound9, ignore.strand=TRUE)])

cat("\n=== Daughterless, scored in the published frame ===\n")
cat(sprintf("  peaks with a 201 bp window fully inside an enhancer : %d\n", sda$n))
cat(sprintf("  Motif2 positive                                    : %d (%.1f%%)\n",
            sda$k, 100*sda$k/sda$n))
cat(sprintf("  mean GC                                             : %.1f%%\n", sda$gc))
cat(sprintf("  vs Trh alone (%.1f%%): OR %.2f [%.2f, %.2f], p = %.3g\n",
            100*ref$k/ref$n, fd["or"], fd["lo"], fd["hi"], fd["p"]))

cat("\n=== effect on the 'unbound enhancers' baseline ===\n")
cat(sprintf("  published (8-factor union) : %5d scanned, %.1f%% positive\n", u8$n, 100*u8$k/u8$n))
cat(sprintf("  including DA (9-factor)    : %5d scanned, %.1f%% positive\n", u9$n, 100*u9$k/u9$n))
cat(sprintf("  change in the reported value soloUnbdPct: %.1f -> %.1f\n",
            100*u8$k/u8$n, 100*u9$k/u9$n))

out <- data.frame(Set="Da", TF_family="bHLH class I hub (beta)",
                  N_scanned=sda$n, N_Motif2_pos=sda$k,
                  Pct_Motif2=round(100*sda$k/sda$n,1), Mean_GC_pct=round(sda$gc,1),
                  OR_vs_Trh_alone=round(unname(fd["or"]),2),
                  CI_lo=round(unname(fd["lo"]),2), CI_hi=round(unname(fd["hi"]),2),
                  p_value=signif(unname(fd["p"]),3))
write.csv(out, file.path(PROJECT_DIR,"data", "chapter4","table4_da_motif2.csv"),
          row.names=FALSE)
cat("\ncsv: data/chapter4/table4_da_motif2.csv\n")
