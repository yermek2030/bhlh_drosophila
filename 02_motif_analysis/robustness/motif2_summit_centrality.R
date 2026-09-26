# CentriMo centrality analysis at MACS3 signal summits.
# data/final_peaks/TRH_IDR0.05.clean.narrowPeak column 10 is the peak midpoint
# (identical to (start+end)/2 for all 1033 rows), not a summit.
# True summits come from the MACS3 summit files of the two Trh replicates
# (09_peaks/TRH/TRH_rep{1,2}_summits.bed in the processing directory,
# in data/final_peaks/).
## Repository root: the environment variable PROJECT_DIR, or the working directory.
PROJECT_DIR <- normalizePath(Sys.getenv("PROJECT_DIR", unset = "."), mustWork = TRUE)
if (!file.exists(file.path(PROJECT_DIR, "README.md")))
  stop("PROJECT_DIR does not point at the repository root: ", PROJECT_DIR)
suppressPackageStartupMessages({library(GenomicRanges); library(Biostrings)
  library(rtracklayer); library(BSgenome.Dmelanogaster.UCSC.dm6)})
OUT <- file.path(PROJECT_DIR,"data/motif2_robustness")
MEME_BIN <- Sys.getenv("MEME_BIN",""); if (nzchar(MEME_BIN)) Sys.setenv(PATH=paste(MEME_BIN,Sys.getenv("PATH"),sep=":"))

peaks <- import(file.path(PROJECT_DIR,"data/chapter3/Trh_Tgo_joint_peaks_filt.bed"), format="BED")
sm <- do.call(c, lapply(c("TRH_rep1_summits.bed","TRH_rep2_summits.bed"), function(f){
  d <- read.delim(file.path(PROJECT_DIR,"data/final_peaks",f), header=FALSE, stringsAsFactors=FALSE)
  GRanges(d$V1, IRanges(d$V2+1L, d$V3), score=d$V5)}))
h <- findOverlaps(peaks, sm)
# strongest summit per joint peak
best <- tapply(seq_along(subjectHits(h)), queryHits(h),
               function(ix) subjectHits(h)[ix][which.max(sm$score[subjectHits(h)[ix]])])
qi <- as.integer(names(best)); si <- as.integer(best)
cat("joint peaks with a MACS3 TRH summit:", length(qi), "of", length(peaks), "\n")
ctr <- start(sm)[si]
cen <- trim(GRanges(seqnames(peaks)[qi], IRanges(ctr-249L, ctr+250L)))
cen <- cen[width(cen)==500L]
sq <- getSeq(BSgenome.Dmelanogaster.UCSC.dm6, cen); names(sq) <- paste0("summit_",seq_along(sq))
writeXStringSet(sq, file.path(OUT,"trh_TRUEsummits_500bp.fa"))
system(sprintf("centrimo --verbosity 1 --oc %s --local %s %s %s >/dev/null 2>&1",
  shQuote(file.path(OUT,"centrimo_TRUEsummits")), shQuote(file.path(OUT,"trh_TRUEsummits_500bp.fa")),
  shQuote(file.path(PROJECT_DIR,"data/chapter3/STRCACTGARARAAA.meme")),
  shQuote(file.path(OUT,"cme.meme"))))
ct <- read.delim(file.path(OUT,"centrimo_TRUEsummits","centrimo.tsv"), comment.char="#", stringsAsFactors=FALSE)
ct <- ct[nzchar(as.character(ct[[1]])),]
print(ct[,c("motif_id","E.value","adj_p.value","bin_location","bin_width","sites_in_bin","total_sites")], row.names=FALSE)
saveRDS(list(n_summits=length(sq), centrimo=ct), file.path(OUT,"motif2_R11b_truesummits.rds"))
