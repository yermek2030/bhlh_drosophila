D <- "data/chapter5"
orient <- c(S110387="swap",S110381="swap",S110396="swap",
            S110400="swap",S110405="direct",S110406="direct",
            S110392="swap",S110408="direct",S110404="direct",
            S110388="swap",S110391="swap",S110390="direct",
            S110393="swap",S110383="direct",S110384="direct")
stage_of <- c(S110387="R11",S110381="R11",S110396="R11",S110400="R12",S110405="R12",S110406="R12",
              S110392="R13",S110408="R13",S110404="R13",S110388="R14",S110391="R14",S110390="R14",
              S110393="R15",S110383="R15",S110384="R15")
fr <- read.delim("data/reference_gene_lists/Fisher2023_per_gene_fold_changes.csv",header=TRUE,sep="\t",dec=",",
                 na.strings=c("","NA","NaN"),stringsAsFactors=FALSE,strip.white=TRUE)
names(fr)<-trimws(names(fr)); names(fr)[names(fr)=="gene_id"]<-"sym"; fr$sym<-trimws(fr$sym)
fr <- fr[!duplicated(fr$sym), c("sym","log2FC_6.8")]
cat("PER-SLIDE orientation test: each slide's M vs Fisher log2FC (6-8h).\n")
cat("Correctly oriented slide => POSITIVE rho.\n\n")
res <- list()
for (st in c("R11","R12","R13","R14","R15")) {
  d <- read.delim(file.path(D,paste0(st,".vsn.limma.fdr.txt.gz")),skip=1,stringsAsFactors=FALSE,check.names=FALSE)
  d$sym <- sub(",.*","",d$GeneName)
  d <- d[!duplicated(d$UniqueID),]
  mc <- grep("\\.M$",names(d),value=TRUE)
  for (cc in mc) {
    sl <- sub("\\.M$","",cc)
    x <- data.frame(sym=d$sym, m=suppressWarnings(as.numeric(d[[cc]])), stringsAsFactors=FALSE)
    x <- merge(x, fr, by="sym"); x <- x[complete.cases(x),]
    ct <- suppressWarnings(cor.test(x$m, x$log2FC_6.8, method="spearman"))
    res[[sl]] <- c(rho=unname(ct$estimate), p=ct$p.value, n=nrow(x))
    cat(sprintf("  %-4s %-9s %-7s rho=%+0.4f  p=%.2g  %s\n", st, sl, orient[sl],
        ct$estimate, ct$p.value, if (ct$estimate<0) "<-- NEGATIVE" else ""))
  }
}
r <- sapply(res,function(z) z["rho"])
cat("\nSummary by table orientation:\n")
cat(sprintf("  swap slides   (n=%d): mean rho %+0.4f  [%s]\n", sum(orient=="swap"), mean(r[names(orient)[orient=="swap"]]),
    paste(sprintf("%+0.2f", r[names(orient)[orient=="swap"]]), collapse=" ")))
cat(sprintf("  direct slides (n=%d): mean rho %+0.4f  [%s]\n", sum(orient=="direct"), mean(r[names(orient)[orient=="direct"]]),
    paste(sprintf("%+0.2f", r[names(orient)[orient=="direct"]]), collapse=" ")))
cat("\nIf signs were NOT reconciled, swap and direct groups would have OPPOSITE mean sign.\n")

# ============================================================
# Checks dye-swap orientation is sign-reconciled for trh-RNAi
# microarray before vsn/limma fit, using per-slide M=log2(Cy5/Cy3)
# vs Fisher 2023 RNA-seq reference across 15 slides.
# Result: 13/15 positive, swap vs direct Wilcoxon p=0.40; negatives
# S110396, S110387 (stage 11). Run: Rscript dye_swap_orientation_check.R
# ============================================================
