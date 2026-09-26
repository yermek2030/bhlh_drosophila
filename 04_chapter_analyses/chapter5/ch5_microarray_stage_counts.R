## ch5_microarray_stage_counts.R
## Stage-wise counts of probes passing Benjamini-Hochberg thresholds in the in-house
## trh RNAi microarray (Chapter 5), and the minimum detectable effect at stage 13.
## Inputs:  data/chapter5/R11.vsn.limma.fdr.txt.gz ... R15.vsn.limma.fdr.txt.gz
##          (limma output after vsn normalisation, one row per TargetID).
## Outputs: results/reported_values.tsv section ch5_microarray_stage_counts.
PROJECT_DIR <- normalizePath(Sys.getenv("PROJECT_DIR", unset = "."), mustWork = TRUE)
if (!file.exists(file.path(PROJECT_DIR, "README.md")))
  stop("PROJECT_DIR does not point at the repository root: ", PROJECT_DIR)

DATA_DIR <- file.path(PROJECT_DIR, "data")
CH5_DIR <- file.path(DATA_DIR, "chapter5")

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

stages <- 11:15
tab <- do.call(rbind, lapply(stages, function(s) {
  d <- read.delim(gzfile(file.path(CH5_DIR, sprintf("R%d.vsn.limma.fdr.txt.gz", s))),
                  comment.char = "#", check.names = FALSE)
  stopifnot("adj.P.Val" %in% names(d), "TargetID" %in% names(d), !anyDuplicated(d$TargetID))
  data.frame(stage = s, probes = nrow(d),
             n05 = sum(d$adj.P.Val < 0.05, na.rm = TRUE),
             n10 = sum(d$adj.P.Val < 0.10, na.rm = TRUE),
             n25 = sum(d$adj.P.Val < 0.25, na.rm = TRUE))
}))
print(tab)
stopifnot(all(tab$probes == 14441L))
st13  <- tab[tab$stage == 13, ]
other <- tab[tab$stage != 13, ]

## Minimum detectable effect at stage 13: the smallest |t| among probes below the
## threshold, converted to a fold change with the median residual standard error
## implied by logFC / t.
d13  <- read.delim(gzfile(file.path(CH5_DIR, "R13.vsn.limma.fdr.txt.gz")), comment.char = "#", check.names = FALSE)
keep <- is.finite(d13$t) & d13$t != 0 & is.finite(d13$logFC)
se     <- median(abs(d13$logFC[keep] / d13$t[keep]))
t_bh   <- min(abs(d13$t[d13$adj.P.Val < 0.05]), na.rm = TRUE)
t_bonf <- min(abs(d13$t[d13$P.Value < 0.05 / nrow(d13)]), na.rm = TRUE)
cat(sprintf("stage 13: |t| >= %.2f at BH 0.05 (%.2f-fold); |t| >= %.2f at Bonferroni (%.2f-fold)\n",
            t_bh, 2^(t_bh * se), t_bonf, 2^(t_bonf * se)))

save_values("ch5_microarray_stage_counts", list(
  nMarrayBhStThirteen      = st13$n05,
  nMarrayBhStThirteenTen   = st13$n10,
  nMarrayBhStThirteenQuart = st13$n25,
  nMarrayBhOtherMax        = max(other$n05),
  nMarrayProbeSet          = st13$probes,
  marrayCritTStThirteen    = sprintf("%.1f", t_bh),
  marrayMdeStThirteen      = sprintf("%.1f", 2^(t_bh * se)),
  marrayMdeBonfStThirteen  = sprintf("%.1f", 2^(t_bonf * se))))
