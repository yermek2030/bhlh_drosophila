## ============================================================
##  ch3_tab_motif2_regression.R
##  Table: regression of Trh/Tgo occupancy on Motif2 dosage, controlling for
##  sequence composition, peak width, TSS distance and CME presence.
##  Input  : data/chapter3/test1_regression_binding.rds (saved model, not refit)
##  Output : results/tables/ch3_motif2_regression_coefficients.tsv
## ============================================================

## Repository root: the environment variable PROJECT_DIR, or the working directory.
PROJECT_DIR <- normalizePath(Sys.getenv("PROJECT_DIR", unset = "."), mustWork = TRUE)
if (!file.exists(file.path(PROJECT_DIR, "README.md")))
  stop("PROJECT_DIR does not point at the repository root: ", PROJECT_DIR)

RESULTS_DIR   <- file.path(PROJECT_DIR, "results")

AX  <- file.path(PROJECT_DIR, "data", "chapter3")
OUT <- file.path(RESULTS_DIR, "tables", "ch3_motif2_regression_coefficients.tsv")
dir.create(dirname(OUT), recursive = TRUE, showWarnings = FALSE)

t1 <- readRDS(file.path(AX, "test1_regression_binding.rds"))

cf  <- summary(t1$full)$coefficients
vif <- t1$vif
av  <- anova(t1$reduced, t1$full)

## The single-coefficient t-test and the nested F-test are the same test on
## one degree of freedom; assert that identity.
stopifnot(abs(cf["n_motif2", "Pr(>|t|)"] - av$`Pr(>F)`[2]) < 1e-12)

## Row order: predictor of interest first, then covariates in the order they
## enter the model. Intercept omitted, it is not interpretable here and
## occupies a row a reader would have to skip.
rows <- c(n_motif2 = "Motif2 dosage (instances per peak)",
          gc_content = "GC content",
          width = "Peak width (bp)",
          log_dist_tss = "log10 distance to nearest TSS",
          has_cme = "CME present")

## Plain numerals, no thousands separator.
num <- function(x, dp) formatC(x, format = "f", digits = dp, big.mark = "")

fmt_p <- function(p) {
  if (p >= 0.001) return(sprintf("%.3f", p))
  e <- floor(log10(p)); m <- p / 10^e
  sprintf("%.2fe%d", m, as.integer(e))
}

coef_tab <- data.frame(
  predictor = unname(rows),
  beta      = vapply(names(rows), function(k) num(cf[k, "Estimate"],   if (k == "width") 4 else 3), ""),
  se        = vapply(names(rows), function(k) num(cf[k, "Std. Error"], if (k == "width") 4 else 3), ""),
  t         = vapply(names(rows), function(k) num(cf[k, "t value"], 2), ""),
  p         = vapply(names(rows), function(k) fmt_p(cf[k, "Pr(>|t|)"]), ""),
  vif       = vapply(names(rows), function(k) num(vif[[k]], 2), ""),
  stringsAsFactors = FALSE)

write.table(coef_tab, OUT, sep = "\t", quote = FALSE, row.names = FALSE)

cat(sprintf("model n            : %d   residual df: %d\n",
            nrow(t1$df), df.residual(t1$full)))
cat(sprintf("adj R2 reduced/full: %.4f / %.4f\n",
            summary(t1$reduced)$adj.r.squared, summary(t1$full)$adj.r.squared))
cat(sprintf("nested F(1,%d)     : %.3f   p = %.4g\n",
            av$Res.Df[2], av$F[2], av$`Pr(>F)`[2]))
cat(sprintf("max VIF            : %.3f\n", max(vif)))
cat(sprintf("written: %s\n", OUT))
