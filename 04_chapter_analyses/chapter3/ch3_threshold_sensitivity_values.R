#!/usr/bin/env Rscript
## ============================================================
## ch3_threshold_sensitivity_values.R: values for Chapter 3 threshold sensitivity and passive-accessibility control section.
## Input: data/chapter3/threshold_sensitivity_and_accessibility.rds (boundary sweep, threshold-free logistic models, A/T-content absorption test, Mantel-Haenszel test, Spearman correlation).
## Output: results/reported_values.tsv
## Figures drawn by chapter3_4_figures.Rmd from same object. Repository root: PROJECT_DIR env var, or working directory.
## ============================================================
PROJECT_DIR <- normalizePath(Sys.getenv("PROJECT_DIR", unset = "."), mustWork = TRUE)
if (!file.exists(file.path(PROJECT_DIR, "README.md")))
  stop("PROJECT_DIR does not point at the repository root: ", PROJECT_DIR)

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

.t4 <- readRDS(file.path(PROJECT_DIR,
         "data/chapter3/threshold_sensitivity_and_accessibility.rds"))

sw <- .t4$sweep                 # 8-row boundary sweep: cut_bp, or, ci, p, n_distal
a1 <- .t4$at_absorb$a1          # motif2 ~ at
b0 <- .t4$at_absorb$b0          # motif2 ~ log10d
b1 <- .t4$at_absorb$b1          # motif2 ~ log10d + at
b2 <- .t4$at_absorb$b2          # motif2 ~ log10d + at + width + has_cme
lrt_p   <- .t4$at_absorb$lrt_p  # likelihood-ratio test for adding A/T content
mh_test <- .t4$mh               # Mantel-Haenszel test across GC-content quartiles
sp_rho  <- .t4$spearman$estimate

## Boundaries retaining fewer than thirty distal peaks are reported but treated
## as underpowered.
MIN_DISTAL <- 30

.f  <- function(x, digits = 2) sprintf(paste0("%.", digits, "f"), x)
.fp <- function(p) {
  if (p >= 0.01) return(sprintf("%.2f", p))
  e <- floor(log10(p))
  sprintf("%.2fe%d", p / 10^e, e)
}

.ok <- sw[sw$n_distal >= MIN_DISTAL, ]

vals <- c(
  sweepORmin      = .f(min(.ok$OR)),
  sweepORmax      = .f(max(.ok$OR)),
  sweepCutMin     = as.character(min(.ok$cut_bp)),          # bp
  sweepCutMax     = as.character(max(.ok$cut_bp) / 1000),   # kb, powered range
  sweepCutTop     = as.character(max(sw$cut_bp) / 1000),    # kb, all 8 cuts
  sweepMinDistal  = as.character(MIN_DISTAL),
  sweepORoneK     = .f(sw$OR[sw$cut_bp == 1000]),
  sweepPoneK      = .fp(sw$p[sw$cut_bp == 1000]),
  sweepORtwoK     = .f(sw$OR[sw$cut_bp == 2000]),
  sweepPtwoK      = .fp(sw$p[sw$cut_bp == 2000]),
  contORperDecade = .f(exp(coef(b0)["log10d"])),
  contPvalue      = .fp(summary(b0)$coefficients["log10d", 4]),
  contSpearmanRho = .f(as.numeric(sp_rho)),
  atOnlyP         = .f(summary(a1)$coefficients["at", 4]),
  atLRTP          = .f(lrt_p),
  contORadjAT     = .f(exp(coef(b1)["log10d"])),
  mhCommonOR      = .f(as.numeric(mh_test$estimate)),
  mhCommonP       = .fp(mh_test$p.value)
)

print(sw)
print(vals)
save_values("ch3", as.list(vals))
