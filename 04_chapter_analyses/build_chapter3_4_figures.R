## Rebuilds all 15 Chapter 3/4 figures. Run: Rscript --vanilla build_chapter3_4_figures.R
## Figures generated from one chunk per figure in 04_chapter_analyses/chapter3_4_figures.Rmd
## Repository root: PROJECT_DIR environment variable, or working directory
PROJECT_DIR <- normalizePath(Sys.getenv("PROJECT_DIR", unset = "."), mustWork = TRUE)
if (!file.exists(file.path(PROJECT_DIR, "README.md")))
  stop("PROJECT_DIR does not point at the repository root: ", PROJECT_DIR)

RMD <- file.path(PROJECT_DIR, "04_chapter_analyses",
                 "chapter3_4_figures.Rmd")
stopifnot(file.exists(RMD))

## Chunks are eval=FALSE (the file is worked through one figure at a time), so
## read the code blocks and evaluate them directly. Do not use knitr::purl() here:
## it honours eval=FALSE by commenting the code out, and the extracted script
## then runs clean while writing no figures at all.
L   <- readLines(RMD, warn = FALSE)
op  <- grep("^```\\{r", L)
cl  <- grep("^```[[:space:]]*$", L)
lab <- sub("^```\\{r[^']*'([^']*)'\\}.*$", "\\1", L[op])
lab[!grepl("'", L[op])] <- "setup"

env <- new.env(parent = globalenv())
n <- 0L
for (i in seq_along(op)) {
  if (grepl("^Rebuild", lab[i])) next
  cat("---", lab[i], "\n")
  eval(parse(text = L[(op[i] + 1):(cl[cl > op[i]][1] - 1)]), envir = env)
  if (lab[i] != "setup") n <- n + 1L
}
cat("rebuilt", n, "figures into results/figures/\n")
