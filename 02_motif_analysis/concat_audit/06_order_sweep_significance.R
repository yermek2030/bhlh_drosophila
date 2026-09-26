## ====
## 06_order_sweep_significance.R: Motif2 significance across STREME background-order sweep
## For each heterodimer and order (0,2,3,4), locates CACTGA motif in streme.txt, reads E-value, counts runs with E<0.05 out of 16
## Input : data/motif2_concat_audit/streme_order_sweep/
## Output: section ch4_order_sweep of results/reported_values.tsv
## ====
PROJECT_DIR <- normalizePath(Sys.getenv("PROJECT_DIR", unset = "."), mustWork = TRUE)
if (!file.exists(file.path(PROJECT_DIR, "README.md")))
  stop("PROJECT_DIR does not point at the repository root: ", PROJECT_DIR)

DATA_DIR      <- file.path(PROJECT_DIR, "data")

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
D <- file.path(DATA_DIR, "motif2_concat_audit", "streme_order_sweep")
r <- do.call(rbind, lapply(c("Trh", "Sim", "Dys", "Sima"), function(h) do.call(rbind, lapply(c(0, 2, 3, 4), function(o) {
  st <- readLines(file.path(D, sprintf("%s_order%d", h, o), "streme.txt")); hd <- grep("^MOTIF ", st)
  cons <- sub("^MOTIF [0-9]+-(\\S+).*$", "\\1", st[hd]); E <- as.numeric(sub(".*E= *([0-9.eE+-]+).*$", "\\1", st[hd + 1]))
  k <- which(grepl("CACTGA|TCAGTG", cons))[1]
  data.frame(set = h, order = o, rank = k, consensus = cons[k], E = E[k]) }))))
print(r, row.names = FALSE)
stopifnot(nrow(r) == 16L, !anyNA(r$rank))
save_values("ch4_order_sweep", list(sweepRunsN = nrow(r), sweepSigN = sum(r$E < 0.05),
                                  sweepTrhSigN = sum(r$E[r$set == "Trh"] < 0.05)))
