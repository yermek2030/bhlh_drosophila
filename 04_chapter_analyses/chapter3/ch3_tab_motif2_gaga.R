#!/usr/bin/env Rscript
# ============================================================
# ch3_tab_motif2_gaga.R
# Co-occurrence of Motif2 and GAGA (Trl) sites in Trh/Tgo joint peaks (Chapter 3).
# Reads FIMO output from 02_motif_analysis/cooccurrence/motif2_vs_gaga.R.
# Outputs: results/tables/ch3_motif2_gaga_2x2.tsv, results/reported_values.tsv (section ch3_motif2_gaga).
# ============================================================
PROJECT_DIR <- normalizePath(Sys.getenv("PROJECT_DIR", unset = "."), mustWork = TRUE)
if (!file.exists(file.path(PROJECT_DIR, "README.md")))
  stop("PROJECT_DIR does not point at the repository root: ", PROJECT_DIR)

DATA_DIR      <- file.path(PROJECT_DIR, "data")
RESULTS_DIR   <- file.path(PROJECT_DIR, "results")

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

read_value <- function(name) {
  tab <- read.delim(VALUES_FILE, colClasses = "character", quote = "", comment.char = "")
  tab$value[match(name, tab$name)]
}

AUX   <- DATA_DIR
N_JOINT <- 614L                      # canonical joint peak count

read_fimo <- function(p) {
  d <- read.delim(p, comment.char = "#", stringsAsFactors = FALSE)
  d[!is.na(d$start) & d$q.value < 0.05, , drop = FALSE]
}
pk <- function(d) as.integer(sub("^peak_", "", d$sequence_name))

## ---- 1. Motif2 instances (conservative q < 0.05 scan) --------------------
m2 <- read_fimo(file.path(AUX, "chapter3/fimo_joint/fimo.tsv"))
m2 <- m2[m2$motif_id == "STRCACTGARARAAA", , drop = FALSE]
m2$peak <- pk(m2)

## ---- 2. GAGA instances (matrices enriched in these peaks, same threshold) -
gg <- read_fimo(file.path(AUX, "motif2_gaga/fimo_gaga/fimo.tsv"))
gg$peak <- pk(gg)

stopifnot(max(m2$peak) <= N_JOINT, max(gg$peak) <= N_JOINT)

## ---- 3. peak-level 2x2 ----------------------------------------------------
m2_pk   <- sort(unique(m2$peak))
gaga_pk <- sort(unique(gg$peak))
both    <- length(intersect(m2_pk, gaga_pk))
m2_only <- length(setdiff(m2_pk, gaga_pk))
gg_only <- length(setdiff(gaga_pk, m2_pk))
neither <- N_JOINT - both - m2_only - gg_only
tab <- matrix(c(both, m2_only, gg_only, neither), 2, 2, byrow = TRUE,
              dimnames = list(c("Motif2+","Motif2-"), c("GAGA+","GAGA-")))
ft  <- fisher.test(tab)

stopifnot(sum(tab) == N_JOINT)                       # table closes on 614
stopifnot(abs(unname(ft$estimate) - 2.74) < 0.01)    # reported odds ratio

## ---- 4. instance-level coincidence and separation -------------------------
ov <- 0L; sep <- numeric(0)
for (i in seq_len(nrow(m2))) {
  g <- gg[gg$peak == m2$peak[i], , drop = FALSE]
  if (!nrow(g)) next
  if (any(g$start <= m2$stop[i] & g$stop >= m2$start[i])) ov <- ov + 1L
  sep <- c(sep, min(abs(g$start - m2$start[i])))
}
stopifnot(ov == 4L)                                  # reported overlap count

## ---- 5. peak-level table -------------------------------------------------------
pct <- function(n, d) sprintf("%.1f", 100 * n / d)
gaga_tab <- data.frame(
  row           = c("Motif2-positive peaks", "Motif2-negative peaks", "Total",
                    "Motif2 carriage (%)"),
  gaga_positive = c(tab[1, 1], tab[2, 1], sum(tab[, 1]), pct(tab[1, 1], sum(tab[, 1]))),
  gaga_negative = c(tab[1, 2], tab[2, 2], sum(tab[, 2]), pct(tab[1, 2], sum(tab[, 2]))),
  total         = c(sum(tab[1, ]), sum(tab[2, ]), sum(tab), pct(sum(tab[1, ]), sum(tab))),
  stringsAsFactors = FALSE)
gaga_out <- file.path(RESULTS_DIR, "tables", "ch3_motif2_gaga_2x2.tsv")
dir.create(dirname(gaga_out), recursive = TRUE, showWarnings = FALSE)
write.table(gaga_tab, gaga_out, sep = "\t", quote = FALSE, row.names = FALSE)

## ---- 6. values quoted in the text (section ch3_motif2_gaga) ------------------------------
## Every value is derived from the contingency table built above.

sci2 <- function(p) {                       # 4.603e-06 -> "4.6e-6"
  e <- floor(log10(p)); m <- p / 10^e
  sprintf("%.1fe%d", m, e)
}
pct1 <- function(n, d) sprintf("%.1f", 100 * n / d)

vals <- list(
  nGagaPosPeaks          = sum(tab[, 1]),
  nGagaNegPeaks          = sum(tab[, 2]),
  gagaMotifTwoOR         = sprintf("%.2f", unname(ft$estimate)),
  gagaMotifTwoP          = sci2(ft$p.value),
  nMotifTwoInstancesGaga = nrow(m2),
  nMotifTwoOverlapGaga   = ov,
  pctMotifTwoOverlapGaga = pct1(ov, nrow(m2)),
  nMotifTwoInGagaNeg     = tab[1, 2],
  pctMotifTwoInGagaNeg   = pct1(tab[1, 2], sum(tab[, 2])),
  pctMotifTwoInGagaPos   = pct1(tab[1, 1], sum(tab[, 1])),
  nMotifTwoInGagaPos     = tab[1,1],
  nGagaPosNoMotifTwo     = tab[2,1],
  nGagaNegNoMotifTwo     = tab[2,2],
  nMotifTwoGagaCoPeak    = length(sep),
  motifTwoGagaSepMedian  = sprintf("%.0f", median(sep)),
  motifTwoGagaSepMin     = sprintf("%.0f", min(sep)),
  motifTwoGagaSepMax     = sprintf("%.0f", max(sep))
)

## Stop if a value differs from one already recorded in results/reported_values.tsv.
for (m in names(vals)) {
  old <- read_value(m); new <- as.character(vals[[m]])
  if (is.na(old) || identical(old, new)) next
  on <- suppressWarnings(as.numeric(old)); nn <- suppressWarnings(as.numeric(new))
  if (!is.na(on) && !is.na(nn) && abs(on - nn) < 1e-9) next
  stop(sprintf("%s: derived '%s' but reported_values.tsv has '%s'.", m, new, old))
}
save_values("ch3_motif2_gaga", vals)

cat("\n2x2 peak-level table (joint peaks = ", N_JOINT, ")\n", sep = "")
print(tab)
cat(sprintf("Fisher OR = %.3f (%.3f-%.3f), p = %.4g\n",
            ft$estimate, ft$conf.int[1], ft$conf.int[2], ft$p.value))
cat(sprintf("instances: %d Motif2; %d (%.1f%%) overlap a GAGA hit\n",
            nrow(m2), ov, 100 * ov / nrow(m2)))
cat(sprintf("co-occurring in one peak: n = %d, median separation %.0f bp (range %.0f-%.0f)\n",
            length(sep), median(sep), min(sep), max(sep)))
cat("written:", gaga_out, "and results/reported_values.tsv (section ch3_motif2_gaga)\n")
