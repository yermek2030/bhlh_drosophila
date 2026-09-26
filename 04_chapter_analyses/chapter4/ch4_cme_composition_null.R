## Composition nulls for CME carriage in four heterodimer joint peak sets.
## Shuffle null: dinucleotide-shuffled peaks (Euler, k=2, universalmotif), B=2000, counts peaks with >=1 ACGTG either strand.
## Random null: width-matched intervals on same chromosome, N-free, B=2000.
## Input: data/chapter4/*_Tgo_joint_annotated.csv (614/355/870/296 peaks).
## Output: data/chapter4/cme_composition_null/cme_composition_null.csv; results/reported_values.tsv section ch4_cme_composition_null.
## Repository root: PROJECT_DIR env var or working directory.
PROJECT_DIR <- normalizePath(Sys.getenv("PROJECT_DIR", unset = "."), mustWork = TRUE)
if (!file.exists(file.path(PROJECT_DIR, "README.md")))
  stop("PROJECT_DIR does not point at the repository root: ", PROJECT_DIR)

DATA_DIR      <- file.path(PROJECT_DIR, "data")
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
suppressPackageStartupMessages({
  library(GenomicRanges); library(Biostrings); library(universalmotif)
  library(BSgenome.Dmelanogaster.UCSC.dm6)
})
OUT <- file.path(CH4_DIR, "cme_composition_null"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
g <- BSgenome.Dmelanogaster.UCSC.dm6
B <- 2000L
canon <- c(Trh = 186L, Sim = 112L, Dys = 282L, Sima = 103L)   # CME-positive peaks (Chapters 3 and 4)
has_cme <- function(x) (vcountPattern("ACGTG", x, fixed = TRUE) + vcountPattern("CACGT", x, fixed = TRUE)) > 0
res <- list()
for (h in names(canon)) {
  df <- read.csv(file.path(CH4_DIR, paste0(h, "_Tgo_joint_annotated.csv")))
  gr <- GRanges(df$seqnames, IRanges(df$start, df$end))
  s  <- getSeq(g, gr)
  obs <- sum(has_cme(s)); stopifnot(obs == canon[[h]])
  set.seed(20260925)
  shuf <- vapply(seq_len(B), function(i) sum(has_cme(shuffle_sequences(s, k = 2, method = "euler"))), numeric(1))
  # random width-matched, chromosome-matched, N-free intervals
  chr <- as.character(seqnames(gr)); w <- width(gr); L <- seqlengths(g)[chr]
  draw <- function() {
    st <- floor(runif(length(w), 1, L - w + 2)); x <- getSeq(g, GRanges(chr, IRanges(st, width = w)))
    bad <- which(letterFrequency(x, "N") > 0)
    while (length(bad)) {
      st[bad] <- floor(runif(length(bad), 1, L[bad] - w[bad] + 2))
      x[bad] <- getSeq(g, GRanges(chr[bad], IRanges(st[bad], width = w[bad])))
      bad <- bad[letterFrequency(x[bad], "N") > 0]
    }
    sum(has_cme(x))
  }
  set.seed(20260926)
  rnd <- vapply(seq_len(B), function(i) draw(), numeric(1))
  p_lo <- (1 + sum(shuf <= obs)) / (1 + B)
  p_rnd <- min(1, 2 * min((1 + sum(rnd <= obs)) / (1 + B), (1 + sum(rnd >= obs)) / (1 + B)))
  res[[h]] <- data.frame(set = h, n = length(gr), obs = obs, shuf_mean = mean(shuf), shuf_sd = sd(shuf),
                         z = (obs - mean(shuf)) / sd(shuf), p_shuf_lower = p_lo,
                         rnd_mean = mean(rnd), p_rnd_two_sided = p_rnd)
  cat(sprintf("%-5s obs %3d | shuffle %.1f (sd %.1f) Z %.2f p %.4f | random %.1f p %.3f\n",
              h, obs, mean(shuf), sd(shuf), res[[h]]$z, p_lo, mean(rnd), p_rnd))
}
R <- do.call(rbind, res); write.csv(R, file.path(OUT, "cme_composition_null.csv"), row.names = FALSE)
f1 <- function(x) formatC(x, format = "f", digits = 1); f2 <- function(x) formatC(x, format = "f", digits = 2)
fp <- function(p) if (p < 0.001) "0.001" else formatC(p, format = "f", digits = 3)
v <- list(cmeNullB = B, cmeShufPmax = fp(max(R$p_shuf_lower)),
          cmeRandPmin = f2(min(R$p_rnd_two_sided)), cmeRandPmax = f2(max(R$p_rnd_two_sided)),
          cmeShufPctExpTrh = f1(100 * R$shuf_mean[R$set == "Trh"] / R$n[R$set == "Trh"]),
          cmeShufFoldTrh = f2(R$obs[R$set == "Trh"] / R$shuf_mean[R$set == "Trh"]))
for (h in names(canon)) {
  v[[paste0("cmeShufExp", h)]] <- f1(R$shuf_mean[R$set == h])
  v[[paste0("cmeShufZ", h)]]   <- f2(R$z[R$set == h])
}
save_values("ch4_cme_composition_null", v)
