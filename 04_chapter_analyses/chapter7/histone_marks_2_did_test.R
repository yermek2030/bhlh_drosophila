# ============================================================
# histone_marks_2_did_test.R
# Difference-in-differences test of histone marks at Motif2-positive vs
# Motif2-negative distal Trh/Tgo joint peaks: DiD = [M2+ (late-early)] - [M2- (late-early)]
# Two-sided Wilcoxon rank-sum on per-peak deltas plus label-permutation null on DiD, N_PERM = 5000
# Reads: data/chapter3/histone_marks_probe_peak.rds
# Writes: data/chapter3/histone_marks_did.rds, data/chapter3/histone_marks_did.csv
# ============================================================

## Repository root: the environment variable PROJECT_DIR, or the working directory.
PROJECT_DIR <- normalizePath(Sys.getenv("PROJECT_DIR", unset = "."), mustWork = TRUE)
if (!file.exists(file.path(PROJECT_DIR, "README.md")))
  stop("PROJECT_DIR does not point at the repository root: ", PROJECT_DIR)
setwd(PROJECT_DIR)

suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
})

set.seed(1234)
N_PERM   <- 5000L
MIN_PROBE <- 2L        # peaks with fewer probes than this are dropped
EARLY    <- "E-8-12"   # pre / early heterodimer activity
LATE     <- "E-12-16h"  # Trh/Tgo ChIP window
OUT_DIR  <- file.path(PROJECT_DIR, "data", "chapter3")

A <- readRDS(file.path(OUT_DIR, "histone_marks_probe_peak.rds"))
peaks <- A$peaks; hits <- A$hits; pct <- A$pct_mat; raw <- A$score_mat
MARKS <- A$mark_levels
stopifnot(length(peaks) == 614L, EARLY %in% A$stage_levels, LATE %in% A$stage_levels)

# ------------------------------------- per-peak mean probe value, per column
qh <- queryHits(hits); sh <- subjectHits(hits)
n_probe <- tabulate(qh, nbins = length(peaks))

peak_means <- function(mat) {
  out <- matrix(NA_real_, nrow = length(peaks), ncol = ncol(mat),
                dimnames = list(peaks$name, colnames(mat)))
  for (j in seq_len(ncol(mat))) {
    s <- vapply(split(mat[sh, j], qh), mean, 0)
    out[as.integer(names(s)), j] <- unname(s)
  }
  out
}
PM_pct <- peak_means(pct)
PM_raw <- peak_means(raw)

# ------------------------------------------------------------- peak universe
keep <- n_probe >= MIN_PROBE & mcols(peaks)$distal
message("distal peaks: ", sum(mcols(peaks)$distal),
        " | with >=", MIN_PROBE, " probes: ", sum(keep))
idx <- which(keep)
m2  <- mcols(peaks)$motif2[idx]
message("  Motif2+: ", sum(m2), "   Motif2-: ", sum(!m2))
stopifnot(sum(m2) >= 20L, sum(!m2) >= 20L)

# ------------------------------------------------------------------ DiD core
did_one <- function(mark, mat) {
  ce <- paste(EARLY, mark, sep = "|"); cl <- paste(LATE, mark, sep = "|")
  stopifnot(ce %in% colnames(mat), cl %in% colnames(mat))
  d  <- mat[idx, cl] - mat[idx, ce]                 # per-peak temporal delta
  obs <- mean(d[m2]) - mean(d[!m2])                 # DiD statistic
  w   <- suppressWarnings(wilcox.test(d[m2], d[!m2], exact = FALSE))
  # label-permutation null on the DiD statistic
  nm2 <- sum(m2)
  perm <- replicate(N_PERM, {
    s <- sample.int(length(d), nm2)
    mean(d[s]) - mean(d[-s])
  })
  p_emp <- (1 + sum(abs(perm) >= abs(obs))) / (N_PERM + 1)
  data.table(mark = mark,
             n_m2_pos = nm2, n_m2_neg = sum(!m2),
             mean_delta_m2_pos = mean(d[m2]),
             mean_delta_m2_neg = mean(d[!m2]),
             did = obs,
             p_wilcox = unname(w$p.value),
             p_empirical = p_emp)
}

res_pct <- rbindlist(lapply(MARKS, did_one, mat = PM_pct))
res_raw <- rbindlist(lapply(MARKS, did_one, mat = PM_raw))
res_pct[, `:=`(scale = "percentile", q_wilcox = p.adjust(p_wilcox, "BH"))]
res_raw[, `:=`(scale = "raw_M",      q_wilcox = p.adjust(p_wilcox, "BH"))]
res <- rbind(res_pct, res_raw)
stopifnot(nrow(res) == 2L * length(MARKS))

# ------------------------- descriptive: absolute percentile by stage and set
stage_tab <- rbindlist(lapply(MARKS, function(mk) rbindlist(lapply(A$stage_levels, function(st) {
  cn <- paste(st, mk, sep = "|")
  if (!cn %in% colnames(PM_pct)) return(NULL)
  v <- PM_pct[idx, cn]
  data.table(mark = mk, stage = st,
             mean_pct_m2_pos = mean(v[m2]), mean_pct_m2_neg = mean(v[!m2]))
}))))

print(res[scale == "percentile", .(mark, did, p_wilcox, q_wilcox, p_empirical)])

saveRDS(list(res = res, stage_tab = stage_tab,
             PM_pct = PM_pct, PM_raw = PM_raw, n_probe = n_probe,
             idx = idx, m2 = m2, early = EARLY, late = LATE,
             n_perm = N_PERM, min_probe = MIN_PROBE),
        file.path(OUT_DIR, "histone_marks_did.rds"))
fwrite(res, file.path(PROJECT_DIR, "data/chapter3/histone_marks_did.csv"))
message("wrote histone_marks_did.rds + .csv")
