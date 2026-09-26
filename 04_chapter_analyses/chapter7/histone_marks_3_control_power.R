# ============================================================
# histone_marks_3_control_power.R
# Positive control: H3K4me3 must be higher at promoter-proximal peaks
# (|dTSS| <= 1 kb) than distal peaks. Power: simulate additive shift on
# Motif2 deltas, record rejection rate at alpha = 0.05 to get minimum
# detectable effect for the achieved sample size.
# Reads: histone_marks_probe_peak.rds, histone_marks_did.rds
# Writes: data/chapter3/histone_marks_control_power.rds, .csv
# ============================================================

## Repository root: the environment variable PROJECT_DIR, or the working directory.
PROJECT_DIR <- normalizePath(Sys.getenv("PROJECT_DIR", unset = "."), mustWork = TRUE)
if (!file.exists(file.path(PROJECT_DIR, "README.md")))
  stop("PROJECT_DIR does not point at the repository root: ", PROJECT_DIR)
setwd(PROJECT_DIR)

suppressPackageStartupMessages({
  library(data.table); library(GenomicRanges)
})
set.seed(1234)
N_SIM   <- 2000L
OUT_DIR <- file.path(PROJECT_DIR, "data", "chapter3")

A <- readRDS(file.path(OUT_DIR, "histone_marks_probe_peak.rds"))
B <- readRDS(file.path(OUT_DIR, "histone_marks_did.rds"))
peaks <- A$peaks; PM <- B$PM_pct; n_probe <- B$n_probe
MARKS <- A$mark_levels; STAGES <- A$stage_levels
LATE  <- B$late

# ============================ (i) POSITIVE CONTROL: proximal vs distal H3K4me3
usable  <- n_probe >= B$min_probe
prox_i  <- which(usable & !mcols(peaks)$distal)
dist_i  <- which(usable &  mcols(peaks)$distal)
message("positive control n: proximal=", length(prox_i), "  distal=", length(dist_i))

ctrl <- rbindlist(lapply(MARKS, function(mk) {
  cn <- paste(LATE, mk, sep = "|")
  vp <- PM[prox_i, cn]; vd <- PM[dist_i, cn]
  w  <- suppressWarnings(wilcox.test(vp, vd, exact = FALSE))
  data.table(mark = mk, stage = LATE,
             n_prox = length(vp), n_distal = length(vd),
             mean_pct_prox = mean(vp), mean_pct_distal = mean(vd),
             diff = mean(vp) - mean(vd), p_wilcox = unname(w$p.value))
}))
ctrl[, q_wilcox := p.adjust(p_wilcox, "BH")]
message("\n-- POSITIVE CONTROL (proximal vs distal, percentile rank, ", LATE, ") --")
print(ctrl[, .(mark, mean_pct_prox, mean_pct_distal, diff, p_wilcox, q_wilcox)])

# second control: all usable peaks vs genome-wide probe median (percentile 0.5)
ctrl2 <- rbindlist(lapply(MARKS, function(mk) {
  cn <- paste(LATE, mk, sep = "|")
  v  <- PM[usable, cn]
  w  <- suppressWarnings(wilcox.test(v, mu = 0.5, exact = FALSE))
  data.table(mark = mk, n_peaks = length(v), mean_pct = mean(v),
             p_vs_genome = unname(w$p.value))
}))
message("\n-- PEAKS vs GENOME-WIDE PROBE DISTRIBUTION (mu = 0.5) --")
print(ctrl2)

# ============================ (ii) ABSOLUTE CONTRAST at the ChIP stage
# Better powered than the DiD (one difference, not two) and it is the contrast
# that maps directly onto "Motif2 marks active distal enhancers".
idx <- B$idx; m2 <- B$m2
abs_tab <- rbindlist(lapply(MARKS, function(mk) rbindlist(lapply(STAGES, function(st) {
  cn <- paste(st, mk, sep = "|"); if (!cn %in% colnames(PM)) return(NULL)
  v <- PM[idx, cn]
  w <- suppressWarnings(wilcox.test(v[m2], v[!m2], exact = FALSE))
  data.table(mark = mk, stage = st, n_m2_pos = sum(m2), n_m2_neg = sum(!m2),
             mean_pct_m2_pos = mean(v[m2]), mean_pct_m2_neg = mean(v[!m2]),
             diff = mean(v[m2]) - mean(v[!m2]), p_wilcox = unname(w$p.value))
}))))
abs_tab[, q_wilcox := p.adjust(p_wilcox, "BH")]
message("\n-- ABSOLUTE: Motif2+ vs Motif2- distal peaks, all stages --")
print(abs_tab[, .(mark, stage, diff, p_wilcox, q_wilcox)])

# ================================================== (iii) POWER / MDE at distal n
# power must be simulated on RESAMPLED datasets. Adding a fixed shift to
# the observed deltas and re-testing is deterministic, it returns the shift
# required to make this dataset significant, not the rejection rate. Here each
# simulation draws a fresh dataset by resampling observed deltas with
# replacement (preserving the empirical distribution and the achieved n), then
# applies the shift to the Motif2+ group.
shifts <- seq(0, 0.30, by = 0.01)

power_one <- function(mk) {
  ce <- paste(B$early, mk, sep = "|"); cl <- paste(B$late, mk, sep = "|")
  d  <- PM[idx, cl] - PM[idx, ce]
  n1 <- sum(m2); n2 <- sum(!m2); pool <- d
  sapply(shifts, function(s) {
    mean(replicate(N_SIM, {
      g1 <- sample(pool, n1, replace = TRUE) + s
      g2 <- sample(pool, n2, replace = TRUE)
      suppressWarnings(wilcox.test(g1, g2, exact = FALSE)$p.value) < 0.05
    }))
  })
}
pw <- rbindlist(lapply(MARKS, function(mk)
  data.table(mark = mk, shift = shifts, power = power_one(mk))))
message("\n-- type I error check (power at shift = 0, should be ~0.05) --")
print(pw[shift == 0])

mde <- pw[power >= 0.80, .(mde_80 = min(shift)), by = mark]
mde_all <- data.table(mark = MARKS)[mde, on = "mark"]
mde_all <- merge(data.table(mark = MARKS), mde, by = "mark", all.x = TRUE)
message("\n-- MINIMUM DETECTABLE DiD at 80% power (percentile-rank units) --")
print(mde_all)

obs <- B$res[scale == "percentile", .(mark, observed_did = did)]
comp <- merge(mde_all, obs, by = "mark")
comp[, ratio_mde_over_obs := mde_80 / abs(observed_did)]
message("\n-- observed DiD vs MDE --")
print(comp)

# ============================== (iv) WHAT SAMPLE SIZE WOULD SUFFICE
# The binding constraint is the number of distal peaks, which is a property of
# the peak set, not of the mark assay: switching from ChIP-chip to ChIP-seq
# reduces measurement noise per peak but cannot increase n. The only route to
# adequate n here is to pool the distal peaks of all four heterodimers, which
# makes the extension a Chapter 4 analysis rather than a Chapter 3 one.
ann_files <- list.files(file.path(PROJECT_DIR, "data/chapter3"),
                        "_joint_annotated\\.csv$", full.names = TRUE)
pool <- rbindlist(lapply(ann_files, function(f) {
  x <- fread(f, showProgress = FALSE)
  data.table(set = sub("_joint_annotated.csv$", "", basename(f)),
             n_total = nrow(x),
             n_distal = sum(abs(x$distanceToTSS) > 1000, na.rm = TRUE),
             n_distal_m2 = sum(abs(x$distanceToTSS) > 1000 &
                               x$Motif2 %in% c(TRUE, "TRUE", 1), na.rm = TRUE))
}))
retention  <- length(idx) / pool[set == "Trh_Tgo", n_distal]   # probe-coverage rate
n_pooled   <- round(sum(pool$n_distal) * retention)
message("\n-- pooled four-heterodimer distal peaks --")
print(pool)
message("pooled distal = ", sum(pool$n_distal), " | expected usable = ", n_pooled,
        " (retention ", sprintf("%.2f", retention), ")")

# power at pooled n, using the observed H3K4me1 delta distribution
d_k4m1 <- PM[idx, paste(B$late, "H3K4Me1", sep = "|")] -
          PM[idx, paste(B$early, "H3K4Me1", sep = "|")]
pw_at <- function(n1, n2, s) mean(replicate(N_SIM, {
  g1 <- sample(d_k4m1, n1, replace = TRUE) + s
  g2 <- sample(d_k4m1, n2, replace = TRUE)
  suppressWarnings(wilcox.test(g1, g2, exact = FALSE)$p.value) < 0.05
}))
h <- round(n_pooled / 2)
sh_p  <- seq(0.005, 0.08, by = 0.005)
p_pool <- vapply(sh_p, function(s) pw_at(h, n_pooled - h, s), 0)
mde_pooled <- min(sh_p[p_pool >= 0.80])
obs_k4m1   <- abs(comp[mark == "H3K4Me1", observed_did])
pool_power_at_obs <- pw_at(h, n_pooled - h, obs_k4m1)
message("MDE(80%) pooled = ", sprintf("%.3f", mde_pooled),
        "  | power at observed ", sprintf("%.3f", obs_k4m1), " = ",
        sprintf("%.2f", pool_power_at_obs))

pooled <- list(pool = pool, retention = retention, n_pooled = n_pooled,
               mde_pooled = mde_pooled, power_at_obs = pool_power_at_obs,
               obs_k4m1 = obs_k4m1)

saveRDS(list(ctrl = ctrl, ctrl2 = ctrl2, abs_tab = abs_tab, power = pw, mde = mde_all,
             comp = comp, pooled = pooled, n_sim = N_SIM, shifts = shifts),
        file.path(OUT_DIR, "histone_marks_control_power.rds"))
fwrite(pool, file.path(PROJECT_DIR, "data/chapter3/histone_marks_pooled_scope.csv"))
fwrite(ctrl,  file.path(PROJECT_DIR, "data/chapter3/histone_marks_control_power.csv"))
fwrite(comp,  file.path(PROJECT_DIR, "data/chapter3/histone_marks_mde.csv"))
fwrite(abs_tab, file.path(PROJECT_DIR, "data/chapter3/histone_marks_absolute.csv"))
message("\nwrote histone_marks_control_power.rds + 4 csv")
