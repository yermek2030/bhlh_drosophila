#!/usr/bin/env Rscript
# ============================================================
# Step 5, record the concatenation-audit values in
# results/reported_values.tsv (safe to re-run; existing values are replaced).
# Run after steps 00-04 have completed.
# ============================================================
## Repository root: the environment variable PROJECT_DIR, or the working directory.
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
setwd(PROJECT_DIR)


OUT <- file.path(PROJECT_DIR, "data", "motif2_concat_audit")
rd  <- function(f) {
  p <- file.path(OUT, f)
  if (!file.exists(p)) stop("missing audit output: ", f, " -- run 00-04 first.")
  read.csv(p, stringsAsFactors = FALSE)
}

ic   <- rd("T0_ic_profile.csv")
mp   <- rd("T2_min_attainable_p.csv")
sl   <- rd("T2c_register_sitelevel.csv")
pc   <- rd("T2c_panelC_values.csv")
or   <- rd("T3_core_only_enhancer_OR.csv")
tt   <- rd("T3_tomtom_4way_core_pairwise.csv")
trim <- rd("T3_core_trim_log.csv")

# --- STREME order sweep: how many order-3/4 runs kept the motif --------------
sw <- file.path(OUT, "streme_order_sweep", "T1_order_sweep_summary.txt")
sweep <- read.table(text = grep("^(Trh|Sim|Dys|Sima) ",
                                readLines(sw), value = TRUE),
                    col.names = c("dimer", "order", "rank", "cons", "E"),
                    stringsAsFactors = FALSE)
hi   <- sweep[sweep$order %in% c(3, 4), ]
kept <- sum(hi$cons != "NOT_RECOVERED")

# Purine tail length after CACTGA, counting the characters {A, G, R}. The tail
# is the same length in every run, so a single value is reported.
tail_after <- function(s) {
  i <- regexpr("CACTGA", s, fixed = TRUE)
  if (i < 0) return(NA_integer_)
  nchar(sub("[^AGR].*$", "", substring(s, i + 6L)))
}
tails <- vapply(sweep$cons[sweep$cons != "NOT_RECOVERED"], tail_after, 0L)
stopifnot(length(unique(tails)) == 1L)   # invariant across all 16 runs
tail_bp <- unname(tails[1])

# Scientific notation, e.g. 2.1e-12
fmt_p <- function(x, d = 1) {
  if (x >= 0.01) return(formatC(x, format = "f", digits = 3))
  e <- floor(log10(x)); m <- round(x / 10^e, d)
  sprintf("%se%d", formatC(m, format = "f", digits = d), e)
}

core_ic <- sum(ic$IC_dm6bg[ic$block == "core_1_9"])
tail_ic <- sum(ic$IC_dm6bg[ic$block == "tail_10_15"])
r15     <- or[or$enh_width_filter == "ge15_matches_published", ]
reg0    <- pc[pc$offset == "register 0", ]
off5    <- pc[pc$offset != "register 0", ]

vals <- list(
  # information content
  motifTwoCoreIC          = sprintf("%.1f", core_ic),
  motifTwoTailIC          = sprintf("%.1f", tail_ic),
  motifTwoTailLowICpos    = as.character(sum(ic$IC_dm6bg[ic$block == "tail_10_15"] < 0.3)),
  # background-order sweep
  motifTwoOrderKept       = as.character(kept),
  motifTwoOrderTested     = as.character(nrow(hi)),
  motifTwoOrderBestE      = fmt_p(min(suppressWarnings(as.numeric(sweep$E)), na.rm = TRUE)),
  motifTwoOrderTailBp     = as.character(tail_bp),
  # T2a threshold feasibility
  motifTwoTailMinP        = fmt_p(mp$minp[mp$pwm == "tail_6"]),
  motifTwoCoreMinP        = fmt_p(mp$minp[mp$pwm == "core_9"]),
  motifTwoFullMinP        = fmt_p(mp$minp[mp$pwm == "full_15"]),
  # T2c site-level register
  motifTwoRegSiteIn       = sprintf("%.1f", sl$delta_inpeak),
  motifTwoRegSiteBg       = sprintf("%.1f", sl$delta_background),
  motifTwoRegSiteN        = format(sl$n_inpeak, big.mark = ""),
  motifTwoRegSiteBgN      = format(sl$n_background, big.mark = ""),
  motifTwoRegSiteP        = fmt_p(sl$wilcox_p),
  motifTwoRegSiteOR       = sprintf("%.2f", sl$tailpos_OR),
  motifTwoRegSiteORP      = fmt_p(sl$tailpos_p),
  motifTwoRegOffOR        = sprintf("%.2f", sl$offreg_OR),
  motifTwoRegOffP         = sprintf("%.2f", sl$offreg_p),
  motifTwoRegPctIn        = sprintf("%.1f", reg0$pct_in_peaks),
  motifTwoRegPctBg        = sprintf("%.1f", reg0$pct_background),
  motifTwoRegOffPctIn     = sprintf("%.1f", off5$pct_in_peaks),
  motifTwoRegOffPctBg     = sprintf("%.1f", off5$pct_background),
  # T3a core-only enhancer enrichment
  motifTwoCoreEnhOR       = sprintf("%.2f", r15$OR),
  motifTwoCoreEnhORloCI   = sprintf("%.2f", r15$CI_lo),
  motifTwoCoreEnhORhiCI   = sprintf("%.2f", r15$CI_hi),
  motifTwoCoreEnhP        = fmt_p(r15$p),
  motifTwoCoreJointPct    = sprintf("%.1f", r15$joint_pct),
  motifTwoCoreEnhPct      = sprintf("%.1f", r15$enh_pct),
  motifTwoCoreJointN      = as.character(r15$joint_pos),
  # T3b pan-heterodimer equivalence on the core
  motifTwoCoreWorstQ      = fmt_p(max(tt$qval)),
  motifTwoCoreNComp       = as.character(nrow(tt)),
  motifTwoCoreWindow      = as.character(nchar(trim$core_consensus[1]))
)

save_values("concat_audit", vals)
for (nm in names(vals)) cat(sprintf("  %-24s %s\n", nm, vals[[nm]]))
