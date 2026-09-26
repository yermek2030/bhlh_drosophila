#!/usr/bin/env Rscript
# ============================================================
# Read every audit output and print one decision table.
# Run last. Nothing here recomputes; it only aggregates and applies the rules.
# ============================================================
## Repository root: the environment variable PROJECT_DIR, or the working directory.
PROJECT_DIR <- normalizePath(Sys.getenv("PROJECT_DIR", unset = "."), mustWork = TRUE)
if (!file.exists(file.path(PROJECT_DIR, "README.md")))
  stop("PROJECT_DIR does not point at the repository root: ", PROJECT_DIR)
setwd(PROJECT_DIR)

OUT <- file.path(PROJECT_DIR, "data", "motif2_concat_audit")
rd  <- function(f) if (file.exists(file.path(OUT, f)))
  read.csv(file.path(OUT, f), stringsAsFactors = FALSE) else NULL

rows <- list()
add  <- function(test, question, result, verdict)
  rows[[length(rows) + 1L]] <<- data.frame(test, question, result, verdict,
                                           stringsAsFactors = FALSE)

# ---- information content --------------------------------------------------
ic <- rd("T0_ic_profile.csv")
if (!is.null(ic)) {
  tail_ic <- ic$IC_dm6bg[ic$block == "tail_10_15"]
  n_low   <- sum(tail_ic < 0.3)
  add("T0 information content",
      "Are tail positions 10-15 below ~0.3 bits (i.e. filler)?",
      sprintf("tail total %.2f bits over 6 positions (core 1-9: %.2f bits); %d/6 tail positions < 0.3 bits",
              sum(tail_ic), sum(ic$IC_dm6bg[ic$block == "core_1_9"]), n_low),
      if (n_low >= 4) "TAIL IS FILLER - trim the PWM to the core."
      else "TAIL IS NOT LOW-IC FILLER - IC alone does not support trimming.")
}

# ---- background-order sweep ----------------------------------------------
sw <- file.path(OUT, "streme_order_sweep", "T1_order_sweep_summary.txt")
if (file.exists(sw)) {
  x  <- readLines(sw); x <- x[!grepl("^DIMER", x) & nzchar(trimws(x))]
  x  <- x[grepl("^(Trh|Sim|Dys|Sima) ", x)]
  cat("\n--- T1 raw ---\n"); writeLines(x)
  f  <- read.table(text = x, col.names = c("dimer", "order", "rank", "cons", "E"),
                   stringsAsFactors = FALSE)
  hi   <- f[f$order %in% c(3, 4), ]
  lost <- sum(hi$cons == "NOT_RECOVERED")
  # Checks for a purine 3' tail after CACTGA in the recovered consensus.
  # Character class {A, G, R}; excludes filler codes W, D, H, V, N, M.
  # All sixteen runs give a 6 bp tail (RARAAA / GAGAAA / GARAAA / AARAAA / AAAAAA).
  tail_after <- function(s) {
    i <- regexpr("CACTGA", s, fixed = TRUE)
    if (i < 0) return(NA_integer_)
    nchar(sub("[^AGR].*$", "", substring(s, i + 6L)))
  }
  f$tail_bp <- vapply(f$cons, tail_after, 0L)
  cat("\ntail length after CACTGA, by order:\n")
  print(f[, c("dimer", "order", "cons", "tail_bp")], row.names = FALSE)
  hi <- f[f$order %in% c(3, 4) & !is.na(f$tail_bp), ]
  add("T1 background order",
      "Does Motif2 (and its 3' tail) survive an order-3/4 background model?",
      sprintf("%d of %d order-3/4 runs lost the motif; in the %d that recovered it the purine tail after CACTGA is %s",
              lost, sum(f$order %in% c(3, 4)), nrow(hi),
              if (min(hi$tail_bp) == max(hi$tail_bp))
                sprintf("%d bp in every run", hi$tail_bp[1])
              else sprintf("%s bp (median %.1f)",
                           paste(range(hi$tail_bp), collapse = "-"),
                           median(hi$tail_bp))),
      if (lost > 0) "LOST under higher-order background - composition artefact. STOP AND RETHINK."
      else if (median(hi$tail_bp) >= 5)
        "SURVIVES with the tail intact - the tail is not a zero-order composition artefact."
      else "MOTIF SURVIVES BUT TAIL SHRINKS - report the core, describe the tail as flanking composition.")
}

# ---- half-scan and register ----------------------------------------------
mp <- rd("T2_min_attainable_p.csv")
if (!is.null(mp)) {
  add("T2a threshold feasibility",
      "Can each sub-PWM even be scanned at p < 1e-4?",
      paste(sprintf("%s: min p = %.2e", mp$pwm, mp$minp), collapse = "; "),
      "The 6-bp tail CANNOT reach p<1e-4. Any 'tail-only: 0 hits' result at that threshold is a construction artefact, not evidence.")
}
sh <- rd("T2_register_sharpness.csv")
if (!is.null(sh)) {
  cat("\n--- T2 mean-profile sharpness (DESCRIPTIVE ONLY) ---\n")
  print(sh, row.names = FALSE)
  add("T2b mean-profile register (descriptive)",
      "Mean tail log-odds at register 0, by set.",
      sprintf("in-peak %.2f (n=%d) vs background %.2f (n=%d)",
              sh$score_at_register0[sh$set == "in_joint_peaks"],
              sh$n[sh$set == "in_joint_peaks"],
              sh$score_at_register0[sh$set == "genome_background"],
              sh$n[sh$set == "genome_background"]),
      paste("DO NOT COMPARE THE Z COLUMNS: n differs 139-fold, so the background",
            "flank SD is far smaller and its Z is inflated by smoothness, not effect.",
            "Use T2c."))
}

# ---- T2c site-level register (the comparable test) ---------------------------
sl <- rd("T2c_register_sitelevel.csv")
if (!is.null(sl)) {
  cat("\n--- T2c site-level register ---\n"); print(sl, row.names = FALSE)
  add("T2c register specificity (site level)",
      "Per-site: is the tail at register 0 specifically enriched in peaks?",
      sprintf("mean delta over own flank: in-peak %.2f vs background %.2f (%.1fx, Mann-Whitney p = %.3g); tail-positive OR = %.2f (p = %.3g); off-register control OR = %.2f (p = %.3g)",
              sl$delta_inpeak, sl$delta_background, sl$delta_ratio, sl$wilcox_p,
              sl$tailpos_OR, sl$tailpos_p, sl$offreg_OR, sl$offreg_p),
      if (sl$tailpos_OR > 2 && sl$tailpos_p < 0.05 && sl$tailpos_OR > 3 * sl$offreg_OR)
        "COMPOSITE ELEMENT - the register is real, position-specific and peak-specific."
      else if (sl$tailpos_OR > 2 && sl$offreg_OR > 2)
        "COMPOSITION ARTEFACT - in-peak sequence is A-rich at every offset, not just register 0."
      else "NO REGISTER STRUCTURE - the 15-bp call looks like a concatenation.")
}

# ---- does the claim survive on the core ----------------------------------
or <- rd("T3_core_only_enhancer_OR.csv")
if (!is.null(or)) {
  r <- or[or$enh_width_filter == "ge15_matches_published", ]
  cat("\n--- T3 core-only enhancer OR ---\n"); print(or, row.names = FALSE)
  add("T3a enhancer enrichment on core",
      "Does OR = 1.28, p = 3.5e-3 survive with the 9-bp core alone?",
      sprintf("core-only OR = %.2f [%.2f-%.2f], p = %.3g (%d/%d peaks vs %d/%d enhancers)",
              r$OR, r$CI_lo, r$CI_hi, r$p, r$joint_pos, r$joint_tot,
              r$enh_pos, r$enh_tot),
      if (r$OR >= 1.15 && r$p < 0.05)
        "SURVIVES - report a 9-bp core; the tail is a flanking compositional feature."
      else "DOES NOT SURVIVE - the enrichment needs the full 15-bp form; keep it and defend it.")
}
tt <- rd("T3_tomtom_4way_core_pairwise.csv")
if (!is.null(tt)) {
  tl <- rd("T3_core_trim_log.csv")
  w  <- if (!is.null(tl)) nchar(tl$core_consensus[1]) else NA_integer_
  add("T3b pan-heterodimer equivalence on core",
      "Does TOMTOM equivalence survive with core-trimmed PWMs?",
      sprintf("%d pairwise comparisons on a %d-bp common core window, worst q = %.3g",
              nrow(tt), w, max(tt$qval)),
      if (max(tt$qval) <= 1e-4)
        "SURVIVES - pan-heterodimer identity does not depend on the tail."
      else "WEAKENS on the core alone - state the dependence explicitly.")
}

verdict <- do.call(rbind, rows)
cat("\n================= MOTIF2 CONCATENATION AUDIT: VERDICT =================\n")
for (i in seq_len(nrow(verdict))) {
  cat("\n[", verdict$test[i], "]\n  Q: ", verdict$question[i],
      "\n  R: ", verdict$result[i], "\n  -> ", verdict$verdict[i], "\n", sep = "")
}
write.csv(verdict, file.path(OUT, "T4_verdict.csv"), row.names = FALSE)
cat("\nWrote", file.path(OUT, "T4_verdict.csv"), "\n")

