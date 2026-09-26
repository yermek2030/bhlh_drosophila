# ============================================================
# histone_marks_4_summary.R
# Stage D: appendix figure and reported values for histone-mark exclusion.
# Reads : histone_marks_probe_peak.rds, histone_marks_did.rds, histone_marks_control_power.rds
# Writes: results/figures/figA_histone_marks_power.png, results/reported_values.tsv (section ch3_histone_marks)
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

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork)
})

OUT_DIR <- file.path(PROJECT_DIR, "data", "chapter3")
FIG_DIR <- file.path(PROJECT_DIR, "results", "figures")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

A <- readRDS(file.path(OUT_DIR, "histone_marks_probe_peak.rds"))
B <- readRDS(file.path(OUT_DIR, "histone_marks_did.rds"))
C <- readRDS(file.path(OUT_DIR, "histone_marks_control_power.rds"))

MARK_LAB <- c(H3K4Me1 = "H3K4me1", H3K27Ac = "H3K27ac",
              H3K4Me3 = "H3K4me3", H3K27Me3 = "H3K27me3")
# CVD-safe: blue/orange opposition, never red/green
PAL <- c("promoter-proximal" = "#0072B2", "distal" = "#E69F00")

base_sz <- 8
thm <- theme_classic(base_size = base_sz) +
  theme(axis.ticks.length = unit(2, "pt"),
        legend.position = "none",
        plot.title = element_text(size = base_sz, hjust = 0, face = "plain"),
        strip.background = element_blank(),
        strip.text = element_text(size = base_sz),
        axis.text = element_text(size = base_sz - 2),
        plot.margin = margin(4, 8, 4, 4))

# ------------------------------------------------ panel a: positive control
ctrl <- as.data.table(C$ctrl)
pa_dat <- rbind(
  ctrl[, .(mark, set = "promoter-proximal", pct = mean_pct_prox, n = n_prox)],
  ctrl[, .(mark, set = "distal",            pct = mean_pct_distal, n = n_distal)]
)
pa_dat[, mark_lab := factor(MARK_LAB[mark], levels = unname(MARK_LAB))]
pa_dat[, set := factor(set, levels = names(PAL))]

# q-values in a fixed row above all bars (avoids per-bar label collisions)
sig_lab <- ctrl[, .(mark, lab = sprintf("q = %.0e", q_wilcox))]
sig_lab[, mark_lab := factor(MARK_LAB[mark], levels = unname(MARK_LAB))]

pa <- ggplot(pa_dat, aes(mark_lab, pct, fill = set)) +
  geom_col(position = position_dodge(width = 0.72), width = 0.62) +
  geom_hline(yintercept = 0.5, linetype = 2, linewidth = 0.3, colour = "grey45") +
  geom_text(data = sig_lab, aes(mark_lab, 1.055, label = lab), inherit.aes = FALSE,
            size = (base_sz - 2) / .pt, colour = "grey20") +
  # direct series labels in the whitespace above bar pair 2
  annotate("text", x = 1.52, y = 0.965, label = "promoter-proximal",
           hjust = 0, size = (base_sz - 1) / .pt, colour = PAL[[1]], fontface = "bold") +
  annotate("text", x = 1.52, y = 0.905, label = "distal",
           hjust = 0, size = (base_sz - 1) / .pt, colour = PAL[[2]], fontface = "bold") +
  # dashed-line identity, in the right margin (clip off) so it collides with nothing
  annotate("text", x = 4.62, y = 0.50, label = "genome-wide\nprobe median",
           hjust = 0, vjust = 0.5, size = (base_sz - 2) / .pt, colour = "grey35",
           lineheight = 0.9) +
  coord_cartesian(xlim = c(0.55, 4.45), ylim = c(0, 1.09), clip = "off") +
  scale_fill_manual(values = PAL) +
  scale_y_continuous(breaks = seq(0, 1, 0.25), expand = expansion(mult = c(0, 0))) +
  labs(title = "All four marks separate promoter-proximal from distal peaks as expected",
       x = NULL, y = "Mean probe percentile rank") +
  thm + theme(plot.margin = margin(4, 62, 4, 4))

# ------------------------------------------------ panel b: power vs observed
pw <- as.data.table(C$power); comp <- as.data.table(C$comp)
pw[, mark_lab := factor(MARK_LAB[mark], levels = unname(MARK_LAB))]
comp[, mark_lab := factor(MARK_LAB[mark], levels = unname(MARK_LAB))]

# power actually achieved at each observed |DiD|, this is the quantity the
# panel title must be true of, so it is computed, not asserted.
comp[, power_at_obs := vapply(seq_len(.N), function(i) {
  p <- pw[mark == comp$mark[i]]
  approx(p$shift, p$power, xout = abs(comp$observed_did[i]))$y
}, 0)]
stopifnot(max(comp$power_at_obs) < 0.80)   # guards the panel-b claim-title
pb_title <- sprintf(
  "No mark reaches 80%% power at its observed Motif2 effect (power %.2f\u2013%.2f)",
  min(comp$power_at_obs), max(comp$power_at_obs))

pb <- ggplot(pw, aes(shift, power, group = mark_lab)) +
  geom_hline(yintercept = 0.80, linetype = 2, linewidth = 0.3, colour = "grey45") +
  geom_line(linewidth = 0.5, colour = "grey30") +
  geom_vline(data = comp, aes(xintercept = abs(observed_did)),
             linewidth = 0.45, colour = "#D55E00") +
  geom_point(data = comp, aes(x = abs(observed_did), y = power_at_obs),
             inherit.aes = FALSE, colour = "#D55E00", size = 1.3) +
  geom_text(data = comp, aes(x = abs(observed_did) + 0.022, y = power_at_obs,
                             label = sprintf("observed %.3f\npower %.2f",
                                             abs(observed_did), power_at_obs)),
            inherit.aes = FALSE, hjust = 0, vjust = 0.4,
            size = (base_sz - 2) / .pt, colour = "#D55E00", lineheight = 0.95) +
  geom_text(data = comp, aes(x = 0.295, y = 0.885,
                             label = sprintf("80%% power at %.2f", mde_80)),
            inherit.aes = FALSE, hjust = 1,
            size = (base_sz - 2) / .pt, colour = "grey20") +
  facet_wrap(~ mark_lab, nrow = 1) +
  scale_x_continuous(limits = c(0, 0.30), breaks = c(0, 0.1, 0.2, 0.3),
                     expand = expansion(mult = c(0.03, 0.05))) +
  scale_y_continuous(limits = c(0, 1.02), breaks = seq(0, 1, 0.25),
                     expand = expansion(mult = c(0.02, 0.02))) +
  labs(title = pb_title,
       x = "Simulated difference in temporal delta (percentile-rank units)",
       y = "Power at \u03b1 = 0.05") +
  thm + theme(panel.spacing.x = unit(9, "pt"))

fig <- pa / pb + plot_layout(heights = c(1, 1)) +
  plot_annotation(tag_levels = "a",
                  theme = theme(plot.margin = margin(2, 2, 2, 2)),
                  tag_prefix = "", tag_suffix = "") &
  theme(plot.tag = element_text(size = base_sz + 2, face = "bold"))

FIG <- file.path(FIG_DIR, "figA_histone_marks_power.png")
ggsave(FIG, plot = fig, width = 190, height = 125, units = "mm", dpi = 300)
message("wrote ", FIG)

# ---------------------------------------------------------- reported values

res_pct <- as.data.table(B$res)[scale == "percentile"]
k4me3   <- ctrl[mark == "H3K4Me3"]
fmt_p   <- function(p) if (p < 1e-4) sprintf("%.1fe%d",
                    p / 10^floor(log10(p)), floor(log10(p))) else sprintf("%.3f", p)

vals <- list(
  histoneNStages        = length(A$stage_levels),
  histoneNMarks         = length(A$mark_levels),
  histoneNProbes        = format(A$n_common, big.mark = ""),
  histoneNProbesDmSix   = format(unname(A$lift_stats["dm6"]), big.mark = ""),
  histoneLiftPct        = sprintf("%.2f", 100 * A$lift_stats["dm6"] / A$lift_stats["dm2"]),
  histonePeaksWithProbe = length(unique(S4Vectors::queryHits(A$hits))),
  histoneNDistalUsable  = length(B$idx),
  histoneNDistalMTwoPos = sum(B$m2),
  histoneNDistalMTwoNeg = sum(!B$m2),
  histoneNProxCtrl      = k4me3$n_prox,
  histoneCtrlKFourMeThreeProx  = sprintf("%.2f", k4me3$mean_pct_prox),
  histoneCtrlKFourMeThreeDist  = sprintf("%.2f", k4me3$mean_pct_distal),
  histoneCtrlKFourMeThreeQ     = fmt_p(k4me3$q_wilcox),
  # min(q) is the best q across the four marks, not the worst. For a null
  # claim the defensible statement is "even the most significant mark reached
  # only q = X", so the value recorded in results/reported_values.tsv is named Best, not Worst.
  histoneDidBestQ       = sprintf("%.2f", min(res_pct$q_wilcox)),
  histoneDidBestPEmp    = sprintf("%.2f", min(res_pct$p_empirical)),
  histoneMdeMin         = sprintf("%.2f", min(comp$mde_80)),
  histoneMdeMax         = sprintf("%.2f", max(comp$mde_80)),
  histonePowerMin       = sprintf("%.2f", min(comp$power_at_obs)),
  histonePowerMax       = sprintf("%.2f", max(comp$power_at_obs)),
  histoneNPerm          = format(B$n_perm, big.mark = ""),
  histoneNSim           = format(C$n_sim, big.mark = ""),
  # pooled four-heterodimer scope (the route to adequate n; Chapter 4 territory)
  histonePooledDistal   = sum(C$pooled$pool$n_distal),
  histonePooledUsable   = C$pooled$n_pooled,
  histonePooledMde      = sprintf("%.3f", C$pooled$mde_pooled),
  histonePooledPower    = sprintf("%.2f", C$pooled$power_at_obs)
)
save_values("ch3_histone_marks", vals)
print(as.data.table(list(name = names(vals), value = unlist(lapply(vals, as.character)))))
