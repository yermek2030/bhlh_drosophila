#!/usr/bin/env Rscript
# ============================================================
# ch3_fig_cme_stratification.R
# Three panels: (A) CME positive/negative split of peak set, (B) both-strand vs forward-strand-only recovery, (C) observed 5-mer rate vs chance expectation for interval length.
# Recomputes counts from peak set and checks against text values; stops on mismatch.
# Writes results/figures/fig3_14_cme_stratification.png
# ============================================================

## Repository root: the environment variable PROJECT_DIR, or the working directory.
PROJECT_DIR <- normalizePath(Sys.getenv("PROJECT_DIR", unset = "."), mustWork = TRUE)
if (!file.exists(file.path(PROJECT_DIR, "README.md")))
  stop("PROJECT_DIR does not point at the repository root: ", PROJECT_DIR)

VALUES_FILE <- file.path(PROJECT_DIR, "results", "reported_values.tsv")

read_value <- function(name) {
  tab <- read.delim(VALUES_FILE, colClasses = "character", quote = "", comment.char = "")
  tab$value[match(name, tab$name)]
}


## Colours
PAL_INK <- c(
  outline    = "#333333",   # bar and Venn outlines
  gridline   = "#DDDDDD",
  annotation = "#666666",   # p-values, n = counts
  threshold  = "#882255",   # significance thresholds, +/-1.96 lines
  fit        = "#4477AA",   # fitted / OLS trend lines
  neutral    = "#BBBBBB",   # non-significant, background sets
  band       = "#EDEDED",   # shaded stage/facet bands and strip backgrounds
  highlight_band = "#4477AA", # tint marking an agreement region (used at low alpha)
  emphasis   = "#CC3311",   # the one highlighted subset in an otherwise neutral plot
## Ink for stars on PAL_DIV marks (Fig 4.9): #000000 chosen, worst-case WCAG 4.34 across the ramp; #1A1A1A (3.60) also passes.
## Not valid for PAL_SEQ (viridis): no single ink stays legible, worst-case <=1.43; nudge labels off the dot instead.
## Exception: labelled geom_tile heatmaps use ink_on_seq() (below), switching ink per tile by measured lightness, worst-case 4.16:1.
  on_mark    = "#000000",

  ## Two-tier ink for ink_on_seq(). `on_mark_dark` is #1A1A1A rather
  ## than the #000000 above: on the LIGHT half of viridis both clear the floor
  ## comfortably, and #1A1A1A matches the ink weight of the surrounding axis text,
  ## so the matrix does not look heavier than the panel around it.
  on_mark_dark  = "#1A1A1A",
  on_mark_light = "#FFFFFF"
)

PAL_BINARY <- c(pos = "#CC6677", neg = "#4477AA")

pal_binary <- function(pos_label, neg_label) {
  setNames(c(PAL_BINARY[["pos"]], PAL_BINARY[["neg"]]), c(pos_label, neg_label))
}


## Bar chart
thesis_text_col <- function(fill) {
  rgb <- grDevices::col2rgb(fill) / 255
  lum <- 0.299 * rgb[1, ] + 0.587 * rgb[2, ] + 0.114 * rgb[3, ]
  ifelse(lum > 0.5, "black", "white")
}

thesis_p_expr <- function(p, digits = 2) {
  if (is.na(p)) return(NA_character_)
  ## The formatted value is QUOTED, so plotmath prints it as formatted rather
  ## than re-formatting it as a number: unquoted, a Fisher p of 1.000 renders
  ## as "1" and 2.10 as "2.1", losing the fixed significant figures that
  ## must be preserved at every mention.
  if (p < 1e-3) {
    parts <- strsplit(formatC(p, format = "e", digits = digits), "e")[[1]]
    sprintf("italic(p)=='%s' %%*%% 10^{%d}", parts[1], as.integer(parts[2]))
  } else {
    ## `digits` governs both branches. "fg" formats to significant digits
    ## rather than decimal places; flag "#" keeps trailing zeros (p = 1
    ## prints 1.0, not 1).
    sprintf("italic(p)=='%s'",
            formatC(signif(p, digits), format = "fg", digits = digits,
                    flag = "#"))
  }
}

thesis_bar <- function(data, x, y, fill = x,
                        value = NULL, value_fmt = "%.1f%%",
                        count = NULL, sig = NULL,
                        facet = NULL, n_lab = NULL, n_lab_size = 3.0,
                        fills = NULL, y_max = NULL, y_lab = NULL,
                        dashed_at = NULL, global_p = NULL, caption_sig = NULL,
                        inside_frac = 0.28, bar_width = 0.62, base_size = 11,
                        title = NULL, facet_free_x = TRUE) {

  d <- as.data.frame(data)
  d$.x <- d[[x]]; d$.y <- d[[y]]; d$.fill <- d[[fill]]

  ## primary value label
  if (!is.null(value)) {
    d$.value <- as.character(d[[value]])
  } else if (!is.na(value_fmt)) {
    d$.value <- sprintf(value_fmt, d$.y)
  } else {
    d$.value <- NA_character_
  }
  d$.count <- if (!is.null(count)) as.character(d[[count]]) else NA_character_
  d$.sig   <- if (!is.null(sig))   as.character(d[[sig]])   else NA_character_

  ## y range + headroom (space for the star that always sits above the bar)
  y_top_data <- max(d$.y, na.rm = TRUE)
  if (is.null(y_max)) y_max <- y_top_data * 1.32
  gap  <- 0.045 * y_max            # one label-line height, scaled to axis
  inside_cut <- inside_frac * y_max

  ## fill -> contrast colour for inside text
  if (!is.null(fills)) {
    missing_fill <- setdiff(unique(as.character(d$.fill)), names(fills))
    if (length(missing_fill))
      stop("thesis_bar(): fill level(s) not found in `fills` (would render grey): ",
           paste(sprintf("'%s'", missing_fill), collapse = ", "),
           "\n  Available in `fills`: ",
           paste(sprintf("'%s'", names(fills)), collapse = ", "))
  }
  fill_levels <- if (!is.null(fills)) names(fills) else unique(as.character(d$.fill))
  fill_cols   <- if (!is.null(fills)) fills else scales::hue_pal()(length(fill_levels))
  names(fill_cols) <- fill_levels
  d$.txtcol <- thesis_text_col(fill_cols[as.character(d$.fill)])

  ## ---- per-bar label geometry ------------------------------------------
  d$.inside <- d$.y >= inside_cut
  has_count <- !is.na(d$.count)
  # value line: inside -> just below bar top; above -> primary value uppermost
  d$.val_y   <- ifelse(d$.inside, d$.y - gap,
                       ifelse(has_count, d$.y + 2 * gap, d$.y + gap))
  d$.val_col <- ifelse(d$.inside, d$.txtcol, "black")
  # count line: inside -> below the value; above -> nearer the bar than value
  d$.cnt_y   <- ifelse(d$.inside, d$.y - 2 * gap, d$.y + gap)
  d$.cnt_col <- ifelse(d$.inside, d$.txtcol, "black")
  # significance always above the bar; clears above-bar labels when present
  above_labels <- (!d$.inside) & (has_count | !is.na(d$.value))
  d$.sig_y    <- d$.y + ifelse(above_labels, 3 * gap, gap * 1.4)

  ## ---- build plot -------------------------------------------------------
  p <- ggplot(d, aes(x = .x, y = .y, fill = .fill)) +
    geom_col(width = bar_width, colour = "black", linewidth = 0.4,
             show.legend = FALSE)

  if (!is.null(dashed_at))
    p <- p + geom_hline(yintercept = dashed_at, linetype = "dashed",
                        colour = PAL_INK[["neutral"]], linewidth = 0.5)

  ## value label
  if (any(!is.na(d$.value)))
    p <- p + geom_text(aes(y = .val_y, label = .value, colour = .val_col),
                       fontface = "bold", size = 3.5, show.legend = FALSE,
                       na.rm = TRUE)
  ## count label
  if (any(!is.na(d$.count)))
    p <- p + geom_text(aes(y = .cnt_y, label = .count, colour = .cnt_col),
                       size = 3.0, show.legend = FALSE, na.rm = TRUE)
  ## per-bar significance (above, spaced)
  if (any(!is.na(d$.sig)))
    p <- p + geom_text(aes(y = .sig_y, label = .sig),
                       fontface = "italic", size = 3.4, colour = PAL_INK[["annotation"]],
                       na.rm = TRUE)

  p <- p +
    scale_colour_identity() +
    scale_fill_manual(values = fill_cols) +
    scale_y_continuous(limits = c(0, y_max), expand = expansion(mult = c(0, 0)))

  ## facets, free_x so each panel shows only its own x categories
  if (!is.null(facet))
    p <- p + facet_wrap(stats::as.formula(paste0("~", facet)), nrow = 1,
                        scales = if (facet_free_x) "free_x" else "fixed")

  ## per-facet n = ... annotation, top-right
  if (!is.null(n_lab)) {
    key <- if (!is.null(facet)) facet else NULL
    nd  <- if (!is.null(key)) d[!duplicated(d[[key]]), ] else d[1, , drop = FALSE]
    nd$.nlab <- as.character(nd[[n_lab]])
    parse_n  <- any(grepl("==|~", nd$.nlab))
    p <- p + geom_text(data = nd,
                       aes(x = Inf, y = Inf, label = .nlab),
                       hjust = 1.08, vjust = 1.6, size = n_lab_size, colour = PAL_INK[["annotation"]],
                       fontface = "italic", parse = parse_n,
                       inherit.aes = FALSE)
  }

  ## single centred p-value, drawn once (as a plot tag), not per-facet,
  ## so a faceted figure shows one central p-value instead of one per panel.
  plot_tag <- NULL
  if (!is.null(global_p)) plot_tag <- parse(text = thesis_p_expr(global_p))[[1]]

  cap <- caption_sig
  p <- p +
    labs(x = NULL, y = y_lab, title = title, caption = cap, tag = plot_tag) +
    theme_classic(base_size = base_size) +
    theme(
      ## BOLD: the panel letter is the one bold element of the house style.
      plot.tag          = element_text(size = base_size + 2, face = "bold"),
      plot.tag.position  = c(0.5, 0.86),
      strip.text       = element_text(face = "bold", size = base_size),
      strip.background = element_rect(fill = NA, colour = "black", linewidth = 0.6),
      axis.title.y     = element_text(size = base_size - 0.5),
      axis.text.x      = element_text(size = base_size - 2, colour = "black"),
      axis.text.y      = element_text(size = base_size - 2, colour = "black"),
      axis.ticks.x     = element_blank(),
      panel.spacing    = unit(1.0, "cm"),
      plot.caption     = element_text(hjust = 0, size = base_size - 2.5,
                                      colour = PAL_INK[["outline"]]),
      plot.title       = element_text(size = base_size, face = "bold"),
      plot.margin      = margin(8, 10, 8, 8)
    )
  p
}

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(Biostrings)
  library(BSgenome.Dmelanogaster.UCSC.dm6)
  library(ggplot2)
  library(patchwork)
})


OUTDIR <- file.path(PROJECT_DIR, "results", "figures")
## ---- 1. values reported in the text ------------------------------------------
## jointPeaks, nCmePos, nCmeNeg and pctCmePos are recorded in
## results/reported_values.tsv by chapter3_trh_tgo_binding.Rmd. The forward-only
## count and the per-position CME probability used for the chance curve are the
## values stated in the Chapter 3 text.
V <- setNames(lapply(c("jointPeaks", "nCmePos", "nCmeNeg", "pctCmePos"),
                     read_value),
              c("jointPeaks", "nCmePos", "nCmeNeg", "pctCmePos"))
if (any(is.na(unlist(V))))
  stop("Run chapter3_trh_tgo_binding.Rmd first: values missing from ",
       VALUES_FILE)
TEXT_VALUES <- list(nCmePosFwdOnly = 115L, cmePfloorObserved = 8.6e-4)

n_total  <- as.integer(V$jointPeaks)
n_pos    <- as.integer(V$nCmePos)
n_neg    <- as.integer(V$nCmeNeg)
n_fwd    <- as.integer(TEXT_VALUES$nCmePosFwdOnly)
pct_pos  <- as.numeric(V$pctCmePos)
p_pos    <- TEXT_VALUES$cmePfloorObserved

## ---- 2. re-derive the split from the peak set and check the values --------
jp <- readRDS(file.path(PROJECT_DIR, "data", "chapter3",
                        "joint_peaks_filt.rds"))
stopifnot(length(jp) == n_total)

seqs   <- getSeq(BSgenome.Dmelanogaster.UCSC.dm6, jp)
hit_f  <- elementNROWS(vmatchPattern("ACGTG", seqs)) > 0   # forward strand
hit_r  <- elementNROWS(vmatchPattern("CACGT", seqs)) > 0   # reverse complement
obs_pos <- sum(hit_f | hit_r)
obs_fwd <- sum(hit_f)

stopifnot(obs_pos == n_pos, obs_fwd == n_fwd,
          n_total - obs_pos == n_neg,
          abs(100 * obs_pos / n_total - pct_pos) < 0.05)

peak_w  <- width(jp)
mean_w  <- mean(peak_w)

## Panel A: CME-positive vs CME-negative split, using shared bar helper (matches Fig 3.6, Fig 3.12).
## Colours: PAL_BINARY via pal_binary() for CME-positive/negative bars, PAL_INK$neutral for discarded forward-only and chance curve/ribbon.
## Observed carriage line uses PAL_BINARY$pos, same as CME-positive bars.
BASE <- 11; ANN <- BASE - 2; TICK <- BASE - 2

col_pos <- PAL_BINARY[["pos"]]
col_neg <- PAL_INK[["neutral"]]
col_ann <- PAL_INK[["annotation"]]
col_out <- PAL_INK[["outline"]]

## panel C keeps a bare theme_classic: it is a curve, not a bar chart, so it
## takes the axis typography of the bar helper without its bar machinery.
thm <- theme_classic(base_size = BASE) +
  theme(axis.title.x = element_text(size = BASE - 0.5),
        axis.title.y = element_text(size = BASE - 0.5),
        axis.text    = element_text(size = TICK, colour = "black"),
        plot.tag     = element_text(size = BASE + 2, face = "bold"),
        legend.position = "none")

dA <- data.frame(
  cls = factor(c("CME-positive", "CME-negative"),
               levels = c("CME-positive", "CME-negative")),
  n   = c(n_pos, n_neg))
dA$Pct   <- 100 * dA$n / n_total
dA$Count <- sprintf("n = %d", dA$n)

pA <- thesis_bar(dA, x = "cls", y = "Pct", fill = "cls",
                 value_fmt = "%.1f%%", count = "Count",
                 fills = pal_binary("CME-positive", "CME-negative"),
                 y_max = 82, base_size = BASE,
                 y_lab = "Percentage of Trh/Tgo peaks (%)") +
  theme(plot.tag.position = "topleft")

## ---- 4. panel B: recovery by scan convention -----------------------------
## Both bars count CME-positive peaks, so the both-strand bar keeps the
## CME-positive ink of panel A; the forward-only scan is the convention this
## dissertation does not use and recedes to the neutral background ink.
dB <- data.frame(
  scan = factor(c("Both strands", "Forward only"),
                levels = c("Both strands", "Forward only")),
  n    = c(n_pos, n_fwd))
dB$Value <- sprintf("%d", dB$n)
## the forward-only bar carries what it costs; the reference bar needs no label
dB$Count <- c(NA, sprintf("%d missed", n_pos - n_fwd))

pB <- thesis_bar(dB, x = "scan", y = "n", fill = "scan",
                 value = "Value", count = "Count",
                 fills = c("Both strands" = col_pos,
                           "Forward only" = col_neg),
                 y_max = 230, base_size = BASE,
                 y_lab = "CME-positive peaks recovered") +
  theme(plot.tag.position = "topleft")

## ---- 5. panel C: chance carriage vs interval width -----------------------
## A 5-bp word at per-position probability p, scanned on both strands of an
## interval of width w, is expected in 1 - (1 - p)^(2(w - 4)) of intervals.
wgrid <- seq(min(peak_w), max(peak_w) + 40, by = 5)
dC <- data.frame(w = wgrid,
                 chance = 1 - (1 - p_pos)^(2 * (wgrid - 4)))
chance_at_mean <- 1 - (1 - p_pos)^(2 * (mean_w - 4))
pct_chance     <- 100 * chance_at_mean
deficit        <- pct_chance - pct_pos

## width at which the chance expectation overtakes the observed rate
w_cross <- 4 + log(1 - pct_pos / 100) / (2 * log(1 - p_pos))

rug_df <- data.frame(w = peak_w)
## the shortfall itself, drawn: chance above, observed below
rib_df  <- dC[100 * dC$chance > pct_pos, ]

pC <- ggplot(dC, aes(w, 100 * chance)) +
  ## 1. the gap between chance and observed, the panel's whole argument
  geom_ribbon(data = rib_df, aes(x = w, ymin = pct_pos, ymax = 100 * chance),
              inherit.aes = FALSE, fill = col_neg, alpha = 0.22) +
  ## 2. chance expectation
  geom_line(colour = col_neg, linewidth = 1.2) +
  ## 3. observed carriage, in the CME-positive ink of panel A
  geom_hline(yintercept = pct_pos, colour = col_pos,
             linetype = "22", linewidth = 0.8) +
  ## 4. the mean peak width, tying the width distribution to the chance level
  annotate("segment", x = mean_w, xend = mean_w, y = 0, yend = pct_chance,
           colour = col_out, linetype = "dotted", linewidth = 0.4) +
  annotate("segment", x = mean_w, xend = mean_w, y = pct_pos, yend = pct_chance,
           colour = col_out, linewidth = 0.5,
           arrow = arrow(ends = "both", type = "closed",
                         length = unit(0.055, "in"))) +
  annotate("point", x = mean_w, y = pct_chance, colour = col_out, size = 2.1) +
  ## 5. observed peak widths
  geom_rug(data = rug_df, aes(x = w), inherit.aes = FALSE, sides = "b",
           alpha = 0.25, colour = col_out, length = unit(0.035, "npc")) +
  ## short labels, each next to the mark it names
  annotate("text", x = max(wgrid), y = 100 * max(dC$chance) - 4.5,
           label = "chance expectation", hjust = 1, vjust = 0,
           size = ANN / .pt, colour = col_ann) +
  annotate("text", x = max(wgrid), y = pct_pos - 4.5,
           label = sprintf("observed %.1f%%", pct_pos), hjust = 1, vjust = 0.5,
           size = ANN / .pt, colour = col_pos, fontface = "bold") +
  annotate("text", x = mean_w + 22, y = (pct_pos + pct_chance) / 2,
           label = sprintf("%.1f points\nbelow chance", deficit),
           hjust = 0, vjust = 0.5, size = ANN / .pt, colour = col_ann,
           lineheight = 0.95) +
  annotate("text", x = mean_w + 22, y = pct_chance + 3.5,
           label = sprintf("chance %.1f%% at the mean peak width (%.0f bp)",
                           pct_chance, mean_w),
           hjust = 0, vjust = 0, size = ANN / .pt, colour = col_ann) +
  annotate("text", x = min(wgrid), y = 7.5, label = "peak widths",
           hjust = 0, vjust = 0.5, size = ANN / .pt, colour = col_ann) +
  scale_y_continuous(limits = c(0, 100), expand = expansion(mult = c(0, 0.02))) +
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.03))) +
  labs(x = "Interval width (bp)",
       y = "Intervals carrying the word (%)") +
  thm

## ---- 6. assemble ---------------------------------------------------------
## No in-panel titles are used: the caption names every panel.
## tag_levels "A": panel letters are capitals throughout.
fig <- (pA | pB) / pC +
  plot_layout(heights = c(1, 1.25)) +
  plot_annotation(tag_levels = "A")

png_path <- file.path(OUTDIR, "fig3_14_cme_stratification.png")
ggsave(png_path, plot = fig, width = 8.6, height = 6.9, dpi = 300, bg = "white")

cat(sprintf("PNG: %s\n", png_path))
cat(sprintf("  CME-positive %d / CME-negative %d of %d (%.1f%%)\n",
            n_pos, n_neg, n_total, pct_pos))
cat(sprintf("  forward-only recovery %d of %d (%.1f%%)\n",
            n_fwd, n_pos, 100 * n_fwd / n_pos))
cat(sprintf("  chance carriage at mean width %.0f bp: %.1f%% (p = %.2e)\n",
            mean_w, 100 * chance_at_mean, p_pos))
