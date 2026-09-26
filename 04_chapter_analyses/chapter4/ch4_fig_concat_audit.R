#!/usr/bin/env Rscript
# ============================================================
# ch4_fig_concat_audit.R: tail composition, register, dispensability for chapter 4.
# (A) STREME sweep, background orders 0,2,3,4, four heterodimers, 16 runs: motif2 recovered in all; order 0 not favoured.
# (B) Genome-wide 9-bp core scan: A-rich tail enriched at register 0 in-peak vs background; 5 bp off-register shows no separation.
# Reads cached outputs, checked against 02_motif_analysis/concat_audit/05_reported_values.R; writes results/figures/ch4_concat_audit.png.
# ============================================================

## Repository root: the environment variable PROJECT_DIR, or the working directory.
PROJECT_DIR <- normalizePath(Sys.getenv("PROJECT_DIR", unset = "."), mustWork = TRUE)
if (!file.exists(file.path(PROJECT_DIR, "README.md")))
  stop("PROJECT_DIR does not point at the repository root: ", PROJECT_DIR)

VALUES_FILE <- file.path(PROJECT_DIR, "results", "reported_values.tsv")


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
  ## Ink for stars printed on PAL_DIV marks in Ch4 4x2 matrix (Fig 4.9); distinct from `annotation` (white panel). Black #000000, worst-case WCAG 4.34 across PAL_DIV; #1A1A1A (3.60) also passes.
  ## Not valid for PAL_SEQ (viridis): no single ink stays legible (best 1.43); nudge the label off the dot instead.
  ## Exception: geom_tile heatmaps use ink_on_seq() (below) to switch ink per tile by measured lightness, worst-case 4.16 across viridis; do not hand-roll with a value threshold.
  on_mark    = "#000000",

  ## Two-tier ink for ink_on_seq(). `on_mark_dark` is #1A1A1A rather
  ## than the #000000 above: on the LIGHT half of viridis both clear the floor
  ## comfortably, and #1A1A1A matches the ink weight of the surrounding axis text,
  ## so the matrix does not look heavier than the panel around it.
  on_mark_dark  = "#1A1A1A",
  on_mark_light = "#FFFFFF"
)

PAL_HETERO <- c(
  "Trh/Tgo"  = "#117733",   # green    L* 44
  "Sim/Tgo"  = "#E69F00",   # amber    L* 71
  "Dys/Tgo"  = "#66CCEE",   # sky      L* 77
  "Sima/Tgo" = "#AA3377"    # magenta  L* 42
)

PAL_SUBUNIT <- c(
  "Trh"  = PAL_HETERO[["Trh/Tgo"]],
  "Sim"  = PAL_HETERO[["Sim/Tgo"]],
  "Dys"  = PAL_HETERO[["Dys/Tgo"]],
  "Sima" = PAL_HETERO[["Sima/Tgo"]],
  "Tgo"  = "#DDCC77"
)

PAL_ARCH <- c(
  alpha_generic = "#4477AA",  # generic alpha-subunit (all four)
  tgo           = PAL_SUBUNIT[["Tgo"]],
## CME/E-box duplex colour: #9BA090, worst-case CIE76 dE 24.0 across normal/deutan/protan/tritan, over four panels.
## Constraints: |dL*| >= 15 from Tgo (L* 82 vs 65), chroma <= 14; picked from 113 
  cme_duplex    = "#9BA090",  # CME / E-box DNA duplex (sage, neutral scaffold)
  motif_two     = "#882255",  # Motif2 element, where drawn
  ink           = "#333333",
  highlight     = "lightblue"   # soft yellow behind the E-box label only
)


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

suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })
## panel B uses the bar chart helper shared with Fig 3.6, 3.12 and 5.6

CA     <- file.path(PROJECT_DIR, "data", "motif2_concat_audit")
OUT    <- file.path(PROJECT_DIR, "results", "figures")
## ---- recorded values (results/reported_values.tsv) ---------------------------
load_values <- function(path = VALUES_FILE) {
  stopifnot(file.exists(path))
  tab <- read.delim(path, colClasses = "character", quote = "", comment.char = "")
  as.list(setNames(tab$value, tab$name))
}
V <- load_values()

## ---- plot theme and layout settings -------------------------------------
BASE <- 11; ANN <- BASE - 2
INK  <- PAL_INK[["outline"]];   ANNC <- PAL_INK[["annotation"]]
NEUT <- PAL_INK[["neutral"]];   M2   <- PAL_ARCH[["motif_two"]]
th_house <- function(base = BASE) {
  theme_classic(base_size = base) +
    theme(axis.title   = element_text(size = base - 0.5),
          axis.text    = element_text(size = base - 2, colour = "black"),
          ## black at theme_classic weight, so panel A's frame matches panel B,
          ## which is drawn by thesis_bar() in the Fig 3.6 / 3.12 / 5.6 style
          axis.line    = element_line(linewidth = 0.5, colour = "black"),
          axis.ticks   = element_line(linewidth = 0.5, colour = "black"),
          panel.grid   = element_blank(),
          plot.title   = element_text(size = base, face = "plain"),
          plot.tag     = element_text(size = base + 2, face = "bold"),
          legend.title = element_blank(),
          legend.text  = element_text(size = base - 2),
          legend.key.height = unit(10, "pt"))
}

## ---- panel A: background-order sweep ------------------------------------
sw_raw <- readLines(file.path(CA, "streme_order_sweep",
                              "T1_order_sweep_summary.txt"), warn = FALSE)
sw <- read.table(text = grep("^(Trh|Sim|Dys|Sima)\\s", sw_raw, value = TRUE),
                 col.names = c("dimer", "order", "rank", "consensus", "E"),
                 stringsAsFactors = FALSE)
sw$E <- as.numeric(sw$E)
sw$hetero <- factor(paste0(sw$dimer, "/Tgo"),
                    levels = c("Trh/Tgo", "Sim/Tgo", "Dys/Tgo", "Sima/Tgo"))
sw$negloge <- -log10(sw$E)

## every run must be CACTGA-anchored (the claim the panel rests on), and the
## high-order arm must match the corresponding values recorded in results/reported_values.tsv
stopifnot(nrow(sw) == 16L,
          all(grepl("CACTGA", sw$consensus)),
          sum(sw$order %in% c(3, 4)) == as.integer(V$motifTwoOrderTested),
          sum(sw$order %in% c(3, 4) & grepl("CACTGA", sw$consensus)) ==
            as.integer(V$motifTwoOrderKept))

pA <- ggplot(sw, aes(factor(order), negloge, colour = hetero, group = hetero)) +
  ## E = 1 is where a run stops being significant at all, and without it a
  ## reader cannot tell which end of a -log10 axis is the good end (Sim/Tgo at
  ## order 4 sits BELOW it, at E = 2). Same role as the AUROC panel's
  ## "0.5 = no discrimination" rule.
  geom_hline(yintercept = 0, linetype = 2, linewidth = 0.4, colour = ANNC) +
  ## x is the LEVEL, not a number: the x scale is discrete, and a numeric
  ## position is rejected ("Discrete value supplied to a continuous scale").
  ## Below the rule, not above it: above, the label would sit in the order-0
  ## column where Trh/Tgo and Sima/Tgo both plot at about 0.7.
  annotate("text", x = "0", y = 0, hjust = 0, vjust = 1.5,
           label = "italic(E)==1", parse = TRUE, size = (ANN - 1) / .pt,
           colour = ANNC) +
  geom_line(linewidth = 0.6, alpha = 0.85) +
  geom_point(aes(shape = rank == 1), size = 2.6) +
  scale_colour_manual(values = PAL_HETERO[levels(sw$hetero)]) +
  ## breaks puts rank 1 first: mapped on a logical, the default order is
  ## FALSE-then-TRUE, so the legend read "rank > 1" before "rank 1"
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1),
                     breaks = c(TRUE, FALSE),
                     labels = c(`TRUE` = "rank 1", `FALSE` = "rank > 1")) +
  ## a discrete scale pads 0.6 of a category at each end by default, which left
  ## a wide gutter between the y axis and the order-0 tick and again past
  ## order~4; halved to 0.3
  scale_x_discrete(expand = expansion(add = 0.3)) +
  scale_y_continuous(expand = expansion(mult = c(0.06, 0.14))) +
## No panel title: caption's (A) sentence already covers it.
## Y axis: quantity with transform in parens; one line only, 38 chars, fits 2.5 in of 3.4 in canvas.
## X label "results rank": rank among motifs STREME returned for that run, broken after STREME for legend width.
  labs(x = "STREME background model order",
       y = expression(Motif2 ~ "significance (" * -log[10] ~ italic(E) * ")"),
       colour = "Heterodimer", shape = "STREME\nresults rank") +
  ## colour key first: colour is the primary key the lines are matched to,
  ## so its legend is ordered before the shape legend.
  guides(colour = guide_legend(order = 1), shape = guide_legend(order = 2)) +
  th_house() +
  ## the local th_house() blanks legend titles; these two legends need theirs,
  ## because "rank 1" and a colour swatch do not say what they are keyed to.
  ## The legend also sat a full 11 pt box-spacing out from the panel plus its
  ## own margin, pulled in to 2 pt, with the keys tightened.
  theme(legend.position = "right",
        legend.title = element_text(size = BASE - 2, face = "bold"),
        ## 1 pt left the "STREME results rank" title sitting on the Sima/Tgo key
        ## above it; the two legends need a visible gap between them
        legend.spacing.y = unit(12, "pt"),
        legend.box.spacing = unit(2, "pt"),
        legend.margin = margin(0, 0, 0, 0),
        legend.key.width = unit(12, "pt"),
        plot.margin = margin(8, 4, 8, 8))

## ---- panel B: register specificity at genomic core occurrences ----------
pc <- read.csv(file.path(CA, "T2c_panelC_values.csv"), stringsAsFactors = FALSE)
sl <- read.csv(file.path(CA, "T2c_register_sitelevel.csv"), stringsAsFactors = FALSE)
stopifnot(nrow(pc) == 2L,
          abs(pc$pct_in_peaks[1]   - as.numeric(V$motifTwoRegPctIn))    < 0.06,
          abs(pc$pct_background[1] - as.numeric(V$motifTwoRegPctBg))    < 0.06,
          abs(pc$OR[1]             - as.numeric(V$motifTwoRegSiteOR))   < 0.01,
          abs(pc$OR[2]             - as.numeric(V$motifTwoRegOffOR))    < 0.01,
          sl$n_inpeak     == as.integer(V$motifTwoRegSiteN),
          sl$n_background == as.integer(V$motifTwoRegSiteBgN))

## Bar chart uses thesis_bar(): two boxed facet strips, bars labelled
## on x axis, bold percentage with denominator beneath, italic panel
## statistic top right. Same two categories in both panels; y axis names the unit, strips give A-tail position, x axis gives peaks/not.
dB <- data.frame(
  offset = factor(rep(c("A-tail in register", "A-tail shifted 5 bp"), each = 2),
                  levels = c("A-tail in register", "A-tail shifted 5 bp")),
  set    = factor(rep(c("In peaks", "Genome-wide"), 2),
                  levels = c("In peaks", "Genome-wide")),
  pct    = c(pc$pct_in_peaks[1], pc$pct_background[1],
             pc$pct_in_peaks[2], pc$pct_background[2]),
  n      = rep(c(sl$n_inpeak, sl$n_background), 2),
  stringsAsFactors = FALSE)
## the denominators, not a derived numerator: both are asserted above against
## the recorded values motifTwoRegSiteN and motifTwoRegSiteBgN. Numbers are
## formatted without a thousands separator, as a plain %d.
dB$cnt <- sprintf("n = %d", dB$n)

## or and p per register as the per-panel annotation. atop() stacks the two
## lines; the p-value is formatted by the shared helper, so its exponent cannot
## drift from the value it formats.
dB$ann <- rep(sprintf("atop(OR~%.2f, %s)", pc$OR,
                      vapply(pc$p, thesis_p_expr, character(1))), each = 2)

pB <- thesis_bar(
  dB, x = "set", y = "pct", fill = "set",
  value_fmt = "%.1f%%", count = "cnt", facet = "offset", n_lab = "ann",
  fills = c("In peaks" = M2, "Genome-wide" = NEUT),
  ## The sites are occurrences of the 9-bp SCANNED core, not of the hexamer, so
  ## the axis names the core the way chapter4 names it at five other sites
  ## ("the CACTGA-anchored core") rather than inventing a third form.
  y_max = 40, y_lab = "CACTGA-anchored core sites with an A-tail (%)")

## One file per panel, lettered as in Fig 4.9, rather than having patchwork
## draw "A"/"B" into the pixels. Panels are stacked in one column, so each
## panel is wide.
png_A <- file.path(OUT, "fig4_14_concat_audit_A.png")
png_B <- file.path(OUT, "fig4_14_concat_audit_B.png")
ggsave(png_A, plot = pA, width = 8.6, height = 3.4, dpi = 300, bg = "white")
ggsave(png_B, plot = pB, width = 8.6, height = 3.6, dpi = 300, bg = "white")

cat(sprintf("PNG: %s\n     %s\n", png_A, png_B))
cat(sprintf("  sweep: %d runs, all CACTGA-anchored; %s of %s order-3/4 runs kept\n",
            nrow(sw), V$motifTwoOrderKept, V$motifTwoOrderTested))
cat(sprintf("  register 0: %.1f%% in-peak vs %.1f%% background, OR %.2f\n",
            pc$pct_in_peaks[1], pc$pct_background[1], pc$OR[1]))
cat(sprintf("  off register: %.1f%% vs %.1f%%, OR %.2f (p = %.2f)\n",
            pc$pct_in_peaks[2], pc$pct_background[2], pc$OR[2], pc$p[2]))
cat(sprintf("  sites: %d in peaks, %d genome-wide\n", sl$n_inpeak, sl$n_background))
