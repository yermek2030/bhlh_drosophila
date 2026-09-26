## ============================================================
## Scans Motif2 (FIMO, fimo_enhancers_genomewide, chromosomal coords) in each
## bHLH-PAS subunit's own peak set (Trh, Tgo, Sim, Dys, Sima) plus the
## Trh/Tgo joint set, same window/threshold as controls, to see if the joint
## set's Motif2 frequency (614 -> 169) reflects bHLH-PAS binding or the
## intersection. Windows must lie fully inside an enhancer (in_enh_fully),
## matched by midpoint; percentages are not comparable to ch4_ctrl_motif2_panel.R.
## Outputs: data/chapter4/table4_singlefactor_motif2.csv,
## results/figures/fig4_15_singlefactor_motif2.png, results/reported_values.tsv
## ============================================================

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
  ## Stars printed on PAL_DIV marks (Fig 4.9) use black #000000, worst-case WCAG contrast 4.34 across the ramp (vs 1.07 white, 2.61 outline grey, 3.60 #1A1A1A).
  ## Not valid for PAL_SEQ (viridis, L* 16-90): no single ink stays legible, so labels sit off the dot instead.
  ## Exception: geom_tile heatmaps use ink_on_seq() (below), switching ink by measured tile lightness, worst-case 4.16:1.
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

suppressPackageStartupMessages({
  library(GenomicRanges); library(GenomeInfoDb); library(Biostrings)
  library(BSgenome.Dmelanogaster.UCSC.dm6); library(ggplot2)
})

HELP <- file.path(PROJECT_DIR, "04_chapter_analyses")

good_chr <- c("chr2L", "chr2R", "chr3L", "chr3R", "chr4", "chrX")
MOTIF_ID <- "STRCACTGARARAAA"
HALF     <- 100L                      # 201 bp window
GENOME   <- BSgenome.Dmelanogaster.UCSC.dm6

## ---- 1. enhancer catalogue: the frame FIMO actually scanned ---------------
enh <- read.table(file.path(PROJECT_DIR, "data/reference_gene_lists",
                            "merged_embryo_5-20_enhancers_dm6.bed"),
                  sep = "\t", header = FALSE, stringsAsFactors = FALSE)
enh <- GRanges(enh[[1]], IRanges(enh[[2]] + 1L, enh[[3]]))
enh <- keepSeqlevels(enh, good_chr, pruning.mode = "coarse")
stopifnot(length(enh) == 35798L)      # == num-sequences in fimo.xml

## ---- 2. Motif2 hits in chromosomal coordinates ---------------------------
fimo <- read.delim(file.path(PROJECT_DIR, "data", "chapter3",
                             "fimo_enhancers_genomewide", "fimo.tsv"),
                   sep = "\t", header = TRUE, comment.char = "#",
                   stringsAsFactors = FALSE)
fimo <- fimo[!is.na(fimo$motif_id) & fimo$motif_id == MOTIF_ID, , drop = FALSE]
stopifnot(nrow(fimo) == 15044L)
mid  <- (fimo$start + fimo$stop) %/% 2L
hits <- GRanges(fimo$sequence_name, IRanges(mid, mid))
hits <- keepSeqlevels(hits, good_chr, pruning.mode = "coarse")

## ---- 3. single-factor peak sets ------------------------------------------
load_peaks <- function(tag) {
  f <- file.path(PROJECT_DIR, "data/final_peaks",
                 sprintf("%s_IDR0.05.final.narrowPeak", tag))
  stopifnot(file.exists(f))
  d  <- read.table(f, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
  gr <- GRanges(d[[1]], IRanges(d[[2]] + 1L, d[[3]]))
  gr <- gr[seqnames(gr) %in% good_chr]
  gr <- keepSeqlevels(gr, good_chr, pruning.mode = "coarse")
  gr[!duplicated(gr)]
}
TAGS <- c("TRH", "SIM", "DYS", "SIMA", "TGO", "VVL", "DA", "TWI", "RIB")
pk   <- lapply(TAGS, load_peaks); names(pk) <- TAGS

## retained joint row, built by the canonical Ch3/Ch4 recipe
o <- findOverlaps(pk$TRH, pk$TGO, minoverlap = 1L, ignore.strand = TRUE)
jt <- pintersect(pk$TRH[queryHits(o)], pk$TGO[subjectHits(o)]); mcols(jt) <- NULL
jt <- keepSeqlevels(jt, good_chr, pruning.mode = "coarse")
jt <- jt[!duplicated(jt)]
stopifnot(length(jt) == 614L)                      # canonical, never regress

## unbound enhancers: no peak from any factor in this panel (DA removes 26
## enhancers from the unbound pool and leaves the rate unchanged at one
## decimal place, so soloUnbdPct does not regress)
bound <- reduce(unlist(GRangesList(unname(pk))), ignore.strand = TRUE)
unb   <- enh[!overlapsAny(enh, bound, ignore.strand = TRUE)]

## ---- 4. score: 201 bp midpoint window, required inside one enhancer -------
score_set <- function(gr) {
  ctr <- (start(gr) + end(gr)) %/% 2L
  win <- GRanges(seqnames(gr), IRanges(ctr - HALF, ctr + HALF))
  keep <- win %within% enh                       # fully contained: in frame
  win  <- win[keep]
  if (!length(win)) return(list(n = 0L, k = 0L, gc = NA_real_))
  k  <- sum(overlapsAny(win, hits, ignore.strand = TRUE))
  sq <- getSeq(GENOME, win)
  gc <- 100 * mean(letterFrequency(sq, "GC", as.prob = TRUE))
  list(n = length(win), k = k, gc = gc)
}

SETS <- list(`Trh/Tgo joint` = jt, `Trh alone` = pk$TRH, `Sim alone` = pk$SIM,
             `Dys alone` = pk$DYS, `Sima alone` = pk$SIMA, `Tgo alone` = pk$TGO,
             Vvl = pk$VVL, Da = pk$DA, Twi = pk$TWI, Rib = pk$RIB,
             `Unbound enhancers` = unb)
FAMILY <- c("bHLH-PAS pair", "bHLH-PAS", "bHLH-PAS", "bHLH-PAS", "bHLH-PAS",
            "bHLH-PAS (beta)", "bHLH, no PAS", "POU-HD", "BTB-ZF",
            "bHLH class I hub", "--")
sc <- lapply(SETS, score_set)

## ---- 5. Fisher against Trh ALONE (the like-for-like reference) ------------
ref <- sc[["Trh alone"]]
tab <- data.frame(
  Set = names(SETS), TF_family = FAMILY,
  N_scanned = vapply(sc, `[[`, integer(1), "n"),
  N_Motif2_pos = vapply(sc, `[[`, integer(1), "k"),
  Mean_GC_pct = round(vapply(sc, `[[`, numeric(1), "gc"), 1),
  stringsAsFactors = FALSE, row.names = NULL)
tab$Pct_Motif2 <- round(100 * tab$N_Motif2_pos / tab$N_scanned, 1)
ft <- Map(function(k, n) {
  if (identical(c(k, n), c(ref$k, ref$n))) return(list(or = 1, p = NA_real_))
  m <- matrix(c(k, n - k, ref$k, ref$n - ref$k), nrow = 2)
  f <- fisher.test(m); list(or = unname(f$estimate), p = f$p.value)
}, tab$N_Motif2_pos, tab$N_scanned)
tab$OR_vs_Trh_alone <- round(vapply(ft, `[[`, numeric(1), "or"), 2)
tab$p_value         <- vapply(ft, `[[`, numeric(1), "p")
tab <- tab[, c("Set","TF_family","N_scanned","N_Motif2_pos","Pct_Motif2",
               "Mean_GC_pct","OR_vs_Trh_alone","p_value")]
print(tab, row.names = FALSE)

bh <- tab$Set %in% c("Trh alone","Sim alone","Dys alone","Sima alone","Tgo alone")
stopifnot(all(tab$p_value[bh & !is.na(tab$p_value)] > 0.05))   # claim check
stopifnot(tab$p_value[tab$Set == "Twi"] < 1e-5,
          tab$p_value[tab$Set == "Rib"] < 1e-4)

dir.create(file.path(PROJECT_DIR, "data", "chapter4"),
           showWarnings = FALSE, recursive = TRUE)
write.csv(tab, file.path(PROJECT_DIR, "data", "chapter4",
                         "table4_singlefactor_motif2.csv"), row.names = FALSE)

## ---- 6. recorded values (letters only, no thousands separators) --------
g <- function(s, col) tab[[col]][tab$Set == s]
vals <- list(
  soloTrhPct   = sprintf("%.1f", g("Trh alone","Pct_Motif2")),
  soloTrhN     = format(g("Trh alone","N_scanned"), trim = TRUE),
  soloSimPct   = sprintf("%.1f", g("Sim alone","Pct_Motif2")),
  soloSimN     = format(g("Sim alone","N_scanned"), trim = TRUE),
  soloDysPct   = sprintf("%.1f", g("Dys alone","Pct_Motif2")),
  soloDysN     = format(g("Dys alone","N_scanned"), trim = TRUE),
  soloSimaPct  = sprintf("%.1f", g("Sima alone","Pct_Motif2")),
  soloSimaN    = format(g("Sima alone","N_scanned"), trim = TRUE),
  soloTgoPct   = sprintf("%.1f", g("Tgo alone","Pct_Motif2")),
  soloTgoN     = format(g("Tgo alone","N_scanned"), trim = TRUE),
  soloTwiPct   = sprintf("%.1f", g("Twi","Pct_Motif2")),
  soloVvlPct   = sprintf("%.1f", g("Vvl","Pct_Motif2")),
  soloRibPct   = sprintf("%.1f", g("Rib","Pct_Motif2")),
  soloDaPct    = sprintf("%.1f", g("Da","Pct_Motif2")),
  soloDaN      = sprintf("%d",   as.integer(g("Da","N_scanned"))),
  soloDaPos    = sprintf("%d",   as.integer(g("Da","N_Motif2_pos"))),
  soloDaGc     = sprintf("%.1f", g("Da","Mean_GC_pct")),
  soloDaOR     = sprintf("%.2f", g("Da","OR_vs_Trh_alone")),
  soloDaP      = sprintf("%.3f", g("Da","p_value")),
  soloUnbdPct  = sprintf("%.1f", g("Unbound enhancers","Pct_Motif2")),
  soloBhlhMin  = sprintf("%.1f", min(tab$Pct_Motif2[bh])),
  soloBhlhMax  = sprintf("%.1f", max(tab$Pct_Motif2[bh])),
  soloJointOR  = sprintf("%.2f", g("Trh/Tgo joint","OR_vs_Trh_alone")),
  soloJointP   = sprintf("%.2f", g("Trh/Tgo joint","p_value")),
  soloVvlP     = sprintf("%.2f", g("Vvl","p_value")),
  soloWindowBp = format(2L * HALF + 1L, trim = TRUE)
)
save_values("ch4_single_factor", vals)

## ---- 7. figure -----------------------------------------------------------
BASE <- 11; ANN <- BASE - 2
fill_for <- c(
  "Trh/Tgo joint"     = PAL_HETERO[["Trh/Tgo"]],
  "Trh alone"         = PAL_SUBUNIT[["Trh"]],
  "Sim alone"         = PAL_SUBUNIT[["Sim"]],
  "Dys alone"         = PAL_SUBUNIT[["Dys"]],
  "Sima alone"        = PAL_SUBUNIT[["Sima"]],
  "Tgo alone"         = PAL_SUBUNIT[["Tgo"]],
  "Vvl"               = PAL_INK[["emphasis"]],   # the one control that matches
  "Da"                = PAL_INK[["neutral"]],
  "Twi"               = PAL_INK[["neutral"]],
  "Rib"               = PAL_INK[["neutral"]],
  "Unbound enhancers" = PAL_INK[["neutral"]])

d <- tab
d$Set <- factor(d$Set, levels = rev(names(SETS)))
d$lab <- ifelse(is.na(d$p_value),
                ## n is a count of things, not a variable
                ## symbol, so it is set UPRIGHT. Plotmath renders a bare token
                ## upright, so dropping italic() is the whole change; p keeps
                ## its italic via thesis_p_expr().
                sprintf("n==%d~'(ref.)'", d$N_scanned),
                sprintf("n==%d*';'~~%s", d$N_scanned,
                        vapply(d$p_value, thesis_p_expr, character(1))))

## Label goes inside bar if long enough to hold it, coloured via thesis_text_col() for contrast; else outside in annotation ink.
## Fit estimated from label glyph count (plotmath stripped) at ~1.2 axis-% per glyph.
## Twi (16.8%) and Rib (14.1%) bars are too short for their exponent labels, so placed outside.
.vis <- nchar(gsub("italic\\(|\\)|\\*|~|'|\\^|==|\\{|\\}", "", d$lab))
d$.inside <- d$Pct_Motif2 >= 1.2 * .vis
d$.lab_y  <- ifelse(d$.inside, d$Pct_Motif2 - 1.2, d$Pct_Motif2 + 1.2)
d$.lab_h  <- ifelse(d$.inside, 1, 0)
d$.lab_c  <- ifelse(d$.inside,
                    thesis_text_col(unname(fill_for[as.character(d$Set)])),
                    PAL_INK[["annotation"]])
p <- ggplot(d, aes(x = Set, y = Pct_Motif2, fill = Set)) +
  geom_col(width = 0.72, colour = PAL_INK[["outline"]], linewidth = 0.25) +
  geom_text(aes(label = lab, y = .lab_y, hjust = .lab_h, colour = I(.lab_c)),
            parse = TRUE, size = ANN / .pt) +
  scale_fill_manual(values = fill_for, guide = "none") +
  ## Axis stops at 50: the tallest bar is ~45.5% and every label sits
  ## inside its bar or ends well short of 50.
  scale_y_continuous(limits = c(0, 50), expand = c(0, 0),
                     breaks = seq(0, 50, 10)) +
  coord_flip() +
  ## No panel title: labs() sets only x and y.
  labs(x = NULL, y = "Peaks carrying Motif2 (%)") +
  theme_classic(base_size = BASE) +
  theme(axis.title = element_text(size = BASE - 0.5),
        axis.text  = element_text(size = BASE - 2),
        plot.title = element_text(size = BASE, hjust = 0),
        axis.line  = element_line(colour = PAL_INK[["outline"]]))

FIG <- file.path(PROJECT_DIR, "results", "figures",
                 "fig4_15_singlefactor_motif2.png")
ggsave(FIG, plot = p, width = 7.0, height = 4.4, dpi = 300, bg = "white")
cat("figure:", FIG, "\n")
cat("csv:", file.path("data", "chapter4",
                      "table4_singlefactor_motif2.csv"), "\n")
