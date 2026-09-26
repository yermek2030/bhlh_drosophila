## ============================================================
##  ch4_ctrl_motif2_panel.R
##  Motif2 carriage at enhancer peaks of Twi, Vvl, Rib, Da; GC/width-matched
##  non-Tgo background; Trh/Tgo reference. 201 bp windows, FIMO p < 1e-4.
##  FASTA via bedtools getfasta -nameOnly; counts unique(sequence_name); uses p not q.
##  Output: data/chapter4/table4_ctrl_motif2_frequency.csv
## ============================================================

## Repository root: the environment variable PROJECT_DIR, or the working directory.
PROJECT_DIR <- normalizePath(Sys.getenv("PROJECT_DIR", unset = "."), mustWork = TRUE)
if (!file.exists(file.path(PROJECT_DIR, "README.md")))
  stop("PROJECT_DIR does not point at the repository root: ", PROJECT_DIR)

good_chr <- c("chr2L", "chr2R", "chr3L", "chr3R", "chr4", "chrX")   # dm6 major arms

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
  library(GenomicRanges)
  library(Biostrings)
  library(Rsamtools)
})

## MEME Suite binaries. Resolution order: $MEME_BIN, then PATH.
## Set MEME_BIN (shell or ~/.Renviron) if the MEME Suite is not on PATH.
MEME_BIN <- Sys.getenv("MEME_BIN", unset = "")
if (!nzchar(MEME_BIN)) MEME_BIN <- dirname(Sys.which("fimo")[[1]])
if (!nzchar(MEME_BIN) || !file.exists(file.path(MEME_BIN, "fimo")))
  stop("MEME Suite not found. Install MEME >= 5.5.5 and set MEME_BIN to its bin/ directory.")
PWM      <- file.path(PROJECT_DIR, "data/chapter3/STRCACTGARARAAA.meme")
DM6      <- file.path(PROJECT_DIR, "reference", "dm6.fa")
ENH_BED  <- file.path(PROJECT_DIR, "data/reference_gene_lists/merged_embryo_5-20_enhancers_dm6.bed")
TRH_BED  <- file.path(PROJECT_DIR, "data/chapter3/Trh_Tgo_joint_peaks_filt.bed")
CTRL_FA  <- file.path(PROJECT_DIR, "data/ctrl_tf/fimo/input_fastas")
CTRL_OUT <- file.path(PROJECT_DIR, "data/ctrl_tf/fimo")
OUT_CSV  <- file.path(PROJECT_DIR, "data/chapter4/table4_ctrl_motif2_frequency.csv")

WIN     <- 100L    # half-window: 201 bp total, matching the control sets
PTHRESH <- 1e-4

stopifnot(file.exists(PWM), file.exists(DM6), file.exists(ENH_BED), file.exists(TRH_BED))

## ---- helpers ----

read_bed <- function(p) {
  d <- read.delim(p, header = FALSE, stringsAsFactors = FALSE)
  GRanges(d[[1]], IRanges(d[[2]] + 1L, d[[3]]))
}

centre_window <- function(gr, win = WIN) {
  mid <- (start(gr) + end(gr)) %/% 2L
  GRanges(seqnames(gr), IRanges(mid - win, mid + win))
}

## Count distinct sequences carrying >=1 hit below the p threshold.
## Check: sequence_name must name a region, not a chromosome.
count_fimo <- function(dir, p_thresh = PTHRESH) {
  f <- file.path(dir, "fimo.tsv")
  if (!file.exists(f)) stop("FIMO output not found: ", f)
  d <- read.delim(f, comment.char = "#", stringsAsFactors = FALSE)
  pc <- if ("p_value" %in% names(d)) "p_value" else "p.value"
  d  <- d[nzchar(d$sequence_name) & !is.na(d[[pc]]), ]
  if (all(grepl("^chr[0-9XYLRUHet]+$", unique(d$sequence_name)))) {
    stop("FIMO sequence_name values are chromosome names, not region ids: ", f,
         "\n  Regenerate the FASTA with `bedtools getfasta -nameOnly`.")
  }
  length(unique(d$sequence_name[d[[pc]] < p_thresh]))
}

run_fimo <- function(fa, outdir) {
  system2(file.path(MEME_BIN, "fimo"),
          c("--thresh", format(PTHRESH, scientific = TRUE),
            "--oc", shQuote(outdir), shQuote(PWM), shQuote(fa)),
          stdout = FALSE, stderr = FALSE)
  outdir
}

n_records <- function(fa) sum(grepl("^>", readLines(fa, warn = FALSE)))

## ---- 1. Trh/Tgo reference, rescanned at the control window ----

enh <- read_bed(ENH_BED)
trh <- read_bed(TRH_BED)
stopifnot(length(trh) == 614L)

trh_enh <- trh[countOverlaps(trh, enh) > 0L]
stopifnot(length(trh_enh) == 169L)     # matches \trhEnhN / \nTrhEnhPeaks

trh_win <- centre_window(trh_enh)
seqs <- scanFa(FaFile(DM6),
               GRanges(seqnames(trh_win), ranges(trh_win)))
names(seqs) <- paste0("trhenh_", seq_along(seqs))
TRH_FA <- file.path(CTRL_FA, "TRH_enh_summits_200bp.fa")
writeXStringSet(seqs, TRH_FA)
stopifnot(all(width(seqs) == 2L * WIN + 1L))

trh_dir <- run_fimo(TRH_FA, file.path(CTRL_OUT, "TRH_motif2"))

## ---- 2. all five rows, identical window and threshold ----

sets <- list(
  list(key = "Trh/Tgo (ref)",      fam = "bHLH-PAS",
       fa = TRH_FA,                                    dir = trh_dir),
  list(key = "TWI",                fam = "bHLH (no PAS)",
       fa = file.path(CTRL_FA, "TWI_enh_summits_200bp.fa"),
       dir = file.path(CTRL_OUT, "TWI_motif2")),
  list(key = "VVL",                fam = "POU-HD",
       fa = file.path(CTRL_FA, "VVL_enh_summits_200bp.fa"),
       dir = file.path(CTRL_OUT, "VVL_motif2")),
  list(key = "RIB",                fam = "BTB-ZF",
       fa = file.path(CTRL_FA, "RIB_enh_summits_200bp.fa"),
       dir = file.path(CTRL_OUT, "RIB_motif2")),
## Da (ENCODE ENCSR416MJQ, IDR 0.05, 923 peaks), class I bHLH hub, control for a different obligate partner.
## Included here as peak overlapping an annotated enhancer, 201 bp window centred on peak, same PWM and threshold.
## n=264, versus 197 in the single-factor figure, which requires the whole window inside one enhancer.
  list(key = "DA",                 fam = "bHLH class I hub",
       fa = file.path(CTRL_FA, "DA_enh_summits_200bp.fa"),
       dir = file.path(CTRL_OUT, "DA_motif2")),
  list(key = "Matched non-Tgo BG", fam = "--",
       fa = file.path(CTRL_FA, "matched_bg_200bp.fa"),
       dir = file.path(CTRL_OUT, "matched_bg_motif2"))
)

## Scan any set whose FIMO output is missing, so that the panel is
## reproducible from the FASTAs alone rather than requiring the control
## directories to already exist.
for (s in sets) {
  if (!file.exists(file.path(s$dir, "fimo.tsv"))) {
    message("scanning ", basename(s$fa), " -> ", basename(s$dir))
    run_fimo(s$fa, s$dir)
  }
}
n_scan <- vapply(sets, function(s) n_records(s$fa),    integer(1))
n_pos  <- vapply(sets, function(s) count_fimo(s$dir),  integer(1))

## Denominators must match the peak-set sizes reported elsewhere.
stopifnot(identical(n_scan, c(169L, 1892L, 258L, 138L, 264L, 614L)))

ref_pos <- n_pos[1]; ref_tot <- n_scan[1]

fisher_vs_ref <- function(pos, tot) {
  ft <- fisher.test(matrix(c(pos, tot - pos,
                             ref_pos, ref_tot - ref_pos),
                           nrow = 2, byrow = TRUE))
  list(or = round(unname(ft$estimate), 2), p = ft$p.value)
}
ft <- lapply(seq_along(sets)[-1], function(i) fisher_vs_ref(n_pos[i], n_scan[i]))

summary_tbl <- data.frame(
  Factor       = vapply(sets, `[[`, character(1), "key"),
  Family       = vapply(sets, `[[`, character(1), "fam"),
  N_scanned    = n_scan,
  N_Motif2_pos = n_pos,
  Pct_Motif2   = round(100 * n_pos / n_scan, 1),
  OR_vs_TrhTgo = c(1.00, vapply(ft, `[[`, numeric(1), "or")),
  p_value      = c(NA,   vapply(ft, `[[`, numeric(1), "p")),
  Window_bp    = 2L * WIN + 1L,
  Threshold    = "p < 1e-4",
  stringsAsFactors = FALSE
)

dir.create(dirname(OUT_CSV), recursive = TRUE, showWarnings = FALSE)
write.csv(summary_tbl, OUT_CSV, row.names = FALSE)
print(summary_tbl, row.names = FALSE)
cat(sprintf("\nwritten: %s\n", OUT_CSV))

## ---- 3. reported values (section ch4_control_panel) ----
## Percentages are formatted with sprintf("%.1f"), not round(), so a value
## like 18.0 keeps its trailing zero and matches its siblings in the table.


fmt_p <- function(p, digits = 2) {
  if (is.na(p)) return(NA_character_)
  s  <- formatC(p, format = "e", digits = digits)
  pr <- strsplit(s, "e")[[1]]
  paste0(pr[1], "e", as.integer(pr[2]))
}
g <- function(k, col) summary_tbl[[col]][summary_tbl$Factor == k]

## Peak-calling yields: IDR 0.05 peaks on the canonical chromosomes.
n_final_peaks <- function(f) {
  x <- read.table(file.path(PROJECT_DIR, "data", "final_peaks", f),
                  sep = "\t", header = FALSE, stringsAsFactors = FALSE)
  as.character(sum(x$V1 %in% good_chr))
}
carry <- list(nTwiPeaks = n_final_peaks("TWI_IDR0.05.final.narrowPeak"),
              nVvlPeaks = n_final_peaks("VVL_IDR0.05.final.narrowPeak"),
              nRibPeaks = n_final_peaks("RIB_IDR0.05.final.narrowPeak"))

ratio <- g("Trh/Tgo (ref)", "Pct_Motif2") / summary_tbl$Pct_Motif2[-1]

ctrl_vals <- c(carry, list(
  nTwiEnhPeaks      = g("TWI","N_scanned"),
  nTwiMotifTwo      = g("TWI","N_Motif2_pos"),
  pctTwiMotifTwo    = sprintf("%.1f", g("TWI","Pct_Motif2")),
  twiMotifTwoOR     = sprintf("%.2f", g("TWI","OR_vs_TrhTgo")),
  twiMotifTwoP      = fmt_p(g("TWI","p_value")),
  nVvlEnhPeaks      = g("VVL","N_scanned"),
  nVvlMotifTwo      = g("VVL","N_Motif2_pos"),
  pctVvlMotifTwo    = sprintf("%.1f", g("VVL","Pct_Motif2")),
  vvlMotifTwoOR     = sprintf("%.2f", g("VVL","OR_vs_TrhTgo")),
  vvlMotifTwoP      = fmt_p(g("VVL","p_value")),
  ## Da: same five values as the other controls
  nDaPeaks          = n_final_peaks("DA_IDR0.05.final.narrowPeak"),
  nDaEnhPeaks       = g("DA","N_scanned"),
  nDaMotifTwo       = g("DA","N_Motif2_pos"),
  pctDaMotifTwo     = sprintf("%.1f", g("DA","Pct_Motif2")),
  daMotifTwoOR      = sprintf("%.2f", g("DA","OR_vs_TrhTgo")),
  daMotifTwoP       = fmt_p(g("DA","p_value")),
  nRibEnhPeaks      = g("RIB","N_scanned"),
  nRibMotifTwo      = g("RIB","N_Motif2_pos"),
  pctRibMotifTwo    = sprintf("%.1f", g("RIB","Pct_Motif2")),
  ribMotifTwoOR     = sprintf("%.2f", g("RIB","OR_vs_TrhTgo")),
  ribMotifTwoP      = fmt_p(g("RIB","p_value")),
  nBgScanned        = g("Matched non-Tgo BG","N_scanned"),
  nBgMotifTwo       = g("Matched non-Tgo BG","N_Motif2_pos"),
  pctBgMotifTwo     = sprintf("%.1f", g("Matched non-Tgo BG","Pct_Motif2")),
  bgMotifTwoOR      = sprintf("%.2f", g("Matched non-Tgo BG","OR_vs_TrhTgo")),
  bgMotifTwoP       = fmt_p(g("Matched non-Tgo BG","p_value")),
  nTrhEnhPeaks      = g("Trh/Tgo (ref)","N_scanned"),
  nTrhMotifTwoEnh   = g("Trh/Tgo (ref)","N_Motif2_pos"),
  pctTrhMotifTwoEnh = sprintf("%.1f", g("Trh/Tgo (ref)","Pct_Motif2")),
  ctrlScanWindowBp  = 2L * WIN + 1L,
  ctrlFoldMin       = sprintf("%.1f", min(ratio)),
  ctrlFoldMax       = sprintf("%.1f", max(ratio))
))

save_values("ch4_control_panel", ctrl_vals)
