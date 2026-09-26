# ============================================================
# ch3_memechip_inventory.R
# Motif inventory for Chapter 3 Trh/Tgo MEME-ChIP run (10 MEME + 7 STREME motifs, E < 0.05 cutoff).
# Parses run directory, tabulates each motif with best database match and disposition.
# Outputs: results/tables/memechip_trh_tgo_inventory.tsv, results/reported_values.tsv
# ============================================================

## Repository root: the environment variable PROJECT_DIR, or the working directory.
PROJECT_DIR <- normalizePath(Sys.getenv("PROJECT_DIR", unset = "."), mustWork = TRUE)
if (!file.exists(file.path(PROJECT_DIR, "README.md")))
  stop("PROJECT_DIR does not point at the repository root: ", PROJECT_DIR)

RESULTS_DIR   <- file.path(PROJECT_DIR, "results")

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


RUN <- file.path(PROJECT_DIR, "data/chapter3/memechip_Trh_Tgo_masked")
stopifnot(dir.exists(RUN))

# ---- 1. MEME motifs --------------------------------------------------------
meme_txt <- readLines(file.path(RUN, "meme_out/meme.txt"), warn = FALSE)
ml <- grep("^MOTIF ", meme_txt, value = TRUE)
meme <- data.frame(
  source    = "MEME",
  rank      = as.integer(sub(".*MEME-([0-9]+)\\s.*", "\\1", ml)),
  consensus = sub("^MOTIF ([^ ]+) .*", "\\1", ml),
  width     = as.integer(sub(".*width =\\s*([0-9]+).*", "\\1", ml)),
  sites     = as.integer(sub(".*sites =\\s*([0-9]+).*", "\\1", ml)),
  evalue    = as.numeric(sub(".*E-value = ([0-9.e+-]+)\\s*$", "\\1", ml)),
  stringsAsFactors = FALSE
)
stopifnot(nrow(meme) == 10, !any(is.na(meme$evalue)))

# ---- 2. STREME motifs ------------------------------------------------------
str_txt <- readLines(file.path(RUN, "streme_out/streme.txt"), warn = FALSE)
sl <- grep("^MOTIF ", str_txt)
lp <- grep("^letter-probability matrix:", str_txt)
lp <- vapply(sl, function(i) lp[which(lp > i)[1]], numeric(1))
streme <- data.frame(
  source    = "STREME",
  rank      = as.integer(sub(".*STREME-([0-9]+)\\s*$", "\\1", str_txt[sl])),
  consensus = sub("^MOTIF [0-9]+-([^ ]+) .*", "\\1", str_txt[sl]),
  width     = as.integer(sub(".*w=\\s*([0-9]+).*", "\\1", str_txt[lp])),
  sites     = as.integer(sub(".*nsites=\\s*([0-9]+).*", "\\1", str_txt[lp])),
  evalue    = as.numeric(sub(".*E=\\s*([0-9.e+-]+)\\s*$", "\\1", str_txt[lp])),
  stringsAsFactors = FALSE
)
stopifnot(nrow(streme) == 7, !any(is.na(streme$evalue)))

# ---- 3. Which motifs MEME-ChIP itself carried forward ----------------------
smry <- read.delim(file.path(RUN, "summary.tsv"), comment.char = "#",
                   stringsAsFactors = FALSE)
smry <- smry[!is.na(smry$MOTIF_INDEX) & nzchar(smry$MOTIF_ID), ]
carried <- sub("^[0-9]+-", "", smry$CONSENSUS)

inv <- rbind(meme, streme)
inv$carried <- inv$consensus %in% carried
SIG_THRESH <- 0.05
inv$sig <- inv$evalue < SIG_THRESH
stopifnot(identical(inv$carried, inv$sig))   # pipeline filter == E < 0.05

# ---- 4. Best TOMTOM database match per motif -------------------------------
best_match <- function(f) {
  if (!file.exists(f)) return(NULL)
  tt <- read.delim(f, comment.char = "#", stringsAsFactors = FALSE)
  tt <- tt[nzchar(tt$Query_ID), ]
  if (!nrow(tt)) return(NULL)
  do.call(rbind, lapply(split(tt, tt$Query_ID), function(d)
    d[which.min(d$q.value), c("Query_ID", "Target_ID", "q.value")]))
}
tt_meme   <- best_match(file.path(RUN, "meme_tomtom_out/tomtom.tsv"))
tt_streme <- best_match(file.path(RUN, "streme_tomtom_out/tomtom.tsv"))
tt <- rbind(tt_meme, tt_streme)
tt$key <- sub("^[0-9]+-", "", tt$Query_ID)
inv$db_target <- tt$Target_ID[match(inv$consensus, tt$key)]
inv$db_q      <- tt$q.value[match(inv$consensus, tt$key)]

cat("\n---- complete inventory ----\n")
print(inv[, c("source", "rank", "consensus", "width", "sites",
              "evalue", "carried", "db_target", "db_q")], row.names = FALSE)

# ---- 5. Disposition --------------------------------------------------------
inv$disposition <- ifelse(
  !inv$sig, "Below threshold",
  ifelse(inv$consensus == "GTRCACTGARARAAA", "Retained (Motif2)",
  ifelse(inv$source == "MEME", "Excluded: low complexity",
         "Not pursued")))

# ---- 6. Number formatting: two significant figures, no thousands separator --
sf2 <- function(v) sprintf(paste0("%.", max(0, 1 - floor(log10(v))), "f"), v)
fmt_eval <- function(x) {
  vapply(x, function(v) {
    if (v >= 0.01 && v < 100) {
      sf2(v)
    } else {
      ex <- floor(log10(v)); ma <- signif(v / 10^ex, 2)
      sprintf("%.1fe%d", ma, ex)
    }
  }, character(1))
}
fmt_int <- function(x) formatC(x, format = "d", big.mark = "")

# ---- 7. Appendix table ------------------------------------------------------
out <- file.path(RESULTS_DIR, "tables", "memechip_trh_tgo_inventory.tsv")
dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
inv_out <- data.frame(tool = inv$source, rank = inv$rank, consensus = inv$consensus,
                      width = fmt_int(inv$width), sites = fmt_int(inv$sites),
                      evalue = fmt_eval(inv$evalue), disposition = inv$disposition,
                      best_db_match = inv$db_target, best_db_q = inv$db_q,
                      stringsAsFactors = FALSE)
inv_out <- inv_out[order(inv_out$tool != "MEME", inv_out$rank), ]
write.table(inv_out, out, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
cat("\nwrote:", out, "\n")

# ---- 8. Values quoted in the text -------------------------------------------
s1 <- inv[inv$source == "STREME" & inv$rank == 1, ]
vals <- list(
  nMemeMotifsReturned   = nrow(meme),
  nMemeMotifsSig        = sum(meme$evalue < SIG_THRESH),
  nMemeMotifsNonSig     = sum(meme$evalue >= SIG_THRESH),
  nStremeMotifsReturned = nrow(streme),
  nStremeMotifsSig      = sum(streme$evalue < SIG_THRESH),
  memeChipSigThresh     = SIG_THRESH,
  memeNextEval          = sf2(min(meme$evalue[meme$evalue >= SIG_THRESH])),
  stremeOneCons         = s1$consensus,
  stremeOneWidth        = fmt_int(s1$width),
  stremeOneSites        = fmt_int(s1$sites),
  stremeOneEval         = sf2(s1$evalue),
  stremeOneBestMatch    = "Mes2",
  stremeOneBestQ        = sf2(s1$db_q)
)
save_values("ch3_memechip", vals)
cat("\n---- values recorded ----\n")
print(unlist(vals))
