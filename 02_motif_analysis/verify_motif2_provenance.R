## ============================================================
## verify_motif2_provenance.R
## Checks Motif2 input facts against Methods: MEME centring length, N counts of masked FASTAs, CME sites masked, PWM site count/E-value/TOMTOM p-value, FIMO/PWM background, permutation resolution and MC SE (5000 draws).
## Output: results/reported_values.tsv section ch2_motif2_inputs
## Run: Rscript 02_motif_analysis/verify_motif2_provenance.R
## ============================================================
PROJECT_DIR <- normalizePath(Sys.getenv("PROJECT_DIR", unset = "."), mustWork = TRUE)
if (!file.exists(file.path(PROJECT_DIR, "README.md")))
  stop("PROJECT_DIR does not point at the repository root: ", PROJECT_DIR)

DATA_DIR      <- file.path(PROJECT_DIR, "data")
CH3_DIR <- file.path(DATA_DIR, "chapter3")
CH4_DIR <- file.path(DATA_DIR, "chapter4")

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
## MEME Suite binaries. Resolution order: $MEME_BIN, then PATH.
MEME_BIN <- Sys.getenv("MEME_BIN", unset = "")
if (!nzchar(MEME_BIN)) MEME_BIN <- dirname(Sys.which("fimo")[[1]])
if (!nzchar(MEME_BIN) || !file.exists(file.path(MEME_BIN, "fimo")))
  stop("MEME Suite not found. Install MEME >= 5.5.5 and set MEME_BIN to its bin/ directory.")
MC <- file.path(CH3_DIR, "memechip_Trh_Tgo_masked")
DG <- file.path(CH4_DIR, "motif2_diagnostic")
n_count <- function(f) { x <- readLines(f); x <- x[!startsWith(x, ">")]; sum(nchar(gsub("[^Nn]", "", x))) }

# MEME-ChIP centring: MEME was run on seqs-centered; read the training-set lengths
mt <- readLines(file.path(MC, "meme_out", "meme.txt"))
stopifnot(any(grepl("^command: meme .*seqs-centered", mt)))
i0 <- grep("^-------------", mt)[1]; i1 <- grep("^\\*{10}", mt); i1 <- i1[i1 > i0][1]
tok <- unlist(strsplit(trimws(mt[(i0 + 1):(i1 - 1)]), "\\s+"))
lens <- as.integer(tok[seq(3, length(tok), by = 3)])
centre_len <- max(lens)
cat(sprintf("MEME training set: %d sequences, length %d-%d\n", length(lens), min(lens), max(lens)))
stopifnot(centre_len == 100L)

raw  <- file.path(CH3_DIR, "Trh_Tgo_joint_peaks_filt_raw.fa")
hard <- file.path(CH3_DIR, "Trh_Tgo_joint_peaks_filt_hardmasked.fa")
stopifnot(identical(unname(tools::md5sum(raw)), unname(tools::md5sum(hard))), n_count(hard) == 0L)
n_rep <- n_count(file.path(CH3_DIR, "Trh_Tgo_joint_peaks_filt_repeatmasked.fa"))
n_cme <- n_count(file.path(DG, "peaks_masked.fa"))
bed   <- read.table(file.path(DG, "cme_positions.bed"), stringsAsFactors = FALSE)
stopifnot(all(grepl("^peak_", bed[[1]])))          # peak-indexed frame, not chromosomal

pwm <- readLines(file.path(CH3_DIR, "STRCACTGARARAAA.meme"))
lp  <- grep("letter-probability", pwm, value = TRUE)
pwm_sites <- as.integer(sub(".*nsites= *([0-9]+).*", "\\1", lp))
pwm_e     <- sub(".*E= *([0-9.eE+-]+).*", "\\1", lp)
tt <- system2(file.path(MEME_BIN, "tomtom"),
              c("-no-ssc", "-text", "-min-overlap", "5", "-dist", "pearson", "-thresh", "1", "-evalue",
                file.path(CH3_DIR, "STRCACTGARARAAA.meme"), file.path(MC, "meme_out", "meme.txt")),
              stdout = TRUE, stderr = FALSE)
tt <- read.delim(text = paste(tt[!startsWith(tt, "#") & nzchar(tt)], collapse = "\n"), stringsAsFactors = FALSE)
tp <- tt$p.value[tt$Target_ID == "GTRCACTGARARAAA"]; stopifnot(length(tp) == 1L)
sci <- function(x, d = 1) formatC(x, format = "e", digits = d)
perm_n <- 5000L                                     # draws per permutation test
# FIMO --nrdb-- background (as scanned) versus the background in the PWM header
tmp <- tempfile(fileext = ".fa"); writeLines(c(">s", "ACGTACGTACGTAAAAATTTTTCACTGAAAAAA"), tmp)
bgx <- tempfile(); dir.create(bgx)
invisible(system2(file.path(MEME_BIN, "fimo"), c("--bgfile", "--nrdb--", "--thresh", "1", "--oc", bgx,
                  file.path(CH3_DIR, "STRCACTGARARAAA.meme"), tmp), stdout = FALSE, stderr = FALSE))
xml <- readLines(file.path(bgx, "fimo.xml"))
nrdb_a <- as.numeric(sub('.*letter="A">([0-9.]+)<.*', "\\1", grep('letter="A"', xml, value = TRUE)[1]))
hdr <- pwm[grep("^Background letter frequencies", pwm) + 1]
pwm_a <- as.numeric(sub("^A ([0-9.]+) .*", "\\1", hdr))
vals <- list(memeChipCentreLen = centre_len, nMaskedBases = n_rep, nCmeMaskedBases = n_cme,
             nCmeMaskedSites = nrow(bed), canonPwmSites = pwm_sites, canonPwmE = sci(as.numeric(pwm_e)),
             canonVsMemeChipTomtomP = sci(tp),
             permPFloor = sci(1 / (perm_n + 1), 0), permMcSe = sci(sqrt(0.01 * 0.99 / perm_n)),
             fimoBgA = formatC(nrdb_a, format = "f", digits = 3), pwmBgA = formatC(pwm_a, format = "f", digits = 3))
print(unlist(vals))
save_values("ch2_motif2_inputs", vals)
