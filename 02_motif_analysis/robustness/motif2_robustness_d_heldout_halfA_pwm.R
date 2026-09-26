## ============================================================
## motif2_robustness_d_heldout_halfA_pwm.R
## Scans half B with Motif2 matrix (rank 4, STREME on half A) and 50 shuffled copies; repeats scan with STRCACTGARARAAA PWM (103/307 hits). FIMO --thresh 1e-4, >=1 hit/sequence.
## Inputs: data/motif2_robustness/streme_halfA/streme.txt, halfB.fa, halfB_shuf.fa, data/chapter3/STRCACTGARARAAA.meme
## Outputs: data/motif2_robustness/heldout_halfA_pwm/, results/reported_values.tsv section ch4_heldout_pwm
## ============================================================
PROJECT_DIR <- normalizePath(Sys.getenv("PROJECT_DIR", unset = "."), mustWork = TRUE)
if (!file.exists(file.path(PROJECT_DIR, "README.md")))
  stop("PROJECT_DIR does not point at the repository root: ", PROJECT_DIR)

DATA_DIR      <- file.path(PROJECT_DIR, "data")
CH3_DIR <- file.path(DATA_DIR, "chapter3")

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
D   <- file.path(DATA_DIR, "motif2_robustness")
OUT <- file.path(D, "heldout_halfA_pwm"); dir.create(OUT, showWarnings = FALSE)
run <- function(tool, args, stdout = TRUE) {
  r <- system2(file.path(MEME_BIN, tool), args, stdout = stdout, stderr = FALSE)
  if (!is.null(attr(r, "status"))) stop(tool, " failed"); r
}
ver <- run("fimo", "--version")
# 1. half-A PWM (STREME rank 4 on half A), and its STREME statistics
st  <- readLines(file.path(D, "streme_halfA/streme.txt"))
hdr <- grep("^MOTIF 4-VGTRCACTGARARAAAW", st); stopifnot(length(hdr) == 1)
lp  <- st[hdr + 1]
halfA_sites <- as.integer(sub(".*nsites= *([0-9]+).*", "\\1", lp))
halfA_E     <- sub(".*E= *([0-9.eE+-]+).*", "\\1", lp)
halfA_E     <- formatC(as.numeric(halfA_E), format = "fg", digits = 2)
writeLines(run("meme-get-motif", c("-id", "4-VGTRCACTGARARAAAW", file.path(D, "streme_halfA/streme.txt"))),
           file.path(OUT, "halfA_motif2.meme"))
# 2. scans
scan_count <- function(meme, fa) {
  x <- run("fimo", c("--bgfile", "--nrdb--", "--thresh", "1e-4", "--text", meme, fa))
  x <- read.delim(text = paste(x, collapse = "\n"), stringsAsFactors = FALSE)
  x <- x[!grepl("^#", x[[1]]) & nzchar(x$sequence_name), ]
  unique(x$sequence_name)
}
fa_B <- file.path(D, "halfB.fa"); fa_S <- file.path(D, "halfB_shuf.fa")
nB <- sum(grepl("^>", readLines(fa_B))); stopifnot(nB == 307L)
null_mean <- function(hits) {
  rep_id <- sub(".*_shuf_([0-9]+)$", "\\1", hits)
  cnt <- table(factor(rep_id, levels = as.character(1:50)))
  c(mean = mean(cnt), max = max(cnt))
}
canon <- file.path(CH3_DIR, "STRCACTGARARAAA.meme")
obs_c <- length(scan_count(canon, fa_B)); nul_c <- null_mean(scan_count(canon, fa_S))
obs_h <- length(scan_count(file.path(OUT, "halfA_motif2.meme"), fa_B))
nul_h <- null_mean(scan_count(file.path(OUT, "halfA_motif2.meme"), fa_S))
cat(sprintf("%s | canonical PWM: %d/307 vs null %.1f | half-A PWM: %d/307 vs null %.1f (max of 50 shuffles %d)\n",
            ver, obs_c, nul_c[["mean"]], obs_h, nul_h[["mean"]], nul_h[["max"]]))
stopifnot(obs_c == 103L)                                   # canonical-PWM scan of half B
write.csv(data.frame(pwm = c("canonical (all 614 peaks)", "half A only"),
                     heldout_hits = c(obs_c, obs_h), n = nB,
                     null_mean = c(nul_c[["mean"]], nul_h[["mean"]]), null_max = c(nul_c[["max"]], nul_h[["max"]])),
          file.path(OUT, "heldout_summary.csv"), row.names = FALSE)
save_values("ch4_heldout_pwm", list(
  motifTwoHalfARank      = 4L,
  motifTwoHalfASites     = halfA_sites,
  motifTwoHalfAE         = halfA_E,
  motifTwoHalfPwmN       = obs_h,
  motifTwoHalfPwmPct     = formatC(100 * obs_h / nB, format = "f", digits = 1),
  motifTwoHalfPwmNull    = formatC(nul_h[["mean"]], format = "f", digits = 1),
  motifTwoHalfPwmNullMax = nul_h[["max"]],
  motifTwoHalfPwmFold    = formatC(obs_h / nul_h[["mean"]], format = "f", digits = 2),
  motifTwoHeldoutShufReps = 50L
))
