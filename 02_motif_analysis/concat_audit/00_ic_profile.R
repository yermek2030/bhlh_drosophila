#!/usr/bin/env Rscript
# ============================================================
# Per-position information content of the Motif2 PWM.
# Question: are positions 10-15 (the ARARAAA tail) low-IC filler?
# Decision rule: a position below ~0.3 bits carries essentially no information.
# ============================================================
## Repository root: the environment variable PROJECT_DIR, or the working directory.
PROJECT_DIR <- normalizePath(Sys.getenv("PROJECT_DIR", unset = "."), mustWork = TRUE)
if (!file.exists(file.path(PROJECT_DIR, "README.md")))
  stop("PROJECT_DIR does not point at the repository root: ", PROJECT_DIR)
setwd(PROJECT_DIR)

OUT <- file.path(PROJECT_DIR, "data", "motif2_concat_audit")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

PWM_FILE <- file.path(PROJECT_DIR, "data", "chapter3",
                      "STRCACTGARARAAA.meme")
stopifnot(file.exists(PWM_FILE))

# --- minimal MEME letter-probability-matrix reader (no package dependency) ---
read_meme_lpm <- function(f) {
  x  <- readLines(f)
  i  <- grep("^letter-probability matrix", x)[1]
  stopifnot(!is.na(i))
  w  <- as.integer(sub(".*w= *([0-9]+).*", "\\1", x[i]))
  rows <- x[(i + 1):(i + w)]
  m <- do.call(rbind, lapply(strsplit(trimws(rows), "[[:space:]]+"), as.numeric))
  colnames(m) <- c("A", "C", "G", "T")
  m
}

p <- read_meme_lpm(PWM_FILE)
stopifnot(nrow(p) == 15L, all(abs(rowSums(p) - 1) < 1e-6))

bg_unif <- c(A = .25,  C = .25,  G = .25,  T = .25)
bg_dm6  <- c(A = .306, C = .194, G = .194, T = .306)   # header of the PWM file

ic <- function(p, bg) apply(p, 1, function(r) sum(ifelse(r > 0, r * log2(r / bg), 0)))

cons <- strsplit("STRCACTGARARAAA", "")[[1]]
ic_tab <- data.frame(
  pos        = seq_len(15),
  consensus  = cons,
  max_prob   = round(apply(p, 1, max), 3),
  IC_uniform = round(ic(p, bg_unif), 3),
  IC_dm6bg   = round(ic(p, bg_dm6),  3),
  block      = ifelse(seq_len(15) <= 9, "core_1_9", "tail_10_15")
)
print(ic_tab, row.names = FALSE)

cat("\n--- block totals (bits) ---\n")
tot <- aggregate(cbind(IC_uniform, IC_dm6bg) ~ block, ic_tab, sum)
print(tot, row.names = FALSE)
cat("TOTAL uniform:", round(sum(ic_tab$IC_uniform), 2),
    " TOTAL dm6-bg:", round(sum(ic_tab$IC_dm6bg), 2), "\n")

cat("\n--- positions below 0.3 bits ---\n")
low_u <- ic_tab$pos[ic_tab$IC_uniform < 0.3]
low_d <- ic_tab$pos[ic_tab$IC_dm6bg   < 0.3]
cat("uniform bg :", if (length(low_u)) paste(low_u, collapse = ",") else "none", "\n")
cat("dm6 AT-rich:", if (length(low_d)) paste(low_d, collapse = ",") else "none", "\n")

cat("\nVERDICT: the tail is low-IC filler ONLY IF most of 10-15 fall below 0.3 bits.\n")

write.csv(ic_tab, file.path(OUT, "T0_ic_profile.csv"), row.names = FALSE)

png(file.path(OUT, "T0_ic_profile.png"), width = 1600, height = 900, res = 180)
op <- par(mar = c(4, 4.2, 2.5, 1))
bp <- barplot(ic_tab$IC_dm6bg, names.arg = ic_tab$consensus,
              col = ifelse(ic_tab$pos <= 9, "grey25", "grey70"),
              ylab = "Information content (bits, dm6 background)",
              xlab = "Motif2 position", ylim = c(0, 2.5),
              main = "Motif2 per-position IC (dark = 5' core, light = 3' tail)")
abline(h = 0.3, lty = 2, col = "red")
text(bp, ic_tab$IC_dm6bg + 0.1, ic_tab$pos, cex = 0.6)
par(op); dev.off()

cat("\nWrote:\n  ", file.path(OUT, "T0_ic_profile.csv"),
    "\n  ", file.path(OUT, "T0_ic_profile.png"), "\n")
