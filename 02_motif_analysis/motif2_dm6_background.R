## ====
## motif2_dm6_background.R
## Rescans Motif2 under FIMO nrdb background (A=T=0.275) and a zero-order dm6 background (fasta-get-markov -m 0, six major arms): (1) 614 joint peaks vs 50
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
suppressPackageStartupMessages({
  library(GenomicRanges); library(Biostrings); library(BSgenome.Dmelanogaster.UCSC.dm6)
})
OUT <- file.path(DATA_DIR, "motif2_dm6_background")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
PWM <- file.path(CH3_DIR, "STRCACTGARARAAA.meme")
std <- c("chr2L", "chr2R", "chr3L", "chr3R", "chr4", "chrX")
run <- function(tool, args) {
  r <- system2(file.path(MEME_BIN, tool), args, stdout = TRUE, stderr = FALSE)
  if (!is.null(attr(r, "status"))) stop(tool, " failed"); r
}

# ---- dm6 zero-order background ------------------------------------------------
BG_DM6 <- file.path(OUT, "dm6_major_arms_order0.bg")
if (!file.exists(BG_DM6)) {
  fa <- file.path(tempdir(), "dm6_major_arms.fa")
  s <- getSeq(BSgenome.Dmelanogaster.UCSC.dm6, std); names(s) <- std
  writeXStringSet(s, fa)
  run("fasta-get-markov", c("-m", "0", "-dna", fa, BG_DM6)); unlink(fa)
}
bgl <- grep("^[ACGT] ", readLines(BG_DM6), value = TRUE)
bgf <- setNames(as.numeric(sub("^[ACGT] +", "", bgl)), substr(bgl, 1, 1))
stopifnot(length(bgf) == 4L, abs(sum(bgf) - 1) < 1e-3)

# ---- scanning -------------------------------------------------------------------
BGS <- c(nrdb = "--nrdb--", dm6 = BG_DM6)
hits <- function(fa, bg) {
  x <- run("fimo", c("--bgfile", bg, "--thresh", "1e-4", "--text", PWM, fa))
  x <- read.delim(text = paste(x, collapse = "\n"), stringsAsFactors = FALSE)
  x <- x[!grepl("^#", x[[1]]) & nzchar(x$sequence_name) & x$p.value < 1e-4, ]
  unique(x$sequence_name)
}
seqnames_fa <- function(fa) sub("^>", "", sub(" .*$", "", grep("^>", readLines(fa), value = TRUE)))
ftab <- function(a, na, b, nb) {
  ft <- fisher.test(matrix(c(a, na - a, b, nb - b), 2))
  c(OR = unname(ft$estimate), lo = ft$conf.int[1], hi = ft$conf.int[2], p = ft$p.value)
}

# (1) joint peaks vs 50 dinucleotide shuffles
R  <- file.path(DATA_DIR, "motif2_robustness")
faJ <- file.path(CH3_DIR, "Trh_Tgo_joint_peaks_filt_raw.fa")
faS <- file.path(R, "joint_dishuf2.fa")
stopifnot(length(seqnames_fa(faJ)) == 614L)
shuf <- lapply(BGS, function(bg) {
  o <- length(hits(faJ, bg)); h <- hits(faS, bg)
  cnt <- as.numeric(table(factor(sub(".*_shuf_([0-9]+)$", "\\1", h), levels = as.character(1:50))))
  c(obs = o, null = mean(cnt), sd = sd(cnt), max = max(cnt), fold = o / mean(cnt), Z = (o - mean(cnt)) / sd(cnt))
})
print(shuf)
stopifnot(shuf$nrdb[["obs"]] == 221)

# (2) joint windows vs unbound enhancer windows; (3) bound vs unbound, pooled
# 201-bp joint-peak windows: centre = (BED start + BED end) %/% 2
jb <- read.delim(file.path(CH3_DIR, "Trh_Tgo_joint_peaks_filt.bed"), header = FALSE)
jc <- (jb$V2 + jb$V3) %/% 2L
jw <- getSeq(BSgenome.Dmelanogaster.UCSC.dm6, GRanges(jb$V1, IRanges(jc - 100L, jc + 100L)))
names(jw) <- paste0("jp", seq_along(jw))
faJW <- file.path(tempdir(), "joint_201bp.fa"); writeXStringSet(jw, faJW)
faEW <- file.path(DATA_DIR, "motif2_attribution", "enh_201bp.fa")
enh_fa <- file.path(CH3_DIR, "enhancers_genomewide.fa")
enh_names <- sub("^>", "", grep("^>", readLines(enh_fa), value = TRUE))
enh_len <- nchar(as.character(readDNAStringSet(enh_fa)))
m <- regmatches(enh_names, regexec("^([^:]+):([0-9]+)-([0-9]+)", enh_names))
cc <- do.call(rbind, lapply(m, function(z) if (length(z) == 4L) z[2:4] else c(NA, NA, NA)))
keep <- enh_len >= 150L
enh_gr <- GRanges(cc[keep, 1], IRanges(as.integer(cc[keep, 2]), as.integer(cc[keep, 3])))
stopifnot(length(enh_gr) == length(seqnames_fa(faEW)))
beds <- c(file.path(CH3_DIR, "Trh_Tgo_joint_peaks_filt.bed"),
          file.path(CH4_DIR, paste0(c("Sim", "Dys", "Sima"), "_Tgo_joint_peaks_filt.bed")))
pooled <- reduce(do.call(c, lapply(beds, function(p) { d <- read.delim(p, header = FALSE); GRanges(d$V1, IRanges(d$V2 + 1L, d$V3)) })))
bound <- countOverlaps(enh_gr, pooled, ignore.strand = TRUE) > 0
enhW <- lapply(BGS, function(bg) {
  hj <- length(hits(faJW, bg)); he <- as.integer(sub("^enh", "", hits(faEW, bg)))
  hf <- logical(length(enh_gr)); hf[he] <- TRUE
  nb <- sum(bound); nu <- sum(!bound); hb <- sum(hf & bound); hu <- sum(hf & !bound)
  list(joint = c(a = hj, n = 614), unbound = c(a = hu, n = nu), bound = c(a = hb, n = nb),
       jv = ftab(hj, 614, hu, nu), bv = ftab(hb, nb, hu, nu))
})
stopifnot(enhW$nrdb$joint[["a"]] == 160, enhW$nrdb$unbound[["a"]] == 2849, enhW$nrdb$unbound[["n"]] == 35410,
          enhW$nrdb$bound[["a"]] == 48, enhW$nrdb$bound[["n"]] == 387)

# (4) symmetric control panel (FASTAs written by ch4_ctrl_motif2_symmetric.R)
C <- file.path(DATA_DIR, "ctrl_tf", "symmetric")
sets <- c("TRH_withTgo", "TWI_withTgo", "VVL_withTgo", "RIB_withTgo", "DA_withTgo",
          "TRH_single", "TWI_single", "VVL_single", "RIB_single", "DA_single",
          "TRH_noTgo", "TWI_noTgo", "VVL_noTgo", "RIB_noTgo", "DA_noTgo", "BG_enhancer_matched", "BG_gcmatched")
sym_tab <- read.csv(file.path(C, "ctrl_symmetric_motif2.csv"), stringsAsFactors = FALSE)
gcs <- lapply(setNames(sets, sets), function(s) {
  x <- readDNAStringSet(file.path(C, paste0(s, ".fa")))
  data.frame(id = names(x), gc = as.numeric(letterFrequency(x, "GC", as.prob = TRUE)), set = s, stringsAsFactors = FALSE) })
G <- do.call(rbind, gcs)
sym <- lapply(names(BGS), function(b) {
  hh <- unlist(lapply(sets, function(s) hits(file.path(C, paste0(s, ".fa")), BGS[[b]])))
  D <- G; D$hit <- D$id %in% hh
  tab <- do.call(rbind, lapply(sets, function(s) { d <- D[D$set == s, ]; data.frame(set = s, n = nrow(d), hits = sum(d$hit)) }))
  ref <- tab[tab$set == "TRH_withTgo", ]
  tab$p_vs_ref <- sapply(seq_len(nrow(tab)), function(i) if (tab$set[i] == "TRH_withTgo") NA else
    ftab(tab$hits[i], tab$n[i], ref$hits, ref$n)[["p"]])
  D$gc10 <- 10 * D$gc; D$set <- relevel(factor(D$set), ref = "TRH_withTgo")
  co <- summary(glm(hit ~ set + gc10, family = binomial, data = D))$coefficients
  tab$adj_or <- sapply(tab$set, function(s) if (s == "TRH_withTgo") 1 else exp(co[paste0("set", s), 1]))
  tab$pct <- 100 * tab$hits / tab$n; tab$bg <- b; tab
})
names(sym) <- names(BGS)
chk <- merge(sym$nrdb, sym_tab[, c("set", "hits")], by = "set", suffixes = c("", "_sym"))
stopifnot(nrow(chk) == length(sets), all(chk$hits == chk$hits_sym),
          abs(sym$nrdb$adj_or[sym$nrdb$set == "BG_gcmatched"] - sym_tab$adj_or[sym_tab$set == "BG_gcmatched"]) < 0.005)
write.csv(do.call(rbind, sym), file.path(OUT, "ctrl_symmetric_nrdb_vs_dm6.csv"), row.names = FALSE)

# ---- summaries and values quoted in the text ----------------------------------------
f1 <- function(x) formatC(x, format = "f", digits = 1); f2 <- function(x) formatC(x, format = "f", digits = 2)
f3 <- function(x) formatC(x, format = "f", digits = 3)
fp <- function(p) if (p < 1e-3) formatC(p, format = "e", digits = 1) else formatC(p, format = "fg", digits = 2, flag = "#")
tg <- function(b) sym[[b]][sym[[b]]$set %in% c("TWI_withTgo", "VVL_withTgo", "RIB_withTgo", "DA_withTgo"), ]
sv <- function(b, s, col) sym[[b]][sym[[b]]$set == s, col]
summ <- data.frame(bg = names(BGS),
  shuf_obs = sapply(shuf, `[[`, "obs"), shuf_null = sapply(shuf, `[[`, "null"), shuf_fold = sapply(shuf, `[[`, "fold"),
  shuf_Z = sapply(shuf, `[[`, "Z"), shuf_max = sapply(shuf, `[[`, "max"),
  joint_pct = sapply(enhW, function(e) 100 * e$joint[["a"]] / e$joint[["n"]]),
  unbound_pct = sapply(enhW, function(e) 100 * e$unbound[["a"]] / e$unbound[["n"]]),
  bound_pct = sapply(enhW, function(e) 100 * e$bound[["a"]] / e$bound[["n"]]),
  jv_OR = sapply(enhW, function(e) e$jv[["OR"]]), jv_lo = sapply(enhW, function(e) e$jv[["lo"]]), jv_hi = sapply(enhW, function(e) e$jv[["hi"]]),
  jv_p = sapply(enhW, function(e) e$jv[["p"]]),
  bv_OR = sapply(enhW, function(e) e$bv[["OR"]]), bv_lo = sapply(enhW, function(e) e$bv[["lo"]]), bv_hi = sapply(enhW, function(e) e$bv[["hi"]]),
  bv_p = sapply(enhW, function(e) e$bv[["p"]]),
  ref_pct = sapply(names(BGS), function(b) sv(b, "TRH_withTgo", "pct")),
  tgo_pct_min = sapply(names(BGS), function(b) min(tg(b)$pct)), tgo_pct_max = sapply(names(BGS), function(b) max(tg(b)$pct)),
  tgo_p_min = sapply(names(BGS), function(b) min(tg(b)$p_vs_ref)),
  gcbg_pct = sapply(names(BGS), function(b) sv(b, "BG_gcmatched", "pct")),
  gcbg_adj = sapply(names(BGS), function(b) sv(b, "BG_gcmatched", "adj_or")),
  bgenh_pct = sapply(names(BGS), function(b) sv(b, "BG_enhancer_matched", "pct")))
print(summ, digits = 3, row.names = FALSE)
write.csv(summ, file.path(OUT, "motif2_nrdb_vs_dm6_summary.csv"), row.names = FALSE)

d6 <- summ[summ$bg == "dm6", ]
save_values("ch4_dm6_background", list(
  dmBgA = f3(bgf[["A"]]), dmBgC = f3(bgf[["C"]]),
  dmBgShufObs = as.integer(d6$shuf_obs), dmBgShufNull = f1(d6$shuf_null), dmBgShufFold = f2(d6$shuf_fold),
  dmBgShufZ = f1(d6$shuf_Z), nrdbShufFold = f2(summ$shuf_fold[summ$bg == "nrdb"]),
  nrdbShufNull = f1(summ$shuf_null[summ$bg == "nrdb"]),
  dmBgJointPct = f1(d6$joint_pct), dmBgUnboundPct = f1(d6$unbound_pct), dmBgBoundPct = f1(d6$bound_pct),
  dmBgJointOR = f2(d6$jv_OR), dmBgJointORlo = f2(d6$jv_lo), dmBgJointORhi = f2(d6$jv_hi),
  dmBgPooledOR = f2(d6$bv_OR), dmBgPooledORlo = f2(d6$bv_lo), dmBgPooledORhi = f2(d6$bv_hi), dmBgPooledP = fp(d6$bv_p),
  dmBgRefPct = f1(d6$ref_pct), dmBgTgoPctMin = f1(d6$tgo_pct_min), dmBgTgoPctMax = f1(d6$tgo_pct_max),
  dmBgTgoPMin = f2(d6$tgo_p_min), dmBgGcbgPct = f1(d6$gcbg_pct), dmBgGcbgAdjOR = f2(d6$gcbg_adj)))
cat(sprintf("dm6 background: A=%.3f C=%.3f G=%.3f T=%.3f\n", bgf[["A"]], bgf[["C"]], bgf[["G"]], bgf[["T"]]))
