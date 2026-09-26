## ====
## ch3_ch4_reference_set_statistics.R: recomputes reference set analysis statistics; writes CSVs to data/reference_sets/ (refset_composition, refset_tissue_contrast, refset_power_mde, alpha_sharing, refset_sharing_5000, refset_pairwise_5000, refset_motif2_pooled, refset_motif2_pooled_stats, refset_motif2_perhet).
## Runtime ~15 min on 12 cores. Batch job: Rscript 04_chapter_analyses/chapter3/ch3_ch4_reference_set_statistics.R. Permutation and power stats depend on mc.cores stream assignment, agree to about 3rd sig fig at n_perm=5000.
## Trh/Tgo joint peaks=614, peak-proximal genes=693; reference set=2053 genes (811 BDGP-annotated + 1242 screen-derived); BDGP-annotated=398 tracheal-only+313 SG-only+100 both; heterodimer joint sets: Trh 614, Sim 355, Dys 870, Sima 296.
## ====

## Repository root: the environment variable PROJECT_DIR, or the working directory.
PROJECT_DIR <- normalizePath(Sys.getenv("PROJECT_DIR", unset = "."), mustWork = TRUE)
if (!file.exists(file.path(PROJECT_DIR, "README.md")))
  stop("PROJECT_DIR does not point at the repository root: ", PROJECT_DIR)
setwd(PROJECT_DIR)

suppressPackageStartupMessages({
  library(GenomicRanges); library(GenomeInfoDb); library(regioneR)
  library(org.Dm.eg.db);  library(AnnotationDbi); library(parallel)
  library(GenomicFeatures)
  library(TxDb.Dmelanogaster.UCSC.dm6.ensGene)
})

AUD   <- file.path(PROJECT_DIR, "data/reference_sets")
A4    <- file.path(PROJECT_DIR, "data/chapter4")
NPERM <- 5000L
## Core count: NCORE_OVERRIDE if set, otherwise detectCores() or the OS core
## count (up to 12 cores, leaving one free). Override with
##     NCORE_OVERRIDE=8 Rscript ch3_ch4_reference_set_statistics.R

.ok  <- function(x) length(x) == 1L && !is.na(x) && x >= 1L
.try <- function(expr) { v <- suppressWarnings(try(expr, silent = TRUE))
                         if (inherits(v, "try-error")) NA_integer_ else as.integer(v)[1] }
.nc <- .try(as.integer(Sys.getenv("NCORE_OVERRIDE", NA)))
if (!.ok(.nc)) .nc <- .try(detectCores(logical = FALSE))
if (!.ok(.nc)) .nc <- .try(detectCores())
if (!.ok(.nc)) .nc <- .try(system(
  if (Sys.info()[["sysname"]] == "Darwin") "sysctl -n hw.physicalcpu" else "nproc",
  intern = TRUE, ignore.stderr = TRUE))
if (!.ok(.nc)) .nc <- 4L                  # last resort: assume a modest laptop
NCORE <- max(1L, min(12L, .nc - 1L))
RNGkind("L'Ecuyer-CMRG"); set.seed(20260727)
cat(sprintf("n_perm = %d, cores = %d\n", NPERM, NCORE))

good_chr <- c("chr2L","chr2R","chr3L","chr3R","chr4","chrX","chrY")
kp <- function(g) keepSeqlevels(g, intersect(good_chr, seqlevels(g)), pruning.mode = "coarse")
ordH <- c("Trh","Sim","Dys","Sima")

## ---- 1. peak sets, rebuilt from the raw IDR narrowPeak files ----------------
read_np <- function(p) {
  d <- read.table(p, sep = "\t", stringsAsFactors = FALSE)
  GRanges(d$V1, IRanges(d$V2 + 1L, d$V3), signalValue = d$V7)
}
trh <- kp(read_np("data/final_peaks/TRH_IDR0.05.clean.narrowPeak"))
tgo <- kp(read_np("data/final_peaks/TGO_IDR0.05.clean.narrowPeak"))
h     <- findOverlaps(trh, tgo, minoverlap = 1L)
joint <- pintersect(trh[queryHits(h)], tgo[subjectHits(h)])
joint <- kp(joint); joint <- joint[!duplicated(joint)]
## The joint peak set contains 614 peaks.
stopifnot(length(joint) == 614L, length(trh) == 1033L, length(tgo) == 1403L)

alph <- lapply(c(Sim = "SIM", Dys = "DYS", Sima = "SIMA"),
               function(x) kp(read_np(sprintf("data/final_peaks/%s_IDR0.05.clean.narrowPeak", x))))
stopifnot(identical(unname(sapply(alph, length)), c(587L, 1342L, 438L)))

## annotated joint sets carry the Motif2 / enhancer / distance columns
load_ann <- function(nm) {
  d <- read.csv(file.path(A4, paste0(nm, "_Tgo_joint_annotated.csv")), stringsAsFactors = FALSE)
  g <- GRanges(d$seqnames, IRanges(d$start, d$end))
  g$Motif2 <- as.logical(d$Motif2)
  g$distal <- abs(d$distanceToTSS) > 1000
  kp(g)
}
J <- lapply(setNames(ordH, ordH), load_ann)
stopifnot(identical(unname(sapply(J, length)), c(614L, 355L, 870L, 296L)))

## ---- 2. TSS universe and peak-proximal gene set -----------------------------
txdb   <- TxDb.Dmelanogaster.UCSC.dm6.ensGene
tss_gr <- resize(suppressMessages(genes(txdb)), width = 1, fix = "start")
## one-to-many ENSEMBL->FLYBASE maps are dropped rather than arbitrarily resolved
fbl <- suppressMessages(mapIds(org.Dm.eg.db, keys = names(tss_gr), keytype = "ENSEMBL",
                               column = "FLYBASE", multiVals = "list"))
fb  <- suppressMessages(mapIds(org.Dm.eg.db, keys = names(tss_gr), keytype = "ENSEMBL",
                               column = "FLYBASE", multiVals = "first"))
fb[lengths(fbl) > 1] <- NA
tss_gr$flybase <- fb
nearby <- unique(na.omit(tss_gr$flybase[subjectHits(
  findOverlaps(joint, tss_gr, maxgap = 1000L, ignore.strand = TRUE))]))
stopifnot(length(nearby) == 693L)

## ---- 3. BDGP tissue partitions of the reference set -------------------------
## Derived from the BDGP in situ annotation dump (insitu_annot.csv) by matching
## the documented query terms in trach_terms.txt / sg_terms.txt; the resulting
## FBgn lists live alongside this script. The ^FBgn[0-9]+$ filter is load-bearing:
## a blank line in ref_bdgp_any.txt (a `comm` artefact) once inflated the count
## from 811 to 812.
rd <- function(f) { v <- unique(trimws(readLines(file.path(AUD, f)))); v[grepl("^FBgn[0-9]+$", v)] }
sets <- list(full        = rd("ref2053.txt"),      bdgp_any = rd("ref_bdgp_any.txt"),
             trach_any   = rd("ref_trach.txt"),    sg_any   = rd("ref_sg.txt"),
             trach_only  = rd("ref_trach_only.txt"), sg_only = rd("ref_sg_only.txt"),
             both        = rd("ref_both.txt"),     screen_only = rd("ref_screen_only.txt"))
stopifnot(lengths(sets)[c("full","trach_only","sg_only","both","screen_only")] ==
            c(2053L, 398L, 313L, 100L, 1242L),
          length(sets$trach_any) == 498L, length(sets$sg_any) == 413L,
          length(sets$bdgp_any) == 811L)
write.csv(data.frame(
  partition = c("trach_only","sg_only","both","screen_only","bdgp_any","full"),
  n_genes   = c(398L, 313L, 100L, 1242L, 811L, 2053L)),
  file.path(AUD, "refset_composition.csv"), row.names = FALSE)

## ---- 4. gene lengths and the length-matched sampler ------------------------
## Gene length is the confounder of record: tracheal-only genes are ~2x longer
## than SG-only genes. lm_null_one() draws controls from the same length decile
## as each peak-proximal gene; it supplies the null of the power analysis
## (section 8) and gene lengths enter the tissue contrast (section 7).
g    <- suppressMessages(genes(txdb))
fbg  <- suppressMessages(mapIds(org.Dm.eg.db, keys = names(g), keytype = "ENSEMBL",
                                column = "FLYBASE", multiVals = "first"))
glen <- setNames(as.numeric(width(g)), fbg); glen <- glen[!is.na(names(glen))]
univ <- names(glen)
dec    <- setNames(cut(glen, quantile(glen, seq(0, 1, 0.1)), include.lowest = TRUE,
                       labels = FALSE), univ)
nb     <- intersect(nearby, univ)
tab_nb <- table(dec[nb])
by_dec <- split(univ, dec[univ])
lm_null_one <- function(gs) {
  smp <- unlist(lapply(names(tab_nb), function(d) sample(by_dec[[d]], tab_nb[[d]])))
  length(intersect(smp, gs))
}

## ---- 5. alpha-subunit sharing of the 614 Trh/Tgo peaks ----------------------
## Two nulls. The genome-randomised one is the loose bound; the PROMOTER-MATCHED
## one is the honest test, because every factor here binds promoters and peak
## density alone would manufacture overlap.
n_other <- rowSums(sapply(alph, function(x) overlapsAny(joint, x)))
write.csv(data.frame(n_other_alpha = 0:3,
                     n_peaks = as.integer(table(factor(n_other, levels = 0:3)))),
          file.path(AUD, "alpha_sharing.csv"), row.names = FALSE)

dm6_len <- c(chr2L = 23513712, chr2R = 25286936, chr3L = 28110227, chr3R = 32079331,
             chr4 = 1348131, chrX = 23542271, chrY = 3667352)
gen   <- GRanges(names(dm6_len), IRanges(1, dm6_len))
proms <- kp(reduce(trim(resize(tss_gr, width = 2001, fix = "center"))))

share_obs <- mean(n_other > 0)
sh_one <- function(mode) {
  r <- if (mode == "genome")
         randomizeRegions(joint, genome = gen, per.chromosome = TRUE, allow.overlaps = TRUE)
       else resampleRegions(joint, universe = proms, per.chromosome = TRUE)
  mean(rowSums(sapply(alph, function(x) overlapsAny(r, x))) > 0)
}
sh_gw <- unlist(mclapply(seq_len(NPERM), function(i) sh_one("genome"),
                         mc.cores = NCORE, mc.set.seed = TRUE))
sh_pr <- unlist(mclapply(seq_len(NPERM), function(i) sh_one("promoter"),
                         mc.cores = NCORE, mc.set.seed = TRUE))
write.csv(data.frame(
  scenario   = c("observed","genome-randomised","promoter-matched"),
  pct_shared = round(100 * c(share_obs, mean(sh_gw), mean(sh_pr)), 2),
  sd         = c(NA, round(100 * sd(sh_gw), 2), round(100 * sd(sh_pr), 2)),
  fold       = c(NA, round(share_obs / mean(sh_gw), 1), round(share_obs / mean(sh_pr), 1)),
  p_emp      = c(NA, signif((sum(sh_gw >= share_obs) + 1) / (NPERM + 1), 4),
                     signif((sum(sh_pr >= share_obs) + 1) / (NPERM + 1), 4)),
  n_perm     = c(NA, NPERM, NPERM)),
  file.path(AUD, "refset_sharing_5000.csv"), row.names = FALSE)
cat("[1/4] sharing nulls done\n")

## ---- 6. pairwise occupancy overlap ------------------------------------------
## Asymmetric BY CONSTRUCTION: the peak sets differ in size, so these are
## directional containment fractions, not a symmetric similarity measure.
pair_null <- function(a) {
  r <- randomizeRegions(J[[a]], genome = gen, per.chromosome = TRUE, allow.overlaps = TRUE)
  vapply(setdiff(ordH, a), function(b) mean(overlapsAny(r, J[[b]])), numeric(1))
}
pn <- lapply(ordH, function(a)
  rowMeans(simplify2array(mclapply(seq_len(NPERM), function(i) pair_null(a),
                                   mc.cores = NCORE, mc.set.seed = TRUE))))
names(pn) <- ordH
obsM <- outer(ordH, ordH, Vectorize(function(a, b)
  if (a == b) NA_real_ else 100 * mean(overlapsAny(J[[a]], J[[b]]))))
dimnames(obsM) <- list(ordH, ordH)
nullM <- obsM; nullM[] <- NA_real_
for (a in ordH) for (b in setdiff(ordH, a)) nullM[a, b] <- 100 * pn[[a]][[b]]
write.csv(data.frame(row = rep(ordH, 4), col = rep(ordH, each = 4),
                     obs = as.vector(obsM), null = as.vector(nullM),
                     fold = round(as.vector(obsM / nullM), 0)),
          file.path(AUD, "refset_pairwise_5000.csv"), row.names = FALSE)
cat("[2/4] pairwise nulls done\n")

## ---- 7. the direct tissue contrast (FIXED-n; reproduces exactly) ------------
## the INFERENTIAL CLAIM. Comparing one significant and one non-significant test
## against a null does not show that the partitions differ. The defensible
## statement is this DIRECT contrast, adjusted for log gene length.
glen_l <- log10(glen)
cc <- rbind(data.frame(fb = sets$trach_only, grp = "trach_only"),
            data.frame(fb = sets$sg_only,    grp = "sg_only"))
cc      <- cc[cc$fb %in% univ, ]
cc$hit  <- cc$fb %in% nearby
cc$L    <- glen_l[cc$fb]
cc$grp  <- factor(cc$grp, levels = c("trach_only", "sg_only"))
m_adj <- glm(hit ~ grp + L, data = cc, family = binomial)
m_raw <- glm(hit ~ grp,     data = cc, family = binomial)
ci    <- exp(confint.default(m_adj)["grpsg_only", ])
## length across the WHOLE universe is significant, but not within this contrast:
## the outcome is a peak in a fixed +/-1 kb TSS window, which does not scale with
## gene length, so adjustment barely moves the estimate (2.337 -> 2.221).
allg  <- data.frame(hit = univ %in% nearby, L = glen_l[univ])
m_all <- glm(hit ~ L, data = allg, family = binomial)

## minimum detectable effect for the tracheal arm, against its own null.
## Inflating the tracheal rate toward the SG rate would make the arms more
## similar and could never reach 80% power, so power is computed against the null.
gs_t  <- intersect(sets$trach_only, univ)
nullv <- unlist(mclapply(seq_len(NPERM), function(i) lm_null_one(gs_t),
                         mc.cores = NCORE, mc.set.seed = TRUE))
crit  <- quantile(nullv, 0.95); nulmu <- mean(nullv)
mde <- NA_real_
for (f in seq(1.00, 4.00, 0.01)) {
  pw <- mean(rbinom(20000, length(gs_t), min(1, f * nulmu / length(gs_t))) >= crit)
  if (pw >= 0.80) { mde <- f; break }
}
ctr <- data.frame(
  OR_raw = round(exp(coef(m_raw)[["grpsg_only"]]), 3),
  OR_adj = round(exp(coef(m_adj)[["grpsg_only"]]), 3),
  lo = round(ci[[1]], 3), hi = round(ci[[2]], 3),
  p  = signif(summary(m_adj)$coefficients["grpsg_only", 4], 3),
  OR_per_10x_length = round(exp(coef(m_adj)[["L"]]), 3),
  p_length = signif(summary(m_adj)$coefficients["L", 4], 3),
  n_trach = sum(cc$grp == "trach_only"), n_sg = sum(cc$grp == "sg_only"),
  hit_trach = sum(cc$hit[cc$grp == "trach_only"]),
  hit_sg = sum(cc$hit[cc$grp == "sg_only"]),
  med_len_trach = round(median(glen[gs_t], na.rm = TRUE)),
  med_len_sg = round(median(glen[intersect(sets$sg_only, univ)], na.rm = TRUE)),
  MDE_trach_arm = mde,
  OR_length_universe = round(exp(coef(m_all)[["L"]]), 3),
  p_length_universe  = signif(summary(m_all)$coefficients["L", 4], 3))
write.csv(ctr, file.path(AUD, "refset_tissue_contrast.csv"), row.names = FALSE)
cat("[3/4] tissue contrast done\n")

## ---- 8. Motif2 equivalence (FIXED-n) + Monte-Carlo power --------------------
## The POOLED test is primary. Per-heterodimer splits are underpowered
## (6-95%), so four non-significant splits are not four confirmations.
excl <- lapply(setNames(ordH, ordH), function(a)
  !overlapsAny(J[[a]], unlist(GRangesList(unname(J[setdiff(ordH, a)])))))
pool <- do.call(rbind, lapply(ordH, function(a)
  data.frame(het = a, excl = excl[[a]], m2 = J[[a]]$Motif2)))
ft <- fisher.test(table(factor(pool$excl, c(TRUE, FALSE)),
                        factor(pool$m2,   c(TRUE, FALSE))))
write.csv(data.frame(class = c("exclusive","shared"),
                     n = c(sum(pool$excl), sum(!pool$excl)),
                     n_motif2 = c(sum(pool$m2 & pool$excl), sum(pool$m2 & !pool$excl)),
                     pct_motif2 = round(100 * c(mean(pool$m2[pool$excl]),
                                                mean(pool$m2[!pool$excl])), 2)),
          file.path(AUD, "refset_motif2_pooled.csv"), row.names = FALSE)
write.csv(data.frame(OR = round(ft$estimate[[1]], 3),
                     lo = round(ft$conf.int[1], 3), hi = round(ft$conf.int[2], 3),
                     p  = signif(ft$p.value, 3)),
          file.path(AUD, "refset_motif2_pooled_stats.csv"), row.names = FALSE)

## Uses 50000 replicates rather than 1000: at lower replicate counts the estimates
## are unstable between runs, and a number entering the text must be stable.
set.seed(20260727)
pw_stable <- function(ne, ns, ps, f = 1.5, reps = 50000L) {
  x <- rbinom(reps, ne, min(1, ps / f)); y <- rbinom(reps, ns, ps)
  hit <- vapply(seq_len(reps), function(i)
    suppressWarnings(fisher.test(matrix(c(x[i], ne - x[i], y[i], ns - y[i]), 2))$p.value) < 0.05,
    logical(1))
  c(power = mean(hit), se = sqrt(mean(hit) * (1 - mean(hit)) / reps))
}
perhet <- do.call(rbind, lapply(ordH, function(a) {
  e <- excl[[a]]; m <- J[[a]]$Motif2
  r <- pw_stable(sum(e), sum(!e), mean(m[!e]))
  data.frame(heterodimer = a, n_peaks = length(e), n_exclusive = sum(e),
             motif2_pct_exclusive = round(100 * mean(m[e]), 2),
             motif2_pct_shared    = round(100 * mean(m[!e]), 2),
             fisher_p = signif(suppressWarnings(fisher.test(table(
               factor(e, c(TRUE, FALSE)), factor(m, c(TRUE, FALSE))))$p.value), 3),
             power_vs_1.5x = round(100 * r[["power"]], 1),
             power_se = round(100 * r[["se"]], 2))
}))
write.csv(perhet, file.path(AUD, "refset_motif2_perhet.csv"), row.names = FALSE)
rp <- pw_stable(sum(pool$excl), sum(!pool$excl), mean(pool$m2[!pool$excl]))
write.csv(data.frame(stat = c("tracheal_arm_MDE_fold","pooled_power_vs_1.5x_pct"),
                     value = c(mde, round(100 * rp[["power"]], 1))),
          file.path(AUD, "refset_power_mde.csv"), row.names = FALSE)
cat("[4/4] Motif2 equivalence and power done\n")

cat("\nNine CSVs written to data/reference_sets/.\n")
