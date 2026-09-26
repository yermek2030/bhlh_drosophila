## ============================================================
##  fig1_2_bhlh_pas_architecture.R
##  Figure 1.2: bHLH-PAS heterodimer schematic with labels, leader lines and pentamer strip.
##  Structure drawing (870 x 1640 px) placed unchanged, rest drawn in base graphics.
##  Input: data/chapter1/fig1_bhlh_pas_architecture_raw.png
##  Output: results/figures/fig1_2_bhlh_pas_architecture.png. Needs: png, ragg
## ============================================================

## Repository root: the environment variable PROJECT_DIR, or the working directory.
PROJECT_DIR <- normalizePath(Sys.getenv("PROJECT_DIR", unset = "."), mustWork = TRUE)
if (!file.exists(file.path(PROJECT_DIR, "README.md")))
  stop("PROJECT_DIR does not point at the repository root: ", PROJECT_DIR)
FIGURES_DIR <- file.path(PROJECT_DIR, "results", "figures")
dir.create(FIGURES_DIR, recursive = TRUE, showWarnings = FALSE)

IMG_PATH <- file.path(PROJECT_DIR, "data", "chapter1", "fig1_bhlh_pas_architecture_raw.png")
OUT_PATH <- file.path(FIGURES_DIR, "fig1_2_bhlh_pas_architecture.png")
stopifnot(file.exists(IMG_PATH))

library(png)

img <- readPNG(IMG_PATH)
H <- dim(img)[1]; W <- dim(img)[2]          # pixel height / width

## ---- Colours and sizes ----------------------------------------------
## All text and lines are drawn in the ink colour used throughout the thesis
## figures; the two subunits are coloured in the drawing itself.
DARK    <- "#333333"
HL_FILL <- "lightblue"   # shading behind the "CME core" label

## Text sizes (cex, relative to the device pointsize). With pointsize 6 at
## 300 dpi, cex 1.0 is about 25 px tall.
CEX      <- 1.55     # label text
CEX_HEAD <- 1.55     # subunit headers
LWD_LEAD <- 2        # leader-line width
PTSIZE   <- 6        # device pointsize

## Label boxes: no border, white fill behind the text.
BOX_BORDER <- NA
BOX_FILL   <- "white"

## Shaded box behind the "CME core" label (fixed size, in px).
HL_EBOX <- TRUE
HL_PAD  <- 12            # padding around the text when HL_W / HL_H are NA
HL_W    <- 460           # box width  (NA = fit the text plus HL_PAD)
HL_H    <- 265           # box height
HL_DX   <- 0             # shift of the box centre (px; +right)
HL_DY   <- -135          # shift of the box centre (px; +down)

## Canvas margins (px) around the drawing.
mL <- 470; mR <- 500; mT <- 200; mB <- 400

## ============================================================
## Labels with leader lines: x_anc/y_anc point on drawing (image px, origin top-left, y down); x_txt/y_txt label box centre (x<0 or >W = margin); adj justification (0 left, .5 centre, 1 right)
## ============================================================
lab <- data.frame(
  text = c(
    "PAS-B\n(cofactor/\nligand binding)",
    "PAS-A\n(primary dimerisation\ninterface)",
    "Helix 2\n(packing helix)",
    "basic region (Helix 1)",
    "PAS-B\n(cofactor/\nligand binding)",
    "PAS-A\n(primary dimerisation\ninterface)",
    "Helix 2\n(packing helix)",
    "basic region (Helix 1)",
    "CME core",
    "specificity-defining\nhalf-site (AC; \u03b1)",
    "permissive half-site\n(GTG; Tgo)"
  ),
  ## Rows 10 and 11 point at the brackets of the pentamer strip below the
  ## duplex (SEQ_* block), so their anchors follow SEQ_Y / SEQ_DX / SEQ_BR_DY.
  x_anc = c(205, 205, 235, 255,  665, 665, 630, 615,  433, 340, 557),
  y_anc = c(200, 640,1040,1330,  200, 640,1040,1330, 1500,1754,1754),
  x_txt = c(-40, -40, -40, -40, W+40,W+40,W+40,W+40, W/2, 60, W-60),
  y_txt = c(200, 640,1040,1330,  200, 640,1040,1330, 1620,1670,1670),
  adj   = c(  1,   1,   1,   1,    0,   0,   0,   0,   .5,  1,    0),
  col   = DARK,
  lead  = TRUE,
  stringsAsFactors = FALSE
)

## Text without leader lines: subunit headers and strand polarity.
free <- data.frame(
  text = c("\u03b1-subunit\n(Trh / Sim / Dys / Sima)",
           "Tgo\n(obligate \u03b2-subunit;\nARNT orthologue)",
           "5\u2032","3\u2032","3\u2032","5\u2032"),
  x    = c( 20, W-20,  -15,  -15, W+15, W+15),
  y    = c(-55,  -55, 1380, 1540, 1375, 1530),
  adj  = c(  0,    1,   .5,   .5,   .5,   .5),
  cex  = c(CEX_HEAD, CEX_HEAD, CEX, CEX, CEX, CEX),
  col  = DARK,
  bold = TRUE,
  stringsAsFactors = FALSE
)

## ---- Boxed label ----------------------------------------------------
draw_boxed <- function(x, y, txt, adj = .5, col = DARK, cex = CEX,
                       bg = BOX_FILL, border = BOX_BORDER, pad = 6, font = 1) {
  lines_v <- strsplit(txt, "\n", fixed = TRUE)[[1]]
  lh <- abs(strheight("Ag", cex = cex)) * 1.35   # abs(): the y axis is inverted
  ws <- max(strwidth(lines_v, cex = cex, font = font))
  n  <- length(lines_v)
  bw <- ws + 2 * pad; bh <- lh * n + 2 * pad
  xl <- x - adj * bw; xr <- xl + bw
  yt <- y - bh / 2;   yb <- y + bh / 2
  brd <- if (identical(border, "match")) col else border
  if (!(is.na(brd) && is.na(bg)))
    rect(xl, yb, xr, yt, col = bg, border = brd, lwd = 2)
  ty <- y - (n - 1) / 2 * lh
  tx <- xl + pad + adj * ws
  for (i in seq_len(n)) {
    text(tx, ty, lines_v[i], adj = c(adj, .5), col = col, cex = cex, font = font)
    ty <- ty + lh
  }
  invisible(c(xl = xl, xr = xr, yt = yt, yb = yb))
}

## ============================================================
##  Render
## ============================================================
ragg::agg_png(OUT_PATH, width = W + mL + mR, height = H + mT + mB,
              units = "px", res = 300, pointsize = PTSIZE, background = "white")
op <- par(mar = c(0, 0, 0, 0), xaxs = "i", yaxs = "i")

plot.new()
plot.window(xlim = c(-mL, W + mR), ylim = c(H + mB, -mT))   # y increases downward

## Shading behind the "CME core" label, drawn first so that it sits behind
## the drawing (the PNG is transparent where there is no structure).
EBOX_ROW <- which(lab$text == "CME core")
stopifnot(length(EBOX_ROW) == 1L)
if (isTRUE(HL_EBOX)) {
  .ln  <- strsplit(lab$text[EBOX_ROW], "\n", fixed = TRUE)[[1]]
  .lh  <- abs(strheight("Ag", cex = CEX)) * 1.35
  .ws  <- max(strwidth(.ln, cex = CEX)); .pad <- 6
  .bw  <- .ws + 2 * .pad; .bh <- .lh * length(.ln) + 2 * .pad
  .xl  <- lab$x_txt[EBOX_ROW] - lab$adj[EBOX_ROW] * .bw; .xr <- .xl + .bw
  .yt  <- lab$y_txt[EBOX_ROW] - .bh / 2; .yb <- lab$y_txt[EBOX_ROW] + .bh / 2
  .cx  <- (.xl + .xr) / 2 + HL_DX; .cy <- (.yt + .yb) / 2 + HL_DY
  .hw  <- if (is.na(HL_W)) (.bw / 2 + HL_PAD) else HL_W / 2
  .hh  <- if (is.na(HL_H)) (.bh / 2 + HL_PAD) else HL_H / 2
  rect(.cx - .hw, .cy + .hh, .cx + .hw, .cy - .hh, col = HL_FILL, border = NA)
}

rasterImage(img, 0, H, W, 0)

## Note on the gap between the monomers: the spacing is schematic; in the
## AlphaFold 3 model (Figure 1.3) the subunits interdigitate.
SHOW_GAP_NOTE <- TRUE
GAP_X    <- W / 2
GAP_Y    <- 850
GAP_HALF <- 95
GAP_TXT  <- "spacing schematic;\nsubunits interdigitate\n(see AlphaFold3\nmodel, Fig. 1.3)"
if (SHOW_GAP_NOTE) {
  arrows(GAP_X - GAP_HALF, GAP_Y, GAP_X + GAP_HALF, GAP_Y,
         code = 3, length = .06, angle = 20, lwd = 2.2, col = DARK)
  gl <- strsplit(GAP_TXT, "\n", fixed = TRUE)[[1]]
  lh <- abs(strheight("Ag", cex = CEX * .88)) * 1.3
  ty <- GAP_Y + lh * 1.4
  for (ln in gl) { text(GAP_X, ty, ln, adj = c(.5, .5), cex = CEX * .88,
                        font = 3, col = DARK); ty <- ty + lh }
}

## ---- Pentamer strip -------------------------------------------------
## The CME pentamer ACGTG, the register scanned throughout the thesis, with
## the half-sites bracketed: AC read by the alpha-subunit, GTG by Tgo
## (Swanson et al. 1995). Rows 10 and 11 of `lab` label the brackets.
SEQ_STRIP <- TRUE
SEQ_TXT   <- c("A", "C", "G", "T", "G")
SEQ_ALPHA <- 1:2          # positions read by the alpha-subunit
SEQ_TGO   <- 3:5          # positions read by Tgo
SEQ_X     <- W / 2        # centre of the strip (px)
SEQ_Y     <- 1800         # baseline of the letter row (px)
SEQ_DX    <- 62           # letter pitch (px)
SEQ_HALF  <- 17           # bracket overhang past the outer letter centre (px)
SEQ_BR_DY <- 46           # bracket offset above the letter row (px)
SEQ_TICK  <- 14           # bracket tick length (px)
SEQ_CEX   <- CEX * 1.05   # letter size
SEQ_LCEX  <- CEX * 0.82   # 5' / 3' marks

if (isTRUE(SEQ_STRIP)) {
  seq_x <- SEQ_X + (seq_along(SEQ_TXT) - (length(SEQ_TXT) + 1) / 2) * SEQ_DX

  draw_bracket <- function(pos, y, dy, col = DARK) {
    x0 <- seq_x[min(pos)] - SEQ_HALF; x1 <- seq_x[max(pos)] + SEQ_HALF
    segments(x0, y, x1, y, col = col, lwd = LWD_LEAD)
    segments(c(x0, x1), y, c(x0, x1), y + dy, col = col, lwd = LWD_LEAD)
    invisible((x0 + x1) / 2)
  }

  draw_bracket(SEQ_ALPHA, SEQ_Y - SEQ_BR_DY, SEQ_TICK)
  draw_bracket(SEQ_TGO,   SEQ_Y - SEQ_BR_DY, SEQ_TICK)

  text(seq_x, SEQ_Y, SEQ_TXT, adj = c(.5, .5), font = 2, cex = SEQ_CEX, col = DARK)
  text(seq_x[1] - SEQ_DX, SEQ_Y, "5\u2032", adj = c(.5, .5), cex = SEQ_LCEX, col = DARK)
  text(seq_x[length(seq_x)] + SEQ_DX, SEQ_Y, "3\u2032", adj = c(.5, .5),
       cex = SEQ_LCEX, col = DARK)
}

## Labels with leader lines
for (i in seq_len(nrow(lab))) {
  bx <- draw_boxed(lab$x_txt[i], lab$y_txt[i], lab$text[i],
                   adj = lab$adj[i], col = lab$col[i], cex = CEX,
                   bg = if (isTRUE(HL_EBOX) && i == EBOX_ROW) NA else BOX_FILL)
  if (isTRUE(lab$lead[i])) {
    ex <- if (lab$x_anc[i] >= lab$x_txt[i]) bx["xr"] else bx["xl"]
    ey <- max(bx["yt"], min(bx["yb"], lab$y_anc[i]))
    segments(ex, ey, lab$x_anc[i], lab$y_anc[i],
             col = lab$col[i], lwd = LWD_LEAD)
    points(lab$x_anc[i], lab$y_anc[i], pch = 20, cex = 1.2, col = lab$col[i])
  }
}

## Headers and strand-polarity marks
for (i in seq_len(nrow(free))) {
  lines_v <- strsplit(free$text[i], "\n", fixed = TRUE)[[1]]
  lh <- abs(strheight("Ag", cex = free$cex[i])) * 1.35
  ty <- free$y[i] - (length(lines_v) - 1) / 2 * lh
  for (ln in lines_v) {
    text(free$x[i], ty, ln, adj = c(free$adj[i], .5),
         col = free$col[i], cex = free$cex[i],
         font = if (isTRUE(free$bold[i])) 2 else 1)
    ty <- ty + lh
  }
}

par(op); dev.off()
cat("written:", OUT_PATH, "\n")
