## The cbr hex logo, drawn by cbr.
##
## A flame burned onto a lattice of 16 x 16 cells turned 30 degrees, so the
## cells run parallel to the sides of the hex: the case cb_burn_index() is
## for. Every cell's colour is its exact coverage. The inner flame is a second
## burn, blended in by its own coverage. Run from the package root; writes
## man/figures/logo.png.

library(cbr)

charcoal <- "#1c1917"
border <- "#c2410c"
lattice <- "#3f3a36"
core <- "#fef3c7"
ember <- grDevices::colorRampPalette(c("#3b0a0a", "#7f1d1d", "#c2410c", "#f97316"))(101)

## a pointy-top hexagon of radius r, and a test for points inside the unit one
hexagon <- function(r = 1) {
  a <- (90 + 60 * 0:6) * pi / 180
  cbind(r * cos(a), r * sin(a))
}
in_hex <- function(x, y) abs(x) <= sqrt(3) / 2 & abs(y) <= 1 - abs(x) / sqrt(3)

## a smooth closed curve through control points (periodic spline)
smooth_ring <- function(pts, cx = 0, cy = 0, sx = 1, sy = 1, n = 200) {
  pts <- rbind(pts, pts[1, ])
  t <- seq_len(nrow(pts))
  x <- stats::spline(t, pts[, 1], n = n, method = "periodic")$y
  y <- stats::spline(t, pts[, 2], n = n, method = "periodic")$y
  m <- cbind(cx + x * sx, cy + y * sy)
  m[n, ] <- m[1, ]
  m
}
as_wkt <- function(m) {
  wk::wkt(sprintf("POLYGON ((%s))", paste(sprintf("%.8f %.8f", m[, 1], m[, 2]), collapse = ", ")))
}

## the flame: round base, curved tip and a small tongue on the left; the
## core sits low in the base
flame_outer <- rbind(c(0.05, 0.64), c(0.20, 0.36), c(0.38, 0.10), c(0.42, -0.22), c(0.28, -0.48),
                     c(0.00, -0.57), c(-0.27, -0.49), c(-0.42, -0.24), c(-0.40, 0.04), c(-0.27, 0.22),
                     c(-0.27, 0.42), c(-0.10, 0.22))
flame_inner <- rbind(c(0.03, 0.18), c(0.17, -0.06), c(0.19, -0.30), c(0.02, -0.43), c(-0.16, -0.33),
                     c(-0.17, -0.08))
outer <- smooth_ring(flame_outer, 0, 0.10, 1.05, 1.0)
inner <- smooth_ring(flame_inner, 0, 0.10, 1.05, 1.0)

## the lattice: n x n cells of size s, columns turned 30 degrees, centred
n <- 16
s <- 1.45 / n
u <- c(cos(pi / 6), sin(pi / 6))
w <- c(sin(pi / 6), -cos(pi / 6))
centre <- c(0, 0.10)
gt <- cb_geotransform(centre - n / 2 * s * u - n / 2 * s * w, res = s, rotation = 30)
cell <- function(c, r) {
  p <- rbind(c(c, r), c(c + 1, r), c(c + 1, r + 1), c(c, r + 1))
  cbind(gt[1] + p[, 1] * gt[2] + p[, 2] * gt[3], gt[4] + p[, 1] * gt[5] + p[, 2] * gt[6])
}

coverage <- function(ring) {
  b <- cb_burn_index(as_wkt(ring), gt, c(n, n))
  cb_as_matrix(cb_materialize(b, background = 0))
}
cov_outer <- coverage(outer)
cov_inner <- coverage(inner)

blend <- function(a, b, t) {
  ca <- grDevices::col2rgb(a) / 255
  cb <- grDevices::col2rgb(b) / 255
  grDevices::rgb(t(ca * (1 - t) + cb * t))
}

draw_logo <- function(stroke) {
  op <- graphics::par(mar = c(0, 0, 0, 0), ljoin = "mitre", lmitre = 10)
  on.exit(graphics::par(op))
  graphics::plot.new()
  graphics::plot.window(xlim = c(-sqrt(3) / 2, sqrt(3) / 2), ylim = c(-1, 1), asp = 1,
                        xaxs = "i", yaxs = "i")
  r <- 0.955
  graphics::polygon(hexagon(r), col = charcoal, border = NA)

  ## the empty lattice, extended past the burn: whole cells inside the hex
  ## and above the label
  for (row in -n:(2 * n)) for (col in -n:(2 * n)) {
    xy <- cell(col, row)
    if (all(in_hex(xy[, 1] * 1.07, xy[, 2] * 1.07)) && min(xy[, 2]) > -0.5) {
      graphics::polygon(xy, border = lattice, lwd = stroke * 0.07)
    }
  }
  ## the burned cells
  for (row in seq_len(n)) for (col in seq_len(n)) {
    v <- cov_outer[row, col]
    if (v <= 0) next
    fill <- ember[round(min(1, v) * 100) + 1]
    if (cov_inner[row, col] > 0) fill <- blend(fill, core, cov_inner[row, col])
    graphics::polygon(cell(col - 1, row - 1), col = fill, border = charcoal, lwd = stroke * 0.09)
  }

  cex <- 0.56 / graphics::strwidth("cbr", family = "mono", font = 2, cex = 1)
  graphics::text(0, -0.65, "cbr", family = "mono", font = 2, col = "white", cex = cex)
  graphics::polygon(hexagon(r), border = border, lwd = stroke)
}

dir.create("man/figures", showWarnings = FALSE, recursive = TRUE)
## 518 x 600 pixels, the usual hex sticker size; the stroke is in line widths
## (1/96 inch) at 150 dpi, about 22 pixels
ragg::agg_png("man/figures/logo.png", width = 518, height = 600, res = 150, bg = "transparent")
draw_logo(stroke = 14)
invisible(grDevices::dev.off())
