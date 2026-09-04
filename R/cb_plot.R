#' Plot a burn
#'
#' Draw a [cb_burn()] in map coordinates (`asp = 1`) with the grid extent as
#' the frame. Two renderings are available:
#'
#' \describe{
#'   \item{`"exact"`}{every record is drawn: runs as rectangles, edge cells as
#'     rectangles shaded by coverage fraction, line cells shaded by length,
#'     points as cells. Cost is one rectangle per record, so this is the
#'     picture of the sparse representation itself.}
#'   \item{`"overview"`}{the coarse mean-coverage raster from [cb_overview()]
#'     drawn with [graphics::image()]. Cost does not depend on grid size.}
#' }
#'
#' `"auto"` (default) picks `"exact"` when the total record count is at most
#' `threshold`, else `"overview"`. Either way the optional marginal strips are
#' exact: covered cells per grid row (right) and per grid column (top), from
#' [summary.cb_burn()].
#'
#' @param x A `cb_burn` object.
#' @param what `"auto"`, `"exact"`, or `"overview"`.
#' @param marginals Draw the row and column marginal strips. The layout is
#'   sized from the extent's aspect ratio so the strips span exactly the map;
#'   leftover device space is left blank.
#' @param col Colour ramp (a vector of colours, low to high coverage) for
#'   `"overview"` and for edge shading in `"exact"`.
#' @param by_id In `"exact"` mode, colour runs and edges by geometry `id`
#'   instead of by coverage.
#' @param tiles Passed to [cb_overview()] for the `"overview"` rendering.
#' @param threshold Record count above which `"auto"` switches to `"overview"`.
#' @param content Draw the content extent (bounding box of touched cells) as a
#'   dashed rectangle.
#' @param main Optional title, drawn in the outer margin above the strips.
#' @param ... Passed to [graphics::plot()] for the main panel.
#'
#' @return `x`, invisibly. The summary computed for the marginals is attached
#'   as attribute `summary`.
#' @export
plot.cb_burn <- function(x, what = c("auto", "exact", "overview"),
                         marginals = TRUE,
                         col = grDevices::hcl.colors(64, "YlOrRd", rev = TRUE),
                         by_id = FALSE, tiles = 256L, threshold = 5e4,
                         content = FALSE, main = NULL, ...) {
  what <- match.arg(what)
  ext <- attr(x, "extent")
  dm <- attr(x, "dimension")
  nrec <- sum(vapply(x[c("runs", "edges", "lines", "points")], nrow, integer(1)))
  if (what == "auto") what <- if (nrec <= threshold) "exact" else "overview"

  s <- summary(x, tiles = tiles)

  if (marginals) {
    op <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(op), add = TRUE)
    marginal_layout(ext, strip = 1 / 6, main = main)
    graphics::par(mar = c(3.5, 3.5, 0.5, 0.5), mgp = c(2, 0.6, 0))
  } else if (!is.null(main)) {
    op <- graphics::par(mar = c(3.5, 3.5, 2.5, 0.5))
    on.exit(graphics::par(op), add = TRUE)
  }

  graphics::plot(NA, xlim = ext[1:2], ylim = ext[3:4], asp = 1,
                 xlab = "", ylab = "", xaxs = "i", yaxs = "i", bty = "n", ...)
  ## asp = 1 widens one axis to fill the panel: the strips must share these
  ## limits, not the extent, to stay aligned with the map
  usr <- graphics::par("usr")

  if (what == "overview") {
    draw_overview(s$overview, col = col)
  } else {
    draw_exact(x, col = col, by_id = by_id)
  }
  graphics::rect(ext[1L], ext[3L], ext[2L], ext[4L], border = "grey30")
  if (content && !is.null(s$content)) {
    ce <- s$content$extent
    graphics::rect(ce[1L], ce[3L], ce[2L], ce[4L], border = "grey30", lty = 2)
  }

  if (marginals) {
    resx <- (ext[2L] - ext[1L]) / dm[1L]
    resy <- (ext[4L] - ext[3L]) / dm[2L]
    xc <- ext[1L] + (seq_len(dm[1L]) - 0.5) * resx
    yc <- ext[4L] - (seq_len(dm[2L]) - 0.5) * resy

    ## top strip: column marginal against x
    graphics::par(mar = c(0.5, 3.5, 0.5, 0.5))
    graphics::plot(NA, xlim = usr[1:2], ylim = c(0, max(1, s$col_marginal)),
                   xaxs = "i", yaxs = "i", axes = FALSE, xlab = "", ylab = "")
    graphics::polygon(c(ext[1L], xc, ext[2L]), c(0, s$col_marginal, 0),
                      col = "grey70", border = NA)
    graphics::axis(2, at = c(0, max(1, s$col_marginal)), las = 1, cex.axis = 0.7,
                   labels = c("", fmt_num(max(s$col_marginal))))

    ## right strip: row marginal against y
    graphics::par(mar = c(3.5, 0.5, 0.5, 0.5))
    graphics::plot(NA, ylim = usr[3:4], xlim = c(0, max(1, s$row_marginal)),
                   xaxs = "i", yaxs = "i", axes = FALSE, xlab = "", ylab = "")
    graphics::polygon(c(0, s$row_marginal, 0), c(ext[4L], yc, ext[3L]),
                      col = "grey70", border = NA)
    graphics::axis(1, at = c(0, max(1, s$row_marginal)), cex.axis = 0.7,
                   labels = c("", fmt_num(max(s$row_marginal))))
  }

  if (!is.null(main)) graphics::title(main = main, outer = marginals)

  attr(x, "summary") <- s
  invisible(x)
}

## Layout sized to the extent's aspect ratio: main panel bottom-left, column
## strip above it, row strip to its right, spacer cells for whatever the device
## has left over so the strips are never longer than the map.
marginal_layout <- function(ext, strip = 1 / 6, main = NULL) {
  ## title lives in the outer margin; layout() divides what is left
  oma_top <- if (is.null(main)) 0 else 2
  graphics::par(oma = c(0, 0, oma_top, 0))
  din <- grDevices::dev.size("in")
  ## shave a little so rounding never asks for more than the device has
  W <- din[1L] * 0.99
  H <- (din[2L] - oma_top * graphics::par("csi")) * 0.99
  aspect <- (ext[4L] - ext[3L]) / (ext[2L] - ext[1L])
  sw <- strip * min(W, H)                  # strip thickness
  ## main panel: as large as fits with the extent's aspect (margins make this
  ## approximate; asp = 1 in the main plot takes care of the remainder)
  wm <- min(W - sw, (H - sw) / aspect)
  hm <- wm * aspect
  widths <- c(wm, sw, max(0.01, W - wm - sw))
  heights <- c(sw, hm, max(0.01, H - hm - sw))
  ## panels: 1 main, 2 column strip, 3 row strip; 0 = empty
  mat <- rbind(c(2L, 0L, 0L),
               c(1L, 3L, 0L),
               c(0L, 0L, 0L))
  graphics::layout(mat, widths = graphics::lcm(widths * 2.54),
                   heights = graphics::lcm(heights * 2.54))
  invisible(NULL)
}

## image() of a read_ds-shaped vector, in map coordinates
draw_overview <- function(ov, col) {
  gis <- attr(ov, "gis")
  ntx <- gis$dim[1L]
  nty <- gis$dim[2L]
  bb <- gis$bbox
  m <- cb_as_matrix(ov)
  ## image() wants z[i, j] at x[i], y[j] with y ascending: transpose and flip rows
  z <- t(m)[, rev(seq_len(nty)), drop = FALSE]
  xs <- seq(bb[1L], bb[3L], length.out = ntx + 1L)
  ys <- seq(bb[2L], bb[4L], length.out = nty + 1L)
  graphics::image(xs, ys, z, col = col, zlim = c(0, max(1, z, na.rm = TRUE)),
                  add = TRUE, useRaster = TRUE)
}

## one rect per record
draw_exact <- function(x, col, by_id) {
  ext <- attr(x, "extent")
  dm <- attr(x, "dimension")
  resx <- (ext[2L] - ext[1L]) / dm[1L]
  resy <- (ext[4L] - ext[3L]) / dm[2L]
  x0 <- function(col) ext[1L] + (col - 1L) * resx
  y1 <- function(row) ext[4L] - (row - 1L) * resy
  ramp <- function(v) col[pmax(1L, pmin(length(col), ceiling(v * length(col))))]
  idcol <- function(id) grDevices::hcl.colors(max(1L, attr(x, "n_geom") %||% max(id)),
                                              "Dark 3")[id]

  r <- x$runs
  if (nrow(r)) {
    graphics::rect(x0(r$col_start), y1(r$row) - resy, x0(r$col_end + 1L), y1(r$row),
                   col = if (by_id) idcol(r$id) else col[length(col)], border = NA)
  }
  e <- x$edges
  if (nrow(e)) {
    graphics::rect(x0(e$col), y1(e$row) - resy, x0(e$col + 1L), y1(e$row),
                   col = if (by_id) with_alpha(idcol(e$id), e$fraction)
                         else ramp(e$fraction), border = NA)
  }
  l <- x$lines
  if (nrow(l)) {
    diag_len <- sqrt(resx^2 + resy^2)
    graphics::rect(x0(l$col), y1(l$row) - resy, x0(l$col + 1L), y1(l$row),
                   col = if (by_id) idcol(l$id) else ramp(l$length / diag_len),
                   border = NA)
  }
  p <- x$points
  if (nrow(p)) {
    graphics::rect(x0(p$col), y1(p$row) - resy, x0(p$col + 1L), y1(p$row),
                   col = if (by_id) idcol(p$id) else col[length(col)], border = NA)
  }
}

## vectorized alpha (adjustcolor() takes a scalar alpha only)
with_alpha <- function(cols, alpha) {
  m <- grDevices::col2rgb(cols)
  grDevices::rgb(m[1L, ], m[2L, ], m[3L, ], alpha = pmin(1, pmax(0, alpha)) * 255,
                 maxColorValue = 255)
}

`%||%` <- function(a, b) if (is.null(a)) b else a
