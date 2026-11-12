#' Build a GDAL-style geotransform
#'
#' A geotransform is the six numbers `c(x0, dx, rx, y0, ry, dy)` that map a
#' pixel position (`col`, `row`, counted from the top-left corner of the grid,
#' 0-based) to world coordinates:
#'
#' ```
#' x = x0 + col * dx + row * rx
#' y = y0 + col * ry + row * dy
#' ```
#'
#' This helper builds one from an origin, cell sizes, a rotation and a shear,
#' so that rotated and sheared grids can be written down without doing the
#' trigonometry by hand. A north-up grid is `rotation = 0, shear = 0`, which
#' gives `c(x0, res, 0, y0, 0, -res)`.
#'
#' @param origin Position of the top-left corner of the grid, `c(x, y)`.
#' @param res Cell size: one value, or `c(width, height)`, in world units
#'   (positive).
#' @param rotation Angle in degrees, counter-clockwise, of the direction in
#'   which columns increase, measured from the +x axis.
#' @param shear Extra angle in degrees by which the direction in which rows
#'   increase is turned from being perpendicular to the columns (0 leaves the
#'   cells rectangular).
#'
#' @return A numeric vector of length 6 in GDAL order.
#'
#' @examples
#' cb_geotransform(c(100, 200), res = 10)
#' cb_geotransform(c(0, 0), res = c(2, 1), rotation = 30)
#' @export
cb_geotransform <- function(origin = c(0, 0), res = 1, rotation = 0, shear = 0) {
  if (length(origin) != 2L || !all(is.finite(origin))) {
    stop("`origin` must be two finite numbers, c(x, y).", call. = FALSE)
  }
  if (!length(res) %in% 1:2 || !all(is.finite(res)) || any(res <= 0)) {
    stop("`res` must be one or two positive finite numbers.", call. = FALSE)
  }
  if (length(rotation) != 1L || !is.finite(rotation) ||
      length(shear) != 1L || !is.finite(shear)) {
    stop("`rotation` and `shear` must be single finite numbers (degrees).", call. = FALSE)
  }
  res <- rep_len(as.numeric(res), 2L)
  th <- rotation * pi / 180
  ph <- (rotation + shear) * pi / 180
  lin <- c(res[1L] * cos(th), res[2L] * sin(ph), res[1L] * sin(th), -res[2L] * cos(ph))
  ## cos(90 degrees) is 6e-17, not 0: keep right angles exact
  lin[abs(lin) < 1e-14 * max(abs(lin))] <- 0
  c(as.numeric(origin[1L]), lin[1L], lin[2L], as.numeric(origin[2L]), lin[3L], lin[4L])
}

#' Burn geometries onto an affine (rotated or sheared) grid
#'
#' Rasterize onto a grid defined by a geotransform rather than by an
#' axis-aligned extent. The geometries are moved into the grid's index space
#' (fractional column and row) by the inverse of the geotransform and burned
#' on a lattice of unit cells; the records are then already in the grid's own
#' `(row, col)`. Because the map is affine, straight edges stay straight and
#' the coverage fractions are exact: each is the fraction of the (possibly
#' rotated, sheared or non-square) cell's area covered by the geometry, up to
#' the single precision of the `edges` table.
#'
#' The result is an ordinary [cb_burn()] object, so [cb_crop()],
#' [cb_combine()], [cb_materialize()] and [summary()] work on it. Its `extent`
#' attribute is the index lattice `c(0, ncol, 0, nrow)` (so [plot()] draws in
#' index space), and the world definition is kept in attributes:
#' \describe{
#'   \item{`geotransform`}{the six numbers, as given.}
#'   \item{`world_extent`}{`c(xmin, xmax, ymin, ymax)` of the grid's outline in
#'     world coordinates (its bounding box, for a rotated grid).}
#'   \item{`crs`}{as given.}
#' }
#' Use [cb_cells()] to get the cells as polygons in world coordinates, for
#' plotting or export.
#'
#' @section Line lengths:
#' The `length` column of `lines` is in world units when the cells are
#' similar to squares (equal sides at right angles, whatever the rotation).
#' For other affine grids the length of a piece depends on its direction, and
#' the crate reports lengths in cell units, which is what is returned, with a
#' warning.
#'
#' @param geom Geometry input, anything [wk::as_wkb()] understands, in the
#'   world coordinates the geotransform maps to.
#' @param geotransform Numeric vector of length 6, `c(x0, dx, rx, y0, ry, dy)`,
#'   as in GDAL; see [cb_geotransform()] to build one. The grid's top-left
#'   corner is `(x0, y0)`.
#' @param dimension Grid dimensions as `c(ncol, nrow)`.
#' @inheritParams cb_burn
#'
#' @return A `cb_burn` object; see [cb_burn()].
#'
#' @examples
#' library(wk)
#' gt <- cb_geotransform(c(0, 6), res = 1, rotation = 20)
#' sq <- wkt("POLYGON ((3 3, 7 3, 7 7, 3 7, 3 3))")
#' b <- cb_burn_index(sq, gt, dimension = c(12, 8))
#' b
#' attr(b, "world_extent")
#' @export
cb_burn_index <- function(geom, geotransform, dimension, coverage = TRUE,
                          crs = NA_character_, parity = FALSE, threads = 1L,
                          id_offset = 0L, relate = FALSE) {
  gt <- check_geotransform(geotransform)
  if (missing(dimension) || length(dimension) != 2L ||
      anyNA(dimension) || any(dimension < 1)) {
    stop("`dimension` must be two positive whole numbers, c(ncol, nrow).", call. = FALSE)
  }
  dimension <- as.integer(dimension)
  ncol <- dimension[1L]
  nrow <- dimension[2L]

  ax <- gt[2L]; bx <- gt[3L]; ay <- gt[5L]; by <- gt[6L]
  idx <- wk::wk_transform(wk::as_wkb(geom), index_transform(gt, nrow))

  out <- cb_burn(idx, extent = c(0, ncol, 0, nrow), dimension = dimension,
                 coverage = coverage, crs = crs, parity = parity,
                 threads = threads, id_offset = id_offset, relate = relate)

  if (nrow(out$lines)) {
    la <- sqrt(ax^2 + ay^2)
    lb <- sqrt(bx^2 + by^2)
    if (abs(la - lb) <= 1e-9 * la && abs(ax * bx + ay * by) <= 1e-9 * la * lb) {
      out$lines$length <- out$lines$length * la
    } else {
      warning("cells are not similar to squares: `lines$length` is in cell ",
              "units, not world units.", call. = FALSE)
    }
  }

  attr(out, "geotransform") <- gt
  attr(out, "world_extent") <- affine_world_extent(gt, dimension)
  out
}

## world -> index space, with y up and the top row at y = nrow, so the
## regular-grid burn on c(0, ncol, 0, nrow) numbers rows from the top:
##   x_idx = (by (X - x0) - bx (Y - y0)) / det
##   y_idx = nrow + (ay (X - x0) - ax (Y - y0)) / det
## the origin is subtracted first so large world coordinates keep precision
index_transform <- function(gt, nrow) {
  ax <- gt[2L]; bx <- gt[3L]; ay <- gt[5L]; by <- gt[6L]
  det <- ax * by - bx * ay
  lin <- matrix(c(by / det, -bx / det, 0,
                  ay / det, -ax / det, 0,
                  0, 0, 1), 3L, 3L, byrow = TRUE)
  wk::wk_affine_compose(
    wk::wk_affine_translate(-gt[1L], -gt[4L]),
    wk::wk_trans_affine(lin),
    wk::wk_affine_translate(0, nrow)
  )
}

check_geotransform <- function(gt) {
  if (!is.numeric(gt) || length(gt) != 6L || !all(is.finite(gt))) {
    stop("`geotransform` must be six finite numbers, c(x0, dx, rx, y0, ry, dy).",
         call. = FALSE)
  }
  gt <- as.numeric(gt)
  det <- gt[2L] * gt[6L] - gt[3L] * gt[5L]
  scale <- sqrt(gt[2L]^2 + gt[5L]^2) * sqrt(gt[3L]^2 + gt[6L]^2)
  if (!is.finite(det) || abs(det) <= 1e-12 * scale || scale == 0) {
    stop("`geotransform` is singular: the column and row directions are parallel.",
         call. = FALSE)
  }
  gt
}

## bounding box of a grid's outline in world coordinates, c(xmin, xmax, ymin, ymax)
affine_world_extent <- function(gt, dimension) {
  col <- c(0, dimension[1L], 0, dimension[1L])
  row <- c(0, 0, dimension[2L], dimension[2L])
  x <- gt[1L] + col * gt[2L] + row * gt[3L]
  y <- gt[4L] + col * gt[5L] + row * gt[6L]
  c(range(x), range(y))
}

## the geotransform of any burn: stored for affine burns, implied by the
## extent for regular ones
burn_geotransform <- function(x) {
  gt <- attr(x, "geotransform")
  if (!is.null(gt)) {
    return(gt)
  }
  ext <- attr(x, "extent")
  dm <- attr(x, "dimension")
  c(ext[1L], (ext[2L] - ext[1L]) / dm[1L], 0, ext[4L], 0, -(ext[4L] - ext[3L]) / dm[2L])
}

#' Cells of a burn as polygons in world coordinates
#'
#' One quadrilateral per record cell, in the grid's world coordinates: the
#' cells of every run expanded, then the edge, line and point cells. For a
#' burn from [cb_burn_index()] the quadrilaterals are rotated or sheared as
#' the grid is; for a regular [cb_burn()] they are rectangles. Intended for
#' plotting (`ggplot2::geom_polygon()`) and for export.
#'
#' @param x A `cb_burn` object.
#' @param max_cells Refuse to build more cells than this (runs can cover very
#'   many cells), to avoid exhausting memory by accident.
#'
#' @return A data frame with four rows per cell and columns `cell` (cell
#'   number, the polygon group), `x`, `y` (the corner, in order around the
#'   cell), `row`, `col`, `id`, `kind` (`run`, `edge`, `line` or `point`) and
#'   `value` (1 for runs and points, the coverage fraction for edges, the
#'   length for lines).
#'
#' @examples
#' b <- cb_burn(wk::wkt("POLYGON ((1.5 1.5, 4.5 1.5, 4.5 3.5, 1.5 1.5))"),
#'              extent = c(0, 6, 0, 5), dimension = c(6, 5))
#' head(cb_cells(b), 8)
#' @export
cb_cells <- function(x, max_cells = 1e5) {
  if (!inherits(x, "cb_burn")) {
    stop("`x` must be a <cb_burn> object.", call. = FALSE)
  }
  gt <- burn_geotransform(x)
  runs <- x$runs
  len <- runs$col_end - runs$col_start + 1L
  nrun <- sum(len)

  tab <- data.frame(
    row = c(rep.int(runs$row, len), x$edges$row, x$lines$row, x$points$row),
    col = c(sequence(len, from = runs$col_start), x$edges$col, x$lines$col, x$points$col),
    id = c(rep.int(runs$id, len), x$edges$id, x$lines$id, x$points$id),
    kind = factor(c(rep("run", nrun), rep("edge", nrow(x$edges)),
                    rep("line", nrow(x$lines)), rep("point", nrow(x$points))),
                  levels = c("run", "edge", "line", "point")),
    value = c(rep(1, nrun), x$edges$fraction, x$lines$length, rep(1, nrow(x$points)))
  )
  n <- nrow(tab)
  if (n > max_cells) {
    stop(sprintf("%d cells exceeds `max_cells` (%d).", n, as.integer(max_cells)),
         call. = FALSE)
  }

  ## corners in pixel coordinates: (col-1, row-1), (col, row-1), (col, row), (col-1, row)
  pc <- rep(tab$col, each = 4L) + rep(c(-1, 0, 0, -1), n)
  pr <- rep(tab$row, each = 4L) + rep(c(-1, -1, 0, 0), n)
  cbind(
    data.frame(cell = rep(seq_len(n), each = 4L),
               x = gt[1L] + pc * gt[2L] + pr * gt[3L],
               y = gt[4L] + pc * gt[5L] + pr * gt[6L]),
    tab[rep(seq_len(n), each = 4L), , drop = FALSE],
    row.names = NULL
  )
}
