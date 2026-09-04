#' Crop a burn to a sub-window of its grid
#'
#' Snap an arbitrary extent to the cell boundaries of a [cb_burn()] grid and
#' subset the burn tables to that window. This is a pure index operation: no
#' geometry is re-rasterized, runs that straddle the window edge are clipped,
#' and row/column indices are re-based so that the result is a valid `cb_burn`
#' on the smaller grid. Use it before [cb_materialize()] to avoid allocating a
#' dense array for the whole parent grid.
#'
#' @param x A `cb_burn` object.
#' @param extent Crop extent as `c(xmin, xmax, ymin, ymax)`, in the CRS of the
#'   grid. It is snapped to the parent grid and intersected with it. If missing,
#'   the burn is cropped to its content: the cells touched by any geometry
#'   (see `content` in [summary.cb_burn()]).
#' @param snap How to align `extent` to cell boundaries. `"out"` (default)
#'   expands to the enclosing cells, `"in"` shrinks to fully contained cells,
#'   `"near"` moves each edge to the nearest cell boundary.
#'
#' @return A `cb_burn` object on the cropped grid. The `extent` and `dimension`
#'   attributes describe the new grid; `crs` and `coverage` are carried through.
#'   The parent grid is recorded in attribute `parent` as a list with
#'   `extent`, `dimension`, and the 0-based cell `offset` (`c(col, row)`) of
#'   the window within it. An error is raised if the snapped window is empty.
#'
#' @examples
#' \dontrun{
#' library(wk)
#' sq <- wkt("POLYGON ((2.5 4.5, 6.5 4.5, 6.5 8.5, 2.5 8.5, 2.5 4.5))")
#' b <- cb_burn(sq, extent = c(0, 10, 0, 10), dimension = c(10, 10))
#' cb_crop(b, c(2.2, 5.1, 4.9, 9))   # snapped to x[2, 6] y[4, 9]
#' }
#' @export
cb_crop <- function(x, extent, snap = c("out", "in", "near")) {
  if (!inherits(x, "cb_burn")) {
    stop("`x` must be a <cb_burn> object.", call. = FALSE)
  }
  if (missing(extent) || is.null(extent)) {
    content <- content_extent(x)
    if (is.null(content)) {
      stop("burn is empty, no content to crop to.", call. = FALSE)
    }
    extent <- content$extent
  }
  if (length(extent) != 4L) {
    stop("`extent` must be a length-4 vector c(xmin, xmax, ymin, ymax).", call. = FALSE)
  }
  snap <- match.arg(snap)
  extent <- as.numeric(extent)

  ext0 <- attr(x, "extent")
  dm0 <- attr(x, "dimension")
  ncol0 <- dm0[1L]
  nrow0 <- dm0[2L]
  resx <- (ext0[2L] - ext0[1L]) / ncol0
  resy <- (ext0[4L] - ext0[3L]) / nrow0

  ## fractional cell-boundary positions of the requested window:
  ## columns count from xmin, rows count from ymax (row 1 is the top row)
  cpos <- (c(extent[1L], extent[2L]) - ext0[1L]) / resx
  rpos <- (ext0[4L] - c(extent[4L], extent[3L])) / resy

  c01 <- snap_bounds(cpos, snap)
  r01 <- snap_bounds(rpos, snap)

  ## intersect with the parent grid (0-based boundary indices in [0, n])
  c0 <- max(c01[1L], 0L); c1 <- min(c01[2L], ncol0)
  r0 <- max(r01[1L], 0L); r1 <- min(r01[2L], nrow0)

  if (c1 <= c0 || r1 <= r0) {
    stop("crop window does not intersect the grid (or is empty after snapping).",
         call. = FALSE)
  }

  ncol1 <- c1 - c0
  nrow1 <- r1 - r0
  ext1 <- c(ext0[1L] + c0 * resx, ext0[1L] + c1 * resx,
            ext0[4L] - r1 * resy, ext0[4L] - r0 * resy)

  in_rows <- function(row) row > r0 & row <= r1
  in_cols <- function(col) col > c0 & col <= c1

  ## runs: clip to the column window, drop those that fall out
  runs <- x$runs
  keep <- in_rows(runs$row) & runs$col_end > c0 & runs$col_start <= c1
  runs <- runs[keep, , drop = FALSE]
  runs$col_start <- pmax(runs$col_start, c0 + 1L) - c0
  runs$col_end <- pmin(runs$col_end, c1) - c0
  runs$row <- runs$row - r0

  cellwise <- function(tbl) {
    keep <- in_rows(tbl$row) & in_cols(tbl$col)
    tbl <- tbl[keep, , drop = FALSE]
    tbl$row <- tbl$row - r0
    tbl$col <- tbl$col - c0
    tbl
  }

  out <- list(
    runs = reset_rownames(runs),
    edges = reset_rownames(cellwise(x$edges)),
    lines = reset_rownames(cellwise(x$lines)),
    points = reset_rownames(cellwise(x$points)),
    notes = x$notes
  )

  structure(
    out,
    extent = ext1,
    dimension = c(ncol1, nrow1),
    crs = attr(x, "crs"),
    coverage = attr(x, "coverage"),
    n_geom = attr(x, "n_geom"),
    parent = list(extent = ext0, dimension = dm0, offset = c(c0, r0)),
    class = "cb_burn"
  )
}

## Snap a pair of fractional boundary positions (lo, hi) to integer cell
## boundaries. A relative tolerance keeps positions that are already on a
## boundary (up to floating point noise) from being pushed out a whole cell.
snap_bounds <- function(pos, snap, tol = 1e-8) {
  lo <- pos[1L]
  hi <- pos[2L]
  b <- switch(snap,
    out = c(floor(lo + tol), ceiling(hi - tol)),
    "in" = c(ceiling(lo - tol), floor(hi + tol)),
    near = round(c(lo, hi))
  )
  as.integer(b)
}

reset_rownames <- function(tbl) {
  attr(tbl, "row.names") <- .set_row_names(length(tbl[[1L]]))
  tbl
}
