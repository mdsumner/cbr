#' Burn geometries onto a regular grid
#'
#' Rasterize polygons, lines, and points onto a regular grid using the
#' `controlledburn` Rust crate. Polygon interiors are returned as horizontal
#' runs of fully covered cells, and (in coverage mode) boundary cells are
#' returned with their exact coverage fraction.
#'
#' @param geom Geometry input. Anything that [wk::as_wkb()] understands, e.g. a
#'   `wk_wkb` vector, an `sf`/`sfc` column, or a list of raw WKB blobs.
#' @param extent Grid extent as `c(xmin, xmax, ymin, ymax)` (the ordering used
#'   by `terra::ext()`). If missing, the bounding box of `geom` is used (via
#'   [wk::wk_bbox()]).
#' @param dimension Grid dimensions as `c(ncol, nrow)`. If missing, `ncol` is
#'   256 and `nrow` is computed from it so that cells are square in the units
#'   of `extent` (at least 1 row). A single value is taken as `ncol`, with
#'   `nrow` computed the same way.
#' @param coverage Logical. If `TRUE` (default) use exact coverage fractions
#'   (boundary cells appear in `edges`). If `FALSE` use the faster cell-centre
#'   ("approx", fasterize-style) rule, which emits runs only.
#' @param crs Optional coordinate reference system as a WKT string (or anything
#'   coercible to character). Carried through to [cb_materialize()] output.
#' @param parity Logical. By default a multi-part geometry gets one record per
#'   cell (parts' coverage is summed, and a cell the parts together fill
#'   becomes a run), and adjacent runs on a row are coalesced. If `TRUE`, give
#'   the records of the C++ 'controlledburn' core instead: one run per fully
#'   covered boundary cell and separate records for each part.
#' @param threads Number of threads, `1` (default) for none, `0` for all
#'   available. The result is identical to the serial burn. Threads only pay
#'   for large inputs (the crate uses at most one per 8192 units of work), and
#'   with more than one thread the input is parsed up front, which uses more
#'   memory.
#'
#' @return An object of class `cb_burn`: a named list of five data frames with
#'   the grid definition stored in attributes (`extent`, `dimension`, `crs`,
#'   `coverage`) and the number of input geometries in `n_geom`. Indices are 1-based and row 1 is the top row; geometry `k`
#'   (input position) gets `id = k`.
#'   \describe{
#'     \item{`runs`}{`row`, `col_start`, `col_end`, `id`: fully covered cells.}
#'     \item{`edges`}{`row`, `col`, `fraction`, `id`: boundary cells, fraction in (0, 1).}
#'     \item{`lines`}{`row`, `col`, `length`, `id`: line length within the cell, CRS units.}
#'     \item{`points`}{`row`, `col`, `id`: one row per point inside the grid.}
#'     \item{`notes`}{`geom_index`, `message`: non-fatal per-geometry problems.}
#'   }
#'
#' @examples
#' \dontrun{
#' library(wk)
#' sq <- wkt("POLYGON ((2.5 4.5, 6.5 4.5, 6.5 8.5, 2.5 8.5, 2.5 4.5))")
#' cb_burn(sq, extent = c(0, 10, 0, 10), dimension = c(10, 10))
#'
#' ## defaults: bbox of the geometry, 256 columns, square cells
#' cb_burn(sq)
#' }
#' @export
cb_burn <- function(geom, extent, dimension, coverage = TRUE, crs = NA_character_,
                    parity = FALSE, threads = 1L) {
  geom <- wk::as_wkb(geom)

  if (missing(extent) || is.null(extent)) {
    extent <- as.numeric(wk::wk_bbox(geom))[c(1L, 3L, 2L, 4L)]
  }
  if (length(extent) != 4L) {
    stop("`extent` must be a length-4 vector c(xmin, xmax, ymin, ymax).", call. = FALSE)
  }
  extent <- as.numeric(extent)
  if (anyNA(extent) || extent[2L] <= extent[1L] || extent[4L] <= extent[3L]) {
    stop("`extent` must be finite with xmax > xmin and ymax > ymin.", call. = FALSE)
  }

  if (missing(dimension) || is.null(dimension)) {
    dimension <- 256L
  }
  if (length(dimension) == 1L) {
    dimension <- c(dimension, default_nrow(extent, dimension))
  }
  if (length(dimension) != 2L) {
    stop("`dimension` must be a length-2 vector c(ncol, nrow).", call. = FALSE)
  }
  dimension <- as.integer(dimension)
  if (anyNA(dimension) || any(dimension < 1L)) {
    stop("`dimension` must be positive integers.", call. = FALSE)
  }

  if (length(threads) != 1L || is.na(threads) || threads < 0) {
    stop("`threads` must be a single non-negative number (0 for all available).",
         call. = FALSE)
  }
  threads <- as.integer(threads)

  blobs <- unclass(geom)

  res <- cb_burn_wkb(
    blobs,
    xmin = extent[1L], ymin = extent[3L],
    xmax = extent[2L], ymax = extent[4L],
    ncol = dimension[1L], nrow = dimension[2L],
    coverage = isTRUE(coverage),
    parity = isTRUE(parity),
    threads = threads
  )

  tables <- lapply(res, function(tbl) {
    class(tbl) <- "data.frame"
    attr(tbl, "row.names") <- .set_row_names(length(tbl[[1L]]))
    tbl
  })

  structure(
    tables,
    extent = extent,
    dimension = dimension,
    crs = if (length(crs)) as.character(crs)[1L] else NA_character_,
    coverage = isTRUE(coverage),
    n_geom = length(blobs),
    class = "cb_burn"
  )
}

## Number of rows that gives (approximately) square cells for `ncol` columns
## over `extent`; never fewer than one row.
default_nrow <- function(extent, ncol) {
  xr <- extent[2L] - extent[1L]
  yr <- extent[4L] - extent[3L]
  max(1L, as.integer(round(ncol * yr / xr)))
}

#' @export
print.cb_burn <- function(x, ...) {
  ext <- attr(x, "extent")
  dm <- attr(x, "dimension")
  cat("<cb_burn>\n")
  cat(sprintf("  grid:   %d col x %d row\n", dm[1L], dm[2L]))
  cat(sprintf("  extent: x[%g, %g] y[%g, %g]\n", ext[1L], ext[2L], ext[3L], ext[4L]))
  cat(sprintf("  mode:   %s\n", if (isTRUE(attr(x, "coverage"))) "coverage" else "approx"))
  parent <- attr(x, "parent")
  if (!is.null(parent)) {
    cat(sprintf("  crop:   offset col %d row %d of %d x %d parent\n",
                parent$offset[1L], parent$offset[2L],
                parent$dimension[1L], parent$dimension[2L]))
  }
  cat("  tables:", paste(sprintf("%s=%d", names(x), vapply(x, nrow, integer(1))),
                         collapse = ", "), "\n")
  invisible(x)
}
