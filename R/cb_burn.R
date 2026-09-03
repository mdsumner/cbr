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
#'   by `terra::ext()`).
#' @param dimension Grid dimensions as `c(ncol, nrow)`.
#' @param coverage Logical. If `TRUE` (default) use exact coverage fractions
#'   (boundary cells appear in `edges`). If `FALSE` use the faster cell-centre
#'   ("approx", fasterize-style) rule, which emits runs only.
#' @param crs Optional coordinate reference system as a WKT string (or anything
#'   coercible to character). Carried through to [cb_materialize()] output.
#'
#' @return An object of class `cb_burn`: a named list of five data frames with
#'   the grid definition stored in attributes (`extent`, `dimension`, `crs`,
#'   `coverage`). Indices are 1-based and row 1 is the top row; geometry `k`
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
#' }
#' @export
cb_burn <- function(geom, extent, dimension, coverage = TRUE, crs = NA_character_) {
  if (length(extent) != 4L) {
    stop("`extent` must be a length-4 vector c(xmin, xmax, ymin, ymax).", call. = FALSE)
  }
  if (length(dimension) != 2L) {
    stop("`dimension` must be a length-2 vector c(ncol, nrow).", call. = FALSE)
  }

  blobs <- unclass(wk::as_wkb(geom))
  extent <- as.numeric(extent)
  dimension <- as.integer(dimension)

  res <- cb_burn_wkb(
    blobs,
    xmin = extent[1L], ymin = extent[3L],
    xmax = extent[2L], ymax = extent[4L],
    ncol = dimension[1L], nrow = dimension[2L],
    coverage = isTRUE(coverage)
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
    class = "cb_burn"
  )
}

#' @export
print.cb_burn <- function(x, ...) {
  ext <- attr(x, "extent")
  dm <- attr(x, "dimension")
  cat("<cb_burn>\n")
  cat(sprintf("  grid:   %d col x %d row\n", dm[1L], dm[2L]))
  cat(sprintf("  extent: x[%g, %g] y[%g, %g]\n", ext[1L], ext[2L], ext[3L], ext[4L]))
  cat(sprintf("  mode:   %s\n", if (isTRUE(attr(x, "coverage"))) "coverage" else "approx"))
  cat("  tables:", paste(sprintf("%s=%d", names(x), vapply(x, nrow, integer(1))),
                         collapse = ", "), "\n")
  invisible(x)
}
