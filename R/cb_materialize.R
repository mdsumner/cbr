#' Materialize a burn as a dense raster
#'
#' Expand the sparse tables of a [cb_burn()] into a dense cell vector in the
#' layout returned by `gdalraster::read_ds()`: a plain vector, one element per
#' cell, in row-major order starting at the top-left cell, carrying a `"gis"`
#' attribute that describes the grid (`type`, `bbox`, `dim`, `srs`,
#' `datatype`). This is a stable format definition, so no dependency on
#' 'gdalraster' is needed to produce or consume it, and the result can be
#' handed straight to `gdalraster::plot_raster()` or wrapped as a matrix with
#' [cb_as_matrix()].
#'
#' Consider [cb_crop()] first when the parent grid is large and the geometry
#' occupies only part of it.
#'
#' @param x A `cb_burn` object.
#' @param value What each cell should hold:
#'   \describe{
#'     \item{`"coverage"`}{(default) summed coverage. Interior run cells
#'       contribute 1, boundary cells their `fraction`, line cells their
#'       `length` (CRS units), points 1 each. Overlapping geometries add, so
#'       values above 1 are possible; use `pmin(., 1)` if a clamped mask is
#'       wanted. Output `datatype` is `"Float32"`.}
#'     \item{`"id"`}{the 1-based input position of the geometry covering the
#'       cell; where geometries overlap the highest `id` wins (later inputs
#'       are burned over earlier ones). Boundary, line, and point cells count
#'       as covered. Output `datatype` is `"Int32"`.}
#'     \item{`"count"`}{the number of distinct geometries with a record in
#'       the cell. Output `datatype` is `"Int32"`.}
#'     \item{`"intersects"`, `"interiors_meet"`, `"covers"`, `"touches"`,
#'       `"covered_by"`}{the number of polygons satisfying that predicate
#'       against the cell, as in [cb_cells_where()]. Needs a `relations` table
#'       (`cb_burn(relate = TRUE)` or [cb_relate()]). `"intersects"` is the
#'       exact all-touched mask: unlike `"count"` it includes cells a polygon
#'       only touches along a wall or at a corner, and cells holding slivers
#'       the burn rounds away. Output `datatype` is `"Int32"`.}
#'   }
#' @param background Value for cells not touched by any geometry. Default `NA`.
#'
#' @section Affine grids:
#' For a burn from [cb_burn_index()] the cells are rotated or sheared, which
#' the `gis` format cannot describe with a bounding box alone. The vector is
#' still in grid order (row-major from the top-left cell, one element per
#' cell), `bbox` is the bounding box of the grid's outline in world
#' coordinates, and an extra element `geotransform` holds the six GDAL
#' numbers that define where each cell is. Consumers that only read `bbox`
#' and `dim` will place the raster approximately but with square, unrotated
#' cells.
#'
#' @return A numeric (`"coverage"`) or integer (every other `value`) vector of
#'   length `ncol * nrow` with attribute `gis`, a list:
#'   `type = "raster"`, `bbox = c(xmin, ymin, xmax, ymax)`,
#'   `dim = c(ncol, nrow, 1)`, `srs` (the burn's `crs`, or `""`), and
#'   `datatype`.
#'
#' @examples
#' \dontrun{
#' library(wk)
#' sq <- wkt("POLYGON ((2.5 4.5, 6.5 4.5, 6.5 8.5, 2.5 8.5, 2.5 4.5))")
#' b <- cb_burn(sq, extent = c(0, 10, 0, 10), dimension = c(10, 10))
#' v <- cb_materialize(b)
#' cb_as_matrix(v)
#' }
#' @export
cb_materialize <- function(x, value = c("coverage", "id", "count", "intersects",
                                        "interiors_meet", "covers", "touches",
                                        "covered_by"),
                           background = NA) {
  if (!inherits(x, "cb_burn")) {
    stop("`x` must be a <cb_burn> object.", call. = FALSE)
  }
  value <- match.arg(value)

  ext <- attr(x, "extent")
  dm <- attr(x, "dimension")
  ncol <- as.integer(dm[1L])
  nrow <- as.integer(dm[2L])
  ncell <- ncol * nrow

  if (value %in% DE9IM_PREDICATES) {
    ## the number of geometries satisfying the predicate in each cell
    cells <- cb_cells_where(x, value)
    out <- rep(as.integer(background), ncell)
    if (nrow(cells)) {
      s <- rowsum(rep.int(1L, nrow(cells)), (cells$row - 1L) * ncol + cells$col,
                  reorder = FALSE)
      out[as.integer(rownames(s))] <- as.integer(s[, 1L])
    }
    return(with_gis(out, x, "Int32"))
  }

  ## expand every table to (cell index, id, weight), cell = (row - 1) * ncol + col
  runs <- x$runs
  len <- runs$col_end - runs$col_start + 1L
  run_cells <- if (length(len)) {
    (rep.int(runs$row, len) - 1L) * ncol + sequence(len, from = runs$col_start)
  } else integer()

  cell <- c(
    run_cells,
    (x$edges$row - 1L) * ncol + x$edges$col,
    (x$lines$row - 1L) * ncol + x$lines$col,
    (x$points$row - 1L) * ncol + x$points$col
  )
  id <- c(rep.int(runs$id, len), x$edges$id, x$lines$id, x$points$id)
  w <- c(
    rep.int(1, sum(len)),
    x$edges$fraction,
    x$lines$length,
    rep.int(1, nrow(x$points))
  )

  if (value == "coverage") {
    out <- rep(as.numeric(background), ncell)
    if (length(cell)) {
      s <- rowsum(w, cell, reorder = FALSE)
      out[as.integer(rownames(s))] <- s[, 1L]
    }
    datatype <- "Float32"
  } else if (value == "id") {
    out <- rep(as.integer(background), ncell)
    if (length(cell)) {
      o <- order(id)
      out[cell[o]] <- id[o]
    }
    datatype <- "Int32"
  } else {
    out <- rep(as.integer(background), ncell)
    if (length(cell)) {
      u <- !duplicated(cbind(cell, id))
      s <- rowsum(rep.int(1L, sum(u)), cell[u], reorder = FALSE)
      out[as.integer(rownames(s))] <- as.integer(s[, 1L])
    }
    datatype <- "Int32"
  }

  with_gis(out, x, datatype)
}

## attach the read_ds-style `gis` attribute describing the burn's grid
with_gis <- function(out, x, datatype) {
  ext <- attr(x, "extent")
  dm <- attr(x, "dimension")
  crs <- attr(x, "crs")
  gt <- attr(x, "geotransform")
  bbox <- c(ext[1L], ext[3L], ext[2L], ext[4L])
  if (!is.null(gt)) {
    we <- attr(x, "world_extent")
    bbox <- c(we[1L], we[3L], we[2L], we[4L])
  }
  gis <- list(
    type = "raster",
    bbox = bbox,
    dim = c(as.integer(dm[1L]), as.integer(dm[2L]), 1),
    srs = if (is.null(crs) || is.na(crs)) "" else crs,
    datatype = datatype
  )
  if (!is.null(gt)) {
    gis$geotransform <- gt
  }
  attr(out, "gis") <- gis
  out
}

#' Convert a materialized raster to a matrix
#'
#' Reshape the vector returned by [cb_materialize()] (or by
#' `gdalraster::read_ds()`) into a matrix with `nrow` rows and `ncol` columns,
#' oriented as the grid appears on a map (row 1 at the top). The `gis`
#' attribute is dropped.
#'
#' @param x A vector with a `gis` attribute holding `dim = c(ncol, nrow, ...)`.
#' @return A matrix.
#' @export
cb_as_matrix <- function(x) {
  gis <- attr(x, "gis")
  if (is.null(gis) || length(gis$dim) < 2L) {
    stop("`x` must carry a `gis` attribute with a `dim` element.", call. = FALSE)
  }
  ncol <- as.integer(gis$dim[1L])
  nrow <- as.integer(gis$dim[2L])
  v <- as.vector(x)
  attributes(v) <- NULL
  ## read_ds order is row-major from the top-left; R fills column-major
  t(matrix(v[seq_len(ncol * nrow)], nrow = ncol, ncol = nrow))
}
