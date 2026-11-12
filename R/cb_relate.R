#' DE-9IM relations of polygons to the cells their boundaries touch
#'
#' A burn says how much of each cell a polygon covers. The DE-9IM matrix says
#' how the polygon and a cell (a closed box) meet: interiors, boundaries and
#' exteriors, each against each. For almost every cell the burn already
#' answers that: a cell the polygon's boundary does not reach is either wholly
#' inside it (a run) or wholly outside (no record). `cb_relate()` computes the
#' matrix for the remaining cells, those whose closed box the boundary
#' touches, and stores them in a sixth table, `relations`. That includes cells
#' the polygon only touches along a wall or at a corner, which have no record
#' in a burn at all.
#'
#' The matrix is read from the geometry, not from the burn's records. The burn
#' treats a coverage within 1e-6 of 0 or 1 as exactly 0 or 1, so a cell holding
#' a sliver of a polygon has no record and a cell holding a sliver of its
#' exterior is a run; the relations see both.
#'
#' @section The relations table:
#' One row for each polygon and in-grid cell whose closed box the polygon's
#' boundary meets, sorted by `id`, `row` and `col`:
#' \describe{
#'   \item{`row`, `col`, `id`}{as in the other tables: 1-based, row 1 at the
#'     top, `id` the position of the geometry plus `id_offset`.}
#'   \item{`matrix`}{the DE-9IM matrix of the polygon against the cell as nine
#'     characters, row-major over the polygon's interior, boundary and
#'     exterior (rows) against the cell's (columns), each `F` or the
#'     dimension `0`, `1` or `2` of the intersection, for example
#'     `"212F11FF2"`.}
#' }
#' Every other cell of a polygon is `"212FF1FF2"` if it is in one of that
#' polygon's runs and `"FF2FF1212"` otherwise, so the runs and the relations
#' together give the matrix of every cell. Lines, points and empty geometries
#' have no relations.
#'
#' @section Assumptions:
#' Polygons are assumed valid (rings that do not cross, holes inside their
#' shell, parts of a multipolygon meeting only at points); the matrices of
#' invalid polygons are not defined. Contacts are decided with exact
#' predicates on the coordinates as given, with no snapping, so a polygon
#' offset from a cell wall by `1e-14` does not touch it. For a burn from
#' [cb_burn_index()] the geometry is moved into index space exactly as it was
#' burned, and the matrix is exact for the rotated or sheared cell; on grids
#' burned in index space by other front ends it is the index-space matrix.
#'
#' @param x A `cb_burn` object that has not been cropped.
#' @param geom The geometries that were burned into `x`, in the same
#'   coordinates and order: anything [wk::as_wkb()] understands. There must be
#'   `n_geom` of them. The burn does not keep its geometry, so it is passed
#'   again; passing different geometry gives meaningless matrices.
#'
#' @return `x`, with a `relations` table. Use `cb_burn(relate = TRUE)` to get
#'   the same from the geometry in one step.
#'
#' @examples
#' library(wk)
#' g <- wkt("POLYGON ((1 1, 3 1, 3 3, 1 3, 1 1))")
#' b <- cb_relate(cb_burn(g, extent = c(0, 4, 0, 4), dimension = c(4, 4)), g)
#' b$relations
#' table(b$relations$matrix)
#' @export
cb_relate <- function(x, geom) {
  if (!inherits(x, "cb_burn")) {
    stop("`x` must be a <cb_burn> object.", call. = FALSE)
  }
  if (!is.null(attr(x, "parent"))) {
    stop("cannot relate a cropped burn: relate before cropping (the relations ",
         "are cropped along with the other tables).", call. = FALSE)
  }
  geom <- wk::as_wkb(geom)
  n <- attr(x, "n_geom")
  if (!is.null(n) && length(geom) != n) {
    stop(sprintf("`geom` has %d geometries but the burn was of %d.", length(geom), n),
         call. = FALSE)
  }

  dm <- attr(x, "dimension")
  ext <- attr(x, "extent")
  gt <- attr(x, "geotransform")
  if (!is.null(gt)) {
    geom <- wk::wk_transform(geom, index_transform(gt, dm[2L]))
  }

  runs <- x$runs
  runs <- data.frame(row = as.integer(runs$row), col_start = as.integer(runs$col_start),
                     col_end = as.integer(runs$col_end), id = as.integer(runs$id))
  rel <- cb_relate_wkb(
    unclass(geom), runs,
    xmin = ext[1L], ymin = ext[3L], xmax = ext[2L], ymax = ext[4L],
    ncol = dm[1L], nrow = dm[2L],
    coverage = isTRUE(attr(x, "coverage")),
    id_offset = burn_id_offset(x)
  )
  x$relations <- as_burn_table(rel)
  x
}
