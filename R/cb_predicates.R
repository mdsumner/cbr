## DE-9IM predicates per cell, from a burn's runs and its relations table.
##
## A cell the polygon's boundary does not reach is wholly inside the polygon
## (in a run: COVERED) or wholly outside (DISJOINT); the relations table holds
## the matrix of every other cell, so the two together answer for any cell.

DE9IM_COVERED <- "212FF1FF2"
DE9IM_DISJOINT <- "FF2FF1212"
DE9IM_PREDICATES <- c("intersects", "interiors_meet", "covers", "touches", "covered_by")

## A named predicate on DE-9IM strings of polygon (rows) against cell
## (columns), vectorised; NA in, NA out.
de9im_is <- function(m, predicate) {
  set <- function(i) substr(m, i, i) != "F"
  switch(predicate,
    intersects = set(1L) | set(2L) | set(4L) | set(5L),
    interiors_meet = set(1L),
    covers = !set(7L) & !set(8L),
    touches = !set(1L) & (set(2L) | set(4L) | set(5L)),
    covered_by = !set(3L) & !set(6L),
    stop(sprintf("unknown predicate `%s`.", predicate), call. = FALSE)
  )
}

check_pattern <- function(pattern) {
  if (!is.character(pattern) || length(pattern) != 1L || is.na(pattern) ||
      !grepl("^[TtFf*012]{9}$", pattern)) {
    stop("`pattern` must be nine characters from T, F, *, 0, 1 and 2, e.g. \"T*F**F***\".",
         call. = FALSE)
  }
  toupper(pattern)
}

## Does each DE-9IM string match a pattern (T: not F, F: F, *: anything,
## a digit: that dimension)?
de9im_matches <- function(m, pattern) {
  p <- strsplit(check_pattern(pattern), "")[[1L]]
  ok <- !is.na(m)
  for (i in seq_len(9L)) {
    v <- substr(m, i, i)
    ok <- ok & switch(p[i], "*" = TRUE, "T" = v != "F", "F" = v == "F", v == p[i])
  }
  ok
}

require_relations <- function(x, call_name) {
  if (!inherits(x, "cb_burn")) {
    stop("`x` must be a <cb_burn> object.", call. = FALSE)
  }
  if (is.null(x$relations)) {
    stop(sprintf("%s needs a `relations` table: burn with `relate = TRUE`, or add one with `cb_relate()`.",
                 call_name), call. = FALSE)
  }
}

burn_n_geom <- function(x) {
  n <- attr(x, "n_geom")
  if (is.null(n)) {
    n <- max(burn_id_offset(x), x$runs$id, x$edges$id, x$lines$id, x$points$id,
             x$relations$id) - burn_id_offset(x)
  }
  n
}

## A number for each (id, cell), unique within the burn.
cell_key <- function(x, id, row, col) {
  dm <- as.numeric(attr(x, "dimension"))
  ncell <- dm[1L] * dm[2L]
  if (burn_n_geom(x) * ncell >= 2^53) {
    stop("too many geometries times cells to key exactly.", call. = FALSE)
  }
  (id - burn_id_offset(x) - 1) * ncell + (row - 1) * dm[1L] + (col - 1)
}

## Is each (id, row, col) in one of the runs? One findInterval over runs
## sorted by (id, row) group and start column.
in_runs <- function(runs, id, row, col, ncol) {
  if (!nrow(runs) || !length(id)) {
    return(rep(FALSE, length(id)))
  }
  groups <- unique(paste(runs$id, runs$row))
  gi <- match(paste(runs$id, runs$row), groups)
  w <- as.numeric(ncol) + 2
  start <- gi * w + runs$col_start
  o <- order(start)
  start <- start[o]
  end <- (gi * w + runs$col_end)[o]
  qg <- match(paste(id, row), groups)
  q <- qg * w + col
  hit <- !is.na(qg)
  i <- integer(length(q))
  i[hit] <- findInterval(q[hit], start)
  hit[hit] <- i[hit] > 0L & q[hit] <= end[pmax(i[hit], 1L)]
  hit
}

## Every (id, cell) whose matrix is not DISJOINT: the relations, plus the run
## cells the boundary does not touch (COVERED) when `runs` is TRUE.
cell_matrices <- function(x, runs = TRUE) {
  rel <- x$relations
  out <- data.frame(row = rel$row, col = rel$col, id = rel$id, matrix = rel$matrix)
  if (runs && nrow(x$runs)) {
    r <- x$runs
    len <- r$col_end - r$col_start + 1L
    cells <- data.frame(row = rep.int(r$row, len), col = sequence(len, from = r$col_start),
                        id = rep.int(r$id, len))
    ## on the cells the boundary touches, the relation is the answer
    touched <- cell_key(x, cells$id, cells$row, cells$col) %in%
      cell_key(x, rel$id, rel$row, rel$col)
    cells <- cells[!touched, , drop = FALSE]
    cells$matrix <- rep(DE9IM_COVERED, nrow(cells))
    out <- rbind(out, cells)
  }
  out <- out[order(out$id, out$row, out$col), , drop = FALSE]
  attr(out, "row.names") <- .set_row_names(nrow(out))
  out
}

#' DE-9IM matrix of any cell
#'
#' Look up the DE-9IM intersection matrix of a polygon against grid cells,
#' for a burn with a `relations` table (see [cb_relate()]). A cell its
#' boundary touches has its matrix in `relations`; any other cell is
#' `"212FF1FF2"` (covered: wholly inside the polygon) if it is in one of the
#' polygon's runs, and `"FF2FF1212"` (disjoint) if not.
#'
#' @param x A `cb_burn` object with a `relations` table.
#' @param row,col Cell indices (1-based, row 1 at the top), recycled to a
#'   common length.
#' @param id Geometry ids, recycled with `row` and `col`. May be omitted when
#'   the burn has a single geometry.
#'
#' @return A character vector of nine-character DE-9IM strings, the polygon's
#'   interior, boundary and exterior (rows) against the cell's (columns). `NA`
#'   for cells outside the grid, ids not in the burn, and ids of lines and
#'   points (which have no polygon matrix).
#'
#' @seealso [cb_cells_where()] to select cells by predicate, and
#'   `cb_materialize(value = "intersects")` and friends for dense masks.
#' @examples
#' g <- wk::wkt("POLYGON ((1 1, 4 1, 4 4, 1 4, 1 1))")
#' b <- cb_burn(g, extent = c(0, 6, 0, 6), dimension = c(6, 6), relate = TRUE)
#' ## a cell deep inside, one below the bottom edge (sharing a wall), one
#' ## touching a corner, and one far away
#' cb_matrix(b, row = c(4, 6, 6, 1), col = c(3, 3, 1, 6))
#' @export
cb_matrix <- function(x, row, col, id) {
  require_relations(x, "cb_matrix()")
  dm <- attr(x, "dimension")
  off <- burn_id_offset(x)
  n_geom <- burn_n_geom(x)
  if (missing(id)) {
    if (n_geom != 1L) {
      stop("`id` is needed when the burn has more than one geometry.", call. = FALSE)
    }
    id <- off + 1L
  }
  n <- max(length(row), length(col), length(id))
  if (!length(row) || !length(col) || !length(id)) {
    return(character())
  }
  if (n %% length(row) || n %% length(col) || n %% length(id)) {
    stop("`row`, `col` and `id` must recycle to a common length.", call. = FALSE)
  }
  row <- rep_len(as.numeric(row), n)
  col <- rep_len(as.numeric(col), n)
  id <- rep_len(as.numeric(id), n)

  ok <- !is.na(row) & !is.na(col) & !is.na(id) &
    row >= 1 & row <= dm[2L] & col >= 1 & col <= dm[1L] &
    row == trunc(row) & col == trunc(col) & id == trunc(id) &
    id > off & id <= off + n_geom &
    !id %in% c(x$lines$id, x$points$id)

  out <- rep(NA_character_, n)
  out[ok] <- DE9IM_DISJOINT
  out[ok][in_runs(x$runs, id[ok], row[ok], col[ok], dm[1L])] <- DE9IM_COVERED
  rel <- x$relations
  m <- match(cell_key(x, id[ok], row[ok], col[ok]), cell_key(x, rel$id, rel$row, rel$col))
  out[ok][!is.na(m)] <- rel$matrix[m[!is.na(m)]]
  out
}

#' Select the cells that satisfy a spatial predicate
#'
#' The cells of a burn where each polygon and the cell (as a closed box)
#' satisfy a named spatial predicate or a DE-9IM pattern, from the burn's runs
#' and its `relations` table (see [cb_relate()]). Exact: the cells only touched
#' along a wall or at a corner, which have no burn record, are included where
#' the predicate holds for them.
#'
#' Predicates are of the polygon against the cell:
#' \describe{
#'   \item{`"intersects"`}{they share at least a point. The cells are the
#'     exact *all-touched* set.}
#'   \item{`"interiors_meet"`}{the polygon covers part of the cell's area.}
#'   \item{`"covers"`}{no point of the cell is outside the polygon.}
#'   \item{`"touches"`}{they share boundary points but no area: cells across
#'     a wall or a corner from the polygon.}
#'   \item{`"covered_by"`}{the cell covers the polygon (the polygon lies within
#'     the closed cell).}
#' }
#'
#' @param x A `cb_burn` object with a `relations` table.
#' @param predicate One of the predicates above. Ignored if `pattern` is given.
#' @param pattern Optional DE-9IM pattern, nine characters from `T` (not
#'   empty), `F` (empty), `*` (anything) and `0`, `1`, `2` (that dimension),
#'   polygon (rows) against cell (columns). A pattern that the matrix of an
#'   untouched cell (`"FF2FF1212"`) matches is refused: those cells are the
#'   rest of the grid and cannot be listed sparsely.
#' @param id Optional geometry ids to keep.
#'
#' @return A data frame with one row per geometry and cell satisfying the
#'   predicate, sorted by `id`, `row` and `col`: `row`, `col`, `id` and
#'   `matrix` (the DE-9IM string).
#'
#' @examples
#' g <- wk::wkt("POLYGON ((1 1, 3 1, 3 3, 1 3, 1 1))")
#' b <- cb_burn(g, extent = c(0, 4, 0, 4), dimension = c(4, 4), relate = TRUE)
#' cb_cells_where(b, "touches")
#' ## the polygon's boundary runs along the cell's boundary
#' cb_cells_where(b, pattern = "****1****")
#' @export
cb_cells_where <- function(x, predicate = DE9IM_PREDICATES, pattern = NULL, id = NULL) {
  require_relations(x, "cb_cells_where()")
  if (!is.null(pattern)) {
    pattern <- check_pattern(pattern)
    if (de9im_matches(DE9IM_DISJOINT, pattern)) {
      stop("`pattern` matches cells the polygon does not touch (\"FF2FF1212\"), which ",
           "are the rest of the grid; use a pattern that excludes them.", call. = FALSE)
    }
    test <- function(m) de9im_matches(m, pattern)
  } else {
    predicate <- match.arg(predicate)
    test <- function(m) de9im_is(m, predicate)
  }
  cells <- cell_matrices(x, runs = test(DE9IM_COVERED))
  keep <- test(cells$matrix)
  if (!is.null(id)) {
    keep <- keep & cells$id %in% id
  }
  cells <- cells[keep, , drop = FALSE]
  attr(cells, "row.names") <- .set_row_names(nrow(cells))
  cells
}

#' Cells where the burn's records and the exact topology disagree
#'
#' The burn treats coverage within 1e-6 of 0 or 1 as exactly 0 or 1. That is
#' right for area, but it means a cell holding a sliver of a polygon can have
#' no record, and a cell holding a sliver of the polygon's exterior can be in a
#' run. The matrices in the `relations` table come from the geometry, not the
#' records, so comparing the two finds those cells.
#'
#' @param x A `cb_burn` object with a `relations` table (see [cb_relate()]).
#'
#' @return A data frame with one row per geometry and disagreeing cell, sorted
#'   by `id`, `row` and `col`: `row`, `col`, `id`, `matrix` and `kind`, one of
#'   \describe{
#'     \item{`"interior"`}{the polygon shares area with the cell, but the burn
#'       has no record there (a sliver of the polygon, rounded to 0).}
#'     \item{`"exterior"`}{the cell is in a run, but the polygon does not cover
#'       it (a sliver of the exterior, rounded to 1).}
#'   }
#'   No rows when the records agree with the topology everywhere, which is
#'   the usual case.
#'
#' @examples
#' ## a square with a spike of area 5e-8 into the next cell: no record there
#' g <- wk::wkt("POLYGON ((1 1, 3 1, 3 2, 3.0000001 2.5, 3 3, 1 3, 1 1))")
#' b <- cb_burn(g, extent = c(0, 4, 0, 4), dimension = c(4, 4), relate = TRUE)
#' cb_slivers(b)
#' @export
cb_slivers <- function(x) {
  require_relations(x, "cb_slivers()")
  rel <- x$relations
  dm <- attr(x, "dimension")
  key <- cell_key(x, rel$id, rel$row, rel$col)
  run <- in_runs(x$runs, rel$id, rel$row, rel$col, dm[1L])
  edge <- key %in% cell_key(x, x$edges$id, x$edges$row, x$edges$col)
  kind <- rep(NA_character_, nrow(rel))
  kind[!run & !edge & de9im_is(rel$matrix, "interiors_meet")] <- "interior"
  kind[run & !de9im_is(rel$matrix, "covers")] <- "exterior"
  keep <- !is.na(kind)
  out <- data.frame(row = rel$row[keep], col = rel$col[keep], id = rel$id[keep],
                    matrix = rel$matrix[keep], kind = kind[keep])
  out <- out[order(out$id, out$row, out$col), , drop = FALSE]
  attr(out, "row.names") <- .set_row_names(nrow(out))
  out
}
