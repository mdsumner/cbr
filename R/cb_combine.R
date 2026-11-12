#' Combine burns of consecutive slices of a collection
#'
#' Concatenate [cb_burn()] results that share a grid and whose geometry ids do
#' not overlap, as when a collection is burned in chunks (one per worker, say)
#' with `id_offset` set to the number of geometries before each chunk. The
#' result is the same, record for record, as burning the whole collection at
#' once: each geometry's records are produced independently, so joining the
#' chunks in id order changes nothing.
#'
#' @param ... `cb_burn` objects, or a single list of them. They may be given in
#'   any order; the tables are joined in order of `id_offset`.
#'
#' @return A `cb_burn` object on the shared grid. Its `id_offset` is the lowest
#'   of the inputs' and `n_geom` spans to the highest id, so if the chunks
#'   leave a gap in the ids those geometries count as empty. `parent` is kept
#'   when all inputs agree on it.
#'
#' @details All inputs must have identical `extent`, `dimension`, `crs` and
#'   `coverage`, and their id ranges `(id_offset, id_offset + n_geom]` must not
#'   overlap; otherwise an error is raised rather than silently merging
#'   different geometries under one id.
#'
#' @examples
#' \dontrun{
#' library(wk)
#' g <- wkt(c("POLYGON ((1 1, 4 1, 4 4, 1 4, 1 1))",
#'            "POLYGON ((5 5, 9 5, 9 9, 5 9, 5 5))"))
#' grid <- list(extent = c(0, 10, 0, 10), dimension = c(10, 10))
#' a <- cb_burn(g[1], extent = grid$extent, dimension = grid$dimension)
#' b <- cb_burn(g[2], extent = grid$extent, dimension = grid$dimension,
#'              id_offset = 1)
#' cb_combine(a, b)
#' }
#' @export
cb_combine <- function(...) {
  xs <- list(...)
  if (length(xs) == 1L && !inherits(xs[[1L]], "cb_burn") && is.list(xs[[1L]])) {
    xs <- xs[[1L]]
  }
  if (!length(xs)) {
    stop("nothing to combine.", call. = FALSE)
  }
  if (!all(vapply(xs, inherits, logical(1), what = "cb_burn"))) {
    stop("all inputs must be <cb_burn> objects.", call. = FALSE)
  }

  first <- xs[[1L]]
  for (a in c("extent", "dimension", "crs", "coverage")) {
    same <- vapply(xs, function(x) identical(attr(x, a), attr(first, a)), logical(1))
    if (!all(same)) {
      stop(sprintf("burns differ in `%s`; they must share one grid.", a), call. = FALSE)
    }
  }

  lo <- vapply(xs, burn_id_offset, integer(1))
  n <- vapply(xs, function(x) {
    attr(x, "n_geom") %||% (max(0L, x$runs$id, x$edges$id, x$lines$id, x$points$id) -
                              burn_id_offset(x))
  }, numeric(1))
  ord <- order(lo)
  xs <- xs[ord]
  lo <- lo[ord]
  hi <- lo + n[ord]
  if (any(lo[-1L] < hi[-length(hi)])) {
    stop("id ranges overlap; burn each chunk with its own `id_offset`.", call. = FALSE)
  }

  tables <- lapply(names(first), function(nm) {
    tbl <- do.call(rbind, lapply(xs, `[[`, nm))
    attr(tbl, "row.names") <- .set_row_names(nrow(tbl))
    tbl
  })
  names(tables) <- names(first)

  parents <- lapply(xs, attr, "parent")
  parent <- if (all(vapply(parents, identical, logical(1), y = parents[[1L]]))) {
    parents[[1L]]
  }

  structure(
    tables,
    extent = attr(first, "extent"),
    dimension = attr(first, "dimension"),
    crs = attr(first, "crs"),
    coverage = attr(first, "coverage"),
    n_geom = as.integer(max(hi) - min(lo)),
    id_offset = as.integer(min(lo)),
    parent = parent,
    class = "cb_burn"
  )
}
