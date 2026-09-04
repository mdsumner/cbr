#' Summarize a burn without materializing it
#'
#' Compute a compact description of a [cb_burn()] object in time proportional
#' to the number of records in its tables, never to the number of grid cells.
#' The pieces are:
#'
#' \describe{
#'   \item{`overview`}{a coarse raster of mean coverage per tile, in the same
#'     `read_ds` layout as [cb_materialize()] (see [cb_overview()]). When the
#'     grid is small enough that one tile is one cell, this is the exact
#'     coverage.}
#'   \item{`row_marginal`, `col_marginal`}{exact covered-cell totals per grid
#'     row and per grid column, at full resolution.}
#'   \item{`content`}{the extent and 1-based row/col range of touched cells, or
#'     `NULL` for an empty burn. [cb_crop()] uses this as its default window.}
#'   \item{`area`}{per-geometry totals: `interior` (run cells), `boundary`
#'     (summed edge fractions), `lines` (summed length), `points` (count).}
#'   \item{`fractions`}{quantiles of the edge coverage fractions.}
#'   \item{`records`}{row counts of the five tables.}
#' }
#'
#' Coverage here means: run cells 1 each, edge cells their fraction, line and
#' point cells 1 each (presence).
#'
#' @param object A `cb_burn` object.
#' @param tiles Target number of tiles along the longer grid axis for the
#'   overview. The tile size in cells is `ceiling(max(dimension) / tiles)`.
#' @param ... Ignored.
#'
#' @return A list of class `summary.cb_burn`.
#' @export
summary.cb_burn <- function(object, tiles = 256L, ...) {
  x <- object
  ext <- attr(x, "extent")
  dm <- as.integer(attr(x, "dimension"))
  ncol <- dm[1L]
  nrow <- dm[2L]

  runs <- x$runs
  len <- runs$col_end - runs$col_start + 1L

  ## ---- row marginal: rowsum by row, O(records) ---------------------------
  row_marginal <- numeric(nrow)
  row_marginal <- accumulate_into(row_marginal, runs$row, len)
  row_marginal <- accumulate_into(row_marginal, x$edges$row, x$edges$fraction)
  row_marginal <- accumulate_into(row_marginal, x$lines$row, 1)
  row_marginal <- accumulate_into(row_marginal, x$points$row, 1)

  ## ---- column marginal: difference array for runs, O(records + ncol) -----
  d <- numeric(ncol + 1L)
  d <- accumulate_into(d, runs$col_start, 1)
  d <- accumulate_into(d, runs$col_end + 1L, -1)
  col_marginal <- cumsum(d)[seq_len(ncol)]
  col_marginal <- accumulate_into(col_marginal, x$edges$col, x$edges$fraction)
  col_marginal <- accumulate_into(col_marginal, x$lines$col, 1)
  col_marginal <- accumulate_into(col_marginal, x$points$col, 1)

  content <- content_extent(x)

  ## ---- per-geometry area -------------------------------------------------
  n <- attr(x, "n_geom")
  if (is.null(n)) {
    n <- max(0L, runs$id, x$edges$id, x$lines$id, x$points$id)
  }
  area <- data.frame(
    id = seq_len(n),
    interior = accumulate_into(numeric(n), runs$id, len),
    boundary = accumulate_into(numeric(n), x$edges$id, x$edges$fraction),
    lines = accumulate_into(numeric(n), x$lines$id, x$lines$length),
    points = accumulate_into(numeric(n), x$points$id, 1)
  )

  fractions <- if (nrow(x$edges)) {
    stats::quantile(x$edges$fraction, c(0, 0.25, 0.5, 0.75, 1), names = FALSE)
  } else rep(NA_real_, 5L)

  structure(
    list(
      extent = ext,
      dimension = dm,
      crs = attr(x, "crs"),
      coverage = attr(x, "coverage"),
      records = vapply(x, nrow, integer(1)),
      n_geom = n,
      n_noted = length(unique(x$notes$geom_index)),
      content = content,
      area = area,
      fractions = fractions,
      row_marginal = row_marginal,
      col_marginal = col_marginal,
      overview = cb_overview(x, tiles = tiles)
    ),
    class = "summary.cb_burn"
  )
}

#' Coarse overview raster of a burn
#'
#' Aggregate a [cb_burn()] onto a grid of square tiles of `block` cells and
#' return the mean coverage per tile (in `[0, 1]` for polygons) in the
#' `read_ds` layout used by [cb_materialize()]. Edge tiles that are cut off by
#' the grid boundary are normalized by their actual cell count. Cost is
#' proportional to the number of records plus the number of tiles: a run that
#' spans several tiles contributes partial lengths to its two end tiles and the
#' fully covered tiles between them are filled with a difference array, so
#' nothing is ever expanded to cells whatever the block size.
#'
#' @param x A `cb_burn` object.
#' @param tiles Target number of tiles along the longer axis; ignored if
#'   `block` is given.
#' @param block Tile size in cells. `block = 1` gives the exact coverage.
#'
#' @return A numeric vector with a `gis` attribute (`datatype = "Float32"`),
#'   and additionally attribute `block` giving the tile size used. Untouched
#'   tiles are `NA`.
#' @export
cb_overview <- function(x, tiles = 256L, block = NULL) {
  ext <- attr(x, "extent")
  dm <- as.integer(attr(x, "dimension"))
  ncol <- dm[1L]
  nrow <- dm[2L]
  b <- if (is.null(block)) max(1L, as.integer(ceiling(max(dm) / tiles))) else as.integer(block)
  ntx <- as.integer(ceiling(ncol / b))
  nty <- as.integer(ceiling(nrow / b))

  tile_index <- function(row, col) {
    ((row - 1L) %/% b) * ntx + ((col - 1L) %/% b) + 1L
  }

  ## runs: a run touches tile columns tc0..tc1 of tile row tr. Its two end
  ## tiles get partial lengths directly; every tile strictly between them is
  ## fully covered (b cells) and is handled with a difference array along the
  ## tile row, so the cost is O(runs + tiles) whatever the block size.
  runs <- x$runs
  tr <- (runs$row - 1L) %/% b
  tc0 <- (runs$col_start - 1L) %/% b
  tc1 <- (runs$col_end - 1L) %/% b
  base <- tr * ntx + 1L
  single <- tc0 == tc1

  acc <- numeric(ntx * nty)
  acc <- accumulate_into(acc, (base + tc0)[single],
                         (runs$col_end - runs$col_start + 1L)[single])
  multi <- !single
  acc <- accumulate_into(acc, (base + tc0)[multi],
                         ((tc0 + 1L) * b - runs$col_start + 1L)[multi])
  acc <- accumulate_into(acc, (base + tc1)[multi],
                         (runs$col_end - tc1 * b)[multi])
  span <- multi & (tc1 - tc0 > 1L)
  if (any(span)) {
    d <- numeric(ntx * nty + 1L)
    d <- accumulate_into(d, (base + tc0 + 1L)[span], b)
    d <- accumulate_into(d, (base + tc1)[span], -b)
    ## cumsum within each tile row: the tile vector is row-major so a tile row
    ## is a column of matrix(d, ntx)
    full <- matrix(d[seq_len(ntx * nty)], ntx, nty)
    full <- apply(full, 2L, cumsum)
    acc <- acc + as.vector(full)
  }
  acc <- accumulate_into(acc, tile_index(x$edges$row, x$edges$col), x$edges$fraction)
  acc <- accumulate_into(acc, tile_index(x$lines$row, x$lines$col), 1)
  acc <- accumulate_into(acc, tile_index(x$points$row, x$points$col), 1)

  ## actual cell count per tile (edge tiles are smaller)
  tw <- pmin(b, ncol - (seq_len(ntx) - 1L) * b)
  th <- pmin(b, nrow - (seq_len(nty) - 1L) * b)
  cells <- as.vector(t(outer(th, tw)))  # row-major: tile row varies slowest

  out <- acc / cells
  out[acc == 0] <- NA_real_

  crs <- attr(x, "crs")
  attr(out, "gis") <- list(
    type = "raster",
    ## the overview grid extends past the parent grid when b does not divide it
    bbox = c(ext[1L], ext[4L] - nty * b * (ext[4L] - ext[3L]) / nrow,
             ext[1L] + ntx * b * (ext[2L] - ext[1L]) / ncol, ext[4L]),
    dim = c(ntx, nty, 1),
    srs = if (is.null(crs) || is.na(crs)) "" else crs,
    datatype = "Float32"
  )
  attr(out, "block") <- b
  out
}

## Extent and row/col range of touched cells, or NULL for an empty burn.
content_extent <- function(x) {
  ext <- attr(x, "extent")
  dm <- attr(x, "dimension")
  runs <- x$runs
  rows <- c(runs$row, x$edges$row, x$lines$row, x$points$row)
  if (!length(rows)) return(NULL)
  cols <- c(runs$col_start, runs$col_end, x$edges$col, x$lines$col, x$points$col)
  rr <- range(rows)
  cc <- range(cols)
  resx <- (ext[2L] - ext[1L]) / dm[1L]
  resy <- (ext[4L] - ext[3L]) / dm[2L]
  list(
    rows = rr, cols = cc,
    extent = c(ext[1L] + (cc[1L] - 1L) * resx, ext[1L] + cc[2L] * resx,
               ext[4L] - rr[2L] * resy, ext[4L] - (rr[1L] - 1L) * resy)
  )
}

## out[idx] += w, summing duplicates in idx; w recycled to length(idx)
accumulate_into <- function(out, idx, w) {
  if (!length(idx)) return(out)
  s <- rowsum(rep_len(as.numeric(w), length(idx)), idx, reorder = FALSE)
  at <- as.integer(rownames(s))
  out[at] <- out[at] + s[, 1L]
  out
}

#' @export
print.summary.cb_burn <- function(x, ascii = FALSE, ...) {
  ext <- x$extent
  dm <- x$dimension
  cat("<cb_burn summary>\n")
  cat(sprintf("  grid:     %d col x %d row, %s mode\n", dm[1L], dm[2L],
              if (isTRUE(x$coverage)) "coverage" else "approx"))
  cat(sprintf("  extent:   x[%g, %g] y[%g, %g]\n", ext[1L], ext[2L], ext[3L], ext[4L]))
  if (is.null(x$content)) {
    cat("  content:  (empty)\n")
  } else {
    ce <- x$content$extent
    cat(sprintf("  content:  cols %d..%d, rows %d..%d  x[%g, %g] y[%g, %g]\n",
                x$content$cols[1L], x$content$cols[2L],
                x$content$rows[1L], x$content$rows[2L],
                ce[1L], ce[2L], ce[3L], ce[4L]))
  }
  cat("  records: ", paste(sprintf("%s=%d", names(x$records), x$records), collapse = " "), "\n")
  touched <- with(x$area, interior > 0 | boundary > 0 | lines > 0 | points > 0)
  cat(sprintf("  geoms:    %d input, %d burned, %d empty, %d with notes\n",
              x$n_geom, sum(touched), sum(!touched), x$n_noted))
  cat(sprintf("  area:     interior %s cells, boundary %s cells", fmt_num(sum(x$area$interior)),
              fmt_num(sum(x$area$boundary))))
  if (!anyNA(x$fractions)) {
    cat(sprintf(" (fraction min %.3g, median %.3g, max %.3g)",
                x$fractions[1L], x$fractions[3L], x$fractions[5L]))
  }
  cat("\n")
  if (sum(x$area$lines) > 0) cat(sprintf("  lines:    %s length units\n", fmt_num(sum(x$area$lines))))
  if (sum(x$area$points) > 0) cat(sprintf("  points:   %s in grid\n", fmt_num(sum(x$area$points))))
  ov <- attr(x$overview, "gis")$dim
  cat(sprintf("  overview: %d x %d tiles of %d cells\n", ov[1L], ov[2L], attr(x$overview, "block")))
  if (ascii) cat(ascii_overview(x$overview, ...), sep = "\n")
  invisible(x)
}

fmt_num <- function(v) format(signif(v, 4), big.mark = ",", scientific = FALSE, trim = TRUE)

## text rendering of an overview at most `width` characters wide
ascii_overview <- function(ov, width = 72L, ramp = " .:-=+*#%@") {
  m <- cb_as_matrix(ov)
  step <- max(1L, ceiling(ncol(m) / width))
  if (step > 1L) {
    ## characters are about twice as tall as wide: pool 2*step rows per step cols
    ci <- (seq_len(ncol(m)) - 1L) %/% step
    ri <- (seq_len(nrow(m)) - 1L) %/% (2L * step)
    m <- t(rowsum(t(rowsum(m, ri, na.rm = TRUE)), ci, na.rm = TRUE)) /
      outer(tabulate(ri + 1L), tabulate(ci + 1L))
  } else {
    ri <- (seq_len(nrow(m)) - 1L) %/% 2L
    m <- rowsum(m, ri, na.rm = TRUE) / tabulate(ri + 1L)
  }
  chars <- strsplit(ramp, "")[[1]]
  lev <- pmin(length(chars), 1L + floor(pmin(m, 1) * (length(chars) - 1L) + 0.5))
  lev[is.na(lev)] <- 1L
  apply(matrix(chars[lev], nrow(m)), 1L, paste, collapse = "")
}
