## Shared by the relate and predicate tests: the matrix of every cell from a
## burn, and the GEOS oracle (sf::st_relate) on each cell outline.

## The matrix of every cell of one geometry: the relations where there are
## any, a run cell is covered, everything else is disjoint
matrix_grid <- function(b, id = 1L) {
  dm <- attr(b, "dimension")
  out <- rep("FF2FF1212", prod(dm))
  r <- b$runs[b$runs$id == id, ]
  len <- r$col_end - r$col_start + 1L
  out[(rep.int(r$row, len) - 1L) * dm[1] + sequence(len, from = r$col_start)] <- "212FF1FF2"
  rel <- b$relations[b$relations$id == id, ]
  out[(rel$row - 1L) * dm[1] + rel$col] <- rel$matrix
  out
}

## the oracle: sf (GEOS) relates the polygon to each cell outline, in
## world coordinates, row-major from the top-left cell
oracle_grid <- function(wkt, gt, dm) {
  p <- sf::st_as_sfc(wkt)
  px <- function(col, row) c(gt[1] + col * gt[2] + row * gt[3], gt[4] + col * gt[5] + row * gt[6])
  cells <- lapply(seq_len(prod(dm)) - 1L, function(i) {
    r <- i %/% dm[1]; c <- i %% dm[1]
    sf::st_polygon(list(rbind(px(c, r), px(c + 1, r), px(c + 1, r + 1), px(c, r + 1), px(c, r))))
  })
  as.vector(sf::st_relate(p, sf::st_sfc(cells)))
}
north_up <- function(extent, dm) {
  c(extent[1], (extent[2] - extent[1]) / dm[1], 0, extent[4], 0, -(extent[4] - extent[3]) / dm[2])
}

## entries of a matrix string: polygon (I, B, E) against cell (I, B, E)
entry <- function(m, i) substr(m, i, i)
covers <- function(m) entry(m, 7) == "F" & entry(m, 8) == "F"
touches <- function(m) entry(m, 1) == "F" & (entry(m, 2) != "F" | entry(m, 4) != "F" | entry(m, 5) != "F")
interiors_meet <- function(m) entry(m, 1) != "F"

## a star about (cx, cy), optionally with a hole, vertices snapped to `snap`
## (0: free floats)
random_star <- function(dm, snap, hole) {
  cx <- runif(1, 3, dm - 3); cy <- runif(1, 3, dm - 3)
  rad <- runif(1, 1.2, 3)
  n <- sample(8:20, 1)
  a <- seq(0, 2 * pi, length.out = n + 1)[-1]
  r <- rad * (0.4 + 0.6 * runif(n))
  ring <- cbind(cx + r * cos(a), cy + r * sin(a))
  snapit <- function(m) if (snap > 0) round(m / snap) * snap else m
  close <- function(m) {
    m <- m[c(TRUE, rowSums(abs(diff(m))) > 0), , drop = FALSE]
    rbind(m, m[1, ])
  }
  ring <- close(snapit(ring))
  rings <- list(ring)
  if (hole) {
    ha <- seq(0, 2 * pi, length.out = 13)[-1]
    h <- cbind(cx + 0.3 * rad * cos(ha), cy + 0.3 * rad * sin(ha))
    rings[[2]] <- close(snapit(h))
  }
  if (any(vapply(rings, nrow, 1L) < 4L)) return(NULL)
  txt <- vapply(rings, function(m) {
    paste0("(", paste(sprintf("%.17g %.17g", m[, 1], m[, 2]), collapse = ", "), ")")
  }, "")
  sprintf("POLYGON (%s)", paste(txt, collapse = ", "))
}
