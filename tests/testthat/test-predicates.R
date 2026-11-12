## One burn holding the awkward cases: the crate book's polygon (edges on
## walls and through corners), a neighbour sharing a wall with it, a polygon
## with a hole, a polygon inside one cell, a line (no polygon matrix), and a
## second polygon nested in the first's area.
fixture <- c(
  "POLYGON ((1 1, 6 1, 6 3, 8.5 5.5, 4 7, 1 4, 1 1))",
  "POLYGON ((6 1, 9 1, 9 3, 6 3, 6 1))",
  "POLYGON ((7 5.5, 9.5 5.5, 9.5 7.5, 7 7.5, 7 5.5), (7.5 6, 8 6, 8 7, 7.5 7, 7.5 6))",
  "POLYGON ((0.2 0.2, 0.8 0.2, 0.8 0.6, 0.2 0.2))",
  "LINESTRING (0.5 7.5, 9.5 0.5)",
  "POLYGON ((2 2, 4 2, 4 5, 2 5, 2 2))"
)
fx_extent <- c(0, 10, 0, 8)
fx_dim <- c(10L, 8L)
fx_polys <- c(1L, 2L, 3L, 4L, 6L)
fx_burn <- function(...) {
  cb_burn(wk::wkt(fixture), extent = fx_extent, dimension = fx_dim, relate = TRUE, ...)
}

## every cell, row-major from the top-left, as the oracle orders them
all_rows <- function(dm) rep(seq_len(dm[2]), each = dm[1])
all_cols <- function(dm) rep(seq_len(dm[1]), dm[2])

## cell outlines in world coordinates, row-major from the top-left
cell_sfc <- function(gt, dm) {
  px <- function(col, row) c(gt[1] + col * gt[2] + row * gt[3], gt[4] + col * gt[5] + row * gt[6])
  sf::st_sfc(lapply(seq_len(prod(dm)) - 1L, function(i) {
    r <- i %/% dm[1]; c <- i %% dm[1]
    sf::st_polygon(list(rbind(px(c, r), px(c + 1, r), px(c + 1, r + 1), px(c, r + 1), px(c, r))))
  }))
}

## "id row col" keys of the cells where a GEOS predicate holds, per polygon
geos_cells <- function(wkts, ids, gt, dm, predicate) {
  cells <- cell_sfc(gt, dm)
  polys <- sf::st_as_sfc(wkts)
  hit <- switch(predicate,
    intersects = sf::st_intersects(polys, cells, sparse = FALSE),
    covers = sf::st_covers(polys, cells, sparse = FALSE),
    touches = sf::st_touches(polys, cells, sparse = FALSE),
    covered_by = sf::st_covered_by(polys, cells, sparse = FALSE),
    interiors_meet = sf::st_relate(polys, cells, pattern = "T********", sparse = FALSE)
  )
  hit_keys(hit, ids, dm)
}

## "id row col" for the TRUE cells of an ids x cells matrix; paste() would
## turn a row with no hits into a bogus "id  " key, so skip those
hit_keys <- function(hit, ids, dm) {
  rows <- all_rows(dm); cols <- all_cols(dm)
  unlist(lapply(seq_along(ids), function(k) {
    w <- which(hit[k, ])
    if (length(w)) paste(ids[k], rows[w], cols[w]) else character()
  }))
}
keys <- function(cells) paste(cells$id, cells$row, cells$col)

test_that("cb_matrix gives GEOS's matrix for every cell and polygon", {
  skip_if_not_installed("sf")
  b <- fx_burn()
  gt <- north_up(fx_extent, fx_dim)
  for (id in fx_polys) {
    got <- cb_matrix(b, all_rows(fx_dim), all_cols(fx_dim), id)
    expect_identical(got, oracle_grid(fixture[id], gt, fx_dim), info = paste("id", id))
  }
  ## a line has no polygon matrix
  expect_true(all(is.na(cb_matrix(b, all_rows(fx_dim), all_cols(fx_dim), 5L))))
})

test_that("cb_matrix is exact where the burn's records are not", {
  skip_if_not_installed("sf")
  slivers <- c(
    "POLYGON ((1 1, 3 1, 3 2, 3.0000001 2.5, 3 3, 1 3, 1 1))",
    "POLYGON ((1 1, 3 1, 3 3, 2.6 3, 2.5 2.9999999, 2.4 3, 1 3, 1 1))"
  )
  gt <- north_up(c(0, 4, 0, 4), c(4L, 4L))
  b <- cb_burn(wk::wkt(slivers), extent = c(0, 4, 0, 4), dimension = c(4, 4), relate = TRUE)
  for (id in 1:2) {
    expect_identical(cb_matrix(b, all_rows(c(4, 4)), all_cols(c(4, 4)), id),
                     oracle_grid(slivers[id], gt, c(4L, 4L)), info = paste("id", id))
  }
})

test_that("cb_matrix matches GEOS on random polygons, with an id offset", {
  skip_if_not_installed("sf")
  set.seed(5501)
  dm <- 16L
  wkts <- character()
  while (length(wkts) < 12L) {
    w <- random_star(dm, snap = c(0, 0.5, 1)[length(wkts) %% 3L + 1L], hole = length(wkts) %% 2L == 0L)
    if (!is.null(w) && isTRUE(sf::st_is_valid(sf::st_as_sfc(w)))) wkts <- c(wkts, w)
  }
  b <- cb_burn(wk::wkt(wkts), extent = c(0, dm, 0, dm), dimension = c(dm, dm),
               relate = TRUE, id_offset = 40)
  gt <- north_up(c(0, dm, 0, dm), c(dm, dm))
  for (k in seq_along(wkts)) {
    expect_identical(cb_matrix(b, all_rows(c(dm, dm)), all_cols(c(dm, dm)), 40L + k),
                     oracle_grid(wkts[k], gt, c(dm, dm)), info = paste("polygon", k))
  }
})

test_that("cb_matrix works on rotated grids and cropped burns", {
  skip_if_not_installed("sf")
  gt <- cb_geotransform(c(100, 200), res = c(2, 1.5), rotation = 30, shear = 10)
  dm <- c(8L, 6L)
  px <- function(col, row) c(gt[1] + col * gt[2] + row * gt[3], gt[4] + col * gt[5] + row * gt[6])
  a <- seq(0, 2 * pi, length.out = 10)[-1]
  set.seed(8812)
  ring <- t(vapply(a, function(t) px(4, 3) + runif(1, 2.5, 7) * c(cos(t), sin(t)), numeric(2)))
  ring <- rbind(ring, ring[1, ])
  w <- sprintf("POLYGON ((%s))", paste(sprintf("%.17g %.17g", ring[, 1], ring[, 2]), collapse = ", "))
  b <- cb_burn_index(wk::wkt(w), gt, dm, relate = TRUE)
  expect_identical(cb_matrix(b, all_rows(dm), all_cols(dm)), oracle_grid(w, gt, dm))

  full <- fx_burn()
  cr <- cb_crop(full, c(2, 8, 1, 6))
  off <- attr(cr, "parent")$offset
  dmc <- attr(cr, "dimension")
  for (id in fx_polys) {
    expect_identical(cb_matrix(cr, all_rows(dmc), all_cols(dmc), id),
                     cb_matrix(full, all_rows(dmc) + off[2], all_cols(dmc) + off[1], id))
  }
})

test_that("cb_matrix recycles, returns NA off the grid, and checks its input", {
  b <- fx_burn()
  m <- cb_matrix(b, row = c(0, 9, 1.5, 4, 4, NA), col = c(1, 1, 1, 11, 4, 4), id = 1)
  expect_equal(m[1:4], rep(NA_character_, 4))
  expect_false(is.na(m[5]))
  expect_true(is.na(m[6]))
  expect_true(all(is.na(cb_matrix(b, 4, 4, c(0, 7, 99)))))
  expect_identical(cb_matrix(b, 4, 1:3, 1), cb_matrix(b, c(4, 4, 4), 1:3, c(1, 1, 1)))
  expect_identical(cb_matrix(b, integer(), 1, 1), character())
  expect_error(cb_matrix(b, 1:2, 1:3, 1), "recycle")
  expect_error(cb_matrix(b, 1, 1), "`id` is needed")
  single <- cb_burn(wk::wkt(fixture[1]), extent = fx_extent, dimension = fx_dim, relate = TRUE)
  expect_identical(cb_matrix(single, 1:8, 3), cb_matrix(b, 1:8, 3, 1))
  expect_error(cb_matrix(cb_burn(wk::wkt(fixture[1]), extent = fx_extent, dimension = fx_dim), 1, 1),
               "relations")
})

test_that("cb_cells_where selects exactly the cells GEOS's predicates hold for", {
  skip_if_not_installed("sf")
  b <- fx_burn()
  gt <- north_up(fx_extent, fx_dim)
  for (p in c("intersects", "interiors_meet", "covers", "touches", "covered_by")) {
    got <- cb_cells_where(b, p)
    expect_setequal(keys(got), geos_cells(fixture[fx_polys], fx_polys, gt, fx_dim, p))
    expect_identical(got$matrix, cb_matrix(b, got$row, got$col, got$id), info = p)
    expect_identical(order(got$id, got$row, got$col), seq_len(nrow(got)))
  }
  ## the polygon inside one cell is the only one a cell covers
  expect_equal(unique(cb_cells_where(b, "covered_by")$id), 4L)
  ## the two neighbours share a wall: each touches the other's cells along it
  touch <- cb_cells_where(b, "touches", id = 2L)
  expect_true(all(touch$id == 2L))
  expect_true(any(touch$col == 6L))
})

test_that("cb_cells_where takes DE-9IM patterns, as GEOS does", {
  skip_if_not_installed("sf")
  b <- fx_burn()
  gt <- north_up(fx_extent, fx_dim)
  cells <- cell_sfc(gt, fx_dim)
  polys <- sf::st_as_sfc(fixture[fx_polys])
  for (pat in c("*T*******", "F***1****", "T*F**F***", "2121F1212", "****0****")) {
    hit <- sf::st_relate(polys, cells, pattern = pat, sparse = FALSE)
    expect_setequal(keys(cb_cells_where(b, pattern = pat)), hit_keys(hit, fx_polys, fx_dim))
  }
  expect_identical(cb_cells_where(b, pattern = "t*f**f***"), cb_cells_where(b, pattern = "T*F**F***"))
  expect_error(cb_cells_where(b, pattern = "FF*******"), "does not touch")
  expect_error(cb_cells_where(b, pattern = "T*F"), "nine characters")
  expect_error(cb_cells_where(b, pattern = "T*F**F**X"), "nine characters")
  expect_error(cb_cells_where(b, "nearby"), "should be one of")
})

test_that("predicate masks count polygons per cell, and intersects is exact all-touched", {
  skip_if_not_installed("sf")
  b <- fx_burn()
  gt <- north_up(fx_extent, fx_dim)
  cells <- cell_sfc(gt, fx_dim)
  polys <- sf::st_as_sfc(fixture[fx_polys])
  oracle <- list(
    intersects = sf::st_intersects(polys, cells, sparse = FALSE),
    covers = sf::st_covers(polys, cells, sparse = FALSE),
    touches = sf::st_touches(polys, cells, sparse = FALSE),
    covered_by = sf::st_covered_by(polys, cells, sparse = FALSE),
    interiors_meet = sf::st_relate(polys, cells, pattern = "T********", sparse = FALSE)
  )
  for (p in names(oracle)) {
    v <- cb_materialize(b, p, background = 0L)
    expect_identical(as.vector(v), as.integer(colSums(oracle[[p]])), info = p)
    expect_equal(attr(v, "gis")$datatype, "Int32")
    expect_equal(attr(v, "gis")$dim, c(10, 8, 1))
    na <- cb_materialize(b, p)
    expect_identical(is.na(as.vector(na)), colSums(oracle[[p]]) == 0)
  }

  ## the crate book's polygon: 54 cells touched, 34 with a record
  one <- cb_burn(wk::wkt(fixture[1]), extent = fx_extent, dimension = fx_dim, relate = TRUE)
  expect_equal(sum(!is.na(cb_materialize(one, "intersects"))), 54L)
  expect_equal(sum(!is.na(cb_materialize(one, "count"))), 34L)

  expect_error(cb_materialize(cb_burn(wk::wkt(fixture[1]), extent = fx_extent, dimension = fx_dim),
                              "intersects"), "relations")
})

test_that("predicate masks keep an affine grid's georeferencing", {
  gt <- cb_geotransform(c(100, 200), res = 2, rotation = 20)
  b <- cb_burn_index(wk::wkt("POLYGON ((104 196, 110 196, 110 190, 104 196))"), gt, c(8, 8),
                     relate = TRUE, crs = "EPSG:32755")
  gis <- attr(cb_materialize(b, "intersects"), "gis")
  want <- attr(cb_materialize(b), "gis")
  want$datatype <- "Int32"
  expect_equal(gis, want)
  expect_equal(gis$geotransform, gt)
})

test_that("cb_slivers finds the cells the burn's tolerance rounds away", {
  skip_if_not_installed("sf")
  slivers <- c(
    spike = "POLYGON ((1 1, 3 1, 3 2, 3.0000001 2.5, 3 3, 1 3, 1 1))",
    notch = "POLYGON ((1 1, 3 1, 3 3, 2.6 3, 2.5 2.9999999, 2.4 3, 1 3, 1 1))",
    plain = "POLYGON ((0.5 0.5, 3.5 0.5, 3.5 1.5, 0.5 0.5))"
  )
  b <- cb_burn(wk::wkt(slivers), extent = c(0, 4, 0, 4), dimension = c(4, 4), relate = TRUE)
  s <- cb_slivers(b)
  expect_named(s, c("row", "col", "id", "matrix", "kind"))
  expect_equal(s$id, c(1L, 2L))
  expect_equal(s$kind, c("interior", "exterior"))
  ## the spike's cell has no record, and the notch's cell is in a run
  expect_equal(c(s$row[1], s$col[1]), c(2L, 4L))
  expect_equal(c(s$row[2], s$col[2]), c(2L, 3L))
  expect_false(any(b$edges$id == 1L & b$edges$row == 2L & b$edges$col == 4L))
  expect_equal(cb_matrix(b, 2, 3, 2), s$matrix[2])

  ## sf agrees that the area involved is real but below the burn's 1e-6
  gt <- north_up(c(0, 4, 0, 4), c(4L, 4L))
  cells <- cell_sfc(gt, c(4L, 4L))
  area_in <- function(w, row, col) {
    i <- sf::st_intersection(sf::st_as_sfc(w), cells[(row - 1) * 4 + col])
    if (length(i)) as.numeric(sf::st_area(i)) else 0
  }
  inside <- area_in(slivers[["spike"]], 2, 4)
  expect_gt(inside, 0)
  expect_lt(inside, 1e-6)
  outside <- 1 - area_in(slivers[["notch"]], 2, 3)
  expect_gt(outside, 0)
  expect_lt(outside, 1e-6)

  ## records agree with the topology for ordinary polygons
  expect_equal(nrow(cb_slivers(fx_burn())), 0L)
  expect_error(cb_slivers(cb_burn(wk::wkt(slivers), extent = c(0, 4, 0, 4), dimension = c(4, 4))),
               "relations")
})
