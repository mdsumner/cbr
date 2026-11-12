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

test_that("matrices of fixed edge cases match GEOS, cell by cell", {
  skip_if_not_installed("sf")
  cases <- c(
    wall_square = "POLYGON ((1 1, 3 1, 3 3, 1 3, 1 1))",
    equal_to_a_cell = "POLYGON ((2 2, 3 2, 3 3, 2 3, 2 2))",
    hole_equals_a_cell = "POLYGON ((0 0, 6 0, 6 6, 0 6, 0 0), (2 2, 3 2, 3 3, 2 3, 2 2))",
    inside_one_cell = "POLYGON ((2.2 2.2, 2.8 2.2, 2.8 2.8, 2.2 2.8, 2.2 2.2))",
    inside_one_cell_on_walls = "POLYGON ((2 2.5, 2.5 2, 3 2.5, 2.5 3, 2 2.5))",
    diamond_corners_on_walls = "POLYGON ((3 1, 5 3, 3 5, 1 3, 3 1))",
    diagonal_through_corners = "POLYGON ((1 1, 4 4, 1 4, 1 1))",
    crossing_the_extent = "POLYGON ((-2 3.5, 2.5 3.5, 2.5 5.5, -2 5.5, -2 3.5))",
    along_the_extent = "POLYGON ((0 0, 2 0, 2 2, 0 2, 0 0))",
    l_shape = "POLYGON ((1 1, 5 1, 5 3, 3 3, 3 5, 1 5, 1 1))",
    hole_touching_the_shell = "POLYGON ((1 1, 6 1, 6 6, 1 6, 1 1), (1 3, 3 2, 3 4, 1 3))",
    multipolygon_at_a_point = "MULTIPOLYGON (((1 1, 3 1, 3 3, 1 3, 1 1)), ((3 3, 5 3, 5 5, 3 5, 3 3)))"
  )
  gt <- north_up(c(0, 8, 0, 8), c(8, 8))
  for (nm in names(cases)) {
    b <- cb_burn(wk::wkt(cases[[nm]]), extent = c(0, 8, 0, 8), dimension = c(8, 8),
                 relate = TRUE)
    expect_identical(matrix_grid(b), oracle_grid(cases[[nm]], gt, c(8, 8)), info = nm)
  }
})

test_that("the cell-centre mode gives the same matrices", {
  skip_if_not_installed("sf")
  w <- "POLYGON ((1 1, 5 1, 5 3, 3 3, 3 5, 1 5, 1 1))"
  b <- cb_burn(wk::wkt(w), extent = c(0, 8, 0, 8), dimension = c(8, 8),
               coverage = FALSE, relate = TRUE)
  expect_identical(matrix_grid(b), oracle_grid(w, north_up(c(0, 8, 0, 8), c(8, 8)), c(8, 8)))
})

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

test_that("random polygons with holes match GEOS, free and on a lattice", {
  skip_if_not_installed("sf")
  set.seed(7319)
  dm <- 16L
  gt <- north_up(c(0, dm, 0, dm), c(dm, dm))
  for (snap in c(0, 0.5, 1)) {
    done <- 0L
    tries <- 0L
    while (done < 25L && tries < 400L) {
      tries <- tries + 1L
      w <- random_star(dm, snap, hole = done %% 2L == 0L)
      if (is.null(w) || !isTRUE(sf::st_is_valid(sf::st_as_sfc(w)))) next
      b <- cb_burn(wk::wkt(w), extent = c(0, dm, 0, dm), dimension = c(dm, dm), relate = TRUE)
      expect_identical(matrix_grid(b), oracle_grid(w, gt, c(dm, dm)),
                       info = sprintf("snap %g, polygon %d", snap, done))
      done <- done + 1L
    }
    expect_gte(done, 15L)
  }
})

test_that("a rotated, sheared grid matches GEOS in world coordinates", {
  skip_if_not_installed("sf")
  gt <- cb_geotransform(c(100, 200), res = c(2, 1.5), rotation = 30, shear = 10)
  dm <- c(8L, 6L)
  px <- function(col, row) c(gt[1] + col * gt[2] + row * gt[3], gt[4] + col * gt[5] + row * gt[6])
  centre <- px(4, 3)
  set.seed(2241)
  for (k in 1:6) {
    n <- 9
    a <- seq(0, 2 * pi, length.out = n + 1)[-1]
    r <- runif(n, 2.5, 7)
    ring <- cbind(centre[1] + r * cos(a), centre[2] + r * sin(a))
    ring <- rbind(ring, ring[1, ])
    w <- sprintf("POLYGON ((%s))", paste(sprintf("%.17g %.17g", ring[, 1], ring[, 2]), collapse = ", "))
    skip_if_not(isTRUE(sf::st_is_valid(sf::st_as_sfc(w))))
    b <- cb_burn_index(wk::wkt(w), gt, dm, relate = TRUE)
    expect_identical(matrix_grid(b), oracle_grid(w, gt, dm), info = paste("polygon", k))
    expect_equal(cb_relate(cb_burn_index(wk::wkt(w), gt, dm), wk::wkt(w)), b)
  }
})

test_that("the crate's documented example is reproduced", {
  ring <- rbind(c(1, 1), c(6, 1), c(6, 3), c(8.5, 5.5), c(4, 7), c(1, 4), c(1, 1))
  g <- wk::wkt(sprintf("POLYGON ((%s))", paste(sprintf("%g %g", ring[, 1], ring[, 2]), collapse = ", ")))
  b <- cb_burn(g, extent = c(0, 10, 0, 8), dimension = c(10, 8), relate = TRUE)
  m <- b$relations$matrix
  expect_equal(nrow(b$relations), 43L)
  expect_equal(c(sum(covers(m)), sum(interiors_meet(m) & !covers(m)), sum(touches(m))),
               c(13L, 10L, 20L))
  ## cells are (row, col) 1-based here: the crate's (7, 2), (6, 2), (0, 3)
  at <- function(row, col) b$relations$matrix[b$relations$row == row & b$relations$col == col]
  expect_equal(at(8, 3), "FF2F11212")
  expect_equal(at(7, 3), "212F11FF2")
  expect_equal(at(1, 4), "FF2F01212")
  ## the exact all-touched set is the runs plus the cells that intersect:
  ## 11 run cells the boundary does not reach, plus the 43 relations
  all_touched <- unique(c(
    paste(rep.int(b$runs$row, b$runs$col_end - b$runs$col_start + 1L),
          sequence(b$runs$col_end - b$runs$col_start + 1L, from = b$runs$col_start)),
    paste(b$relations$row, b$relations$col)
  ))
  expect_length(all_touched, 54L)
})

test_that("relations come from the burn call or from cb_relate, identically", {
  g <- wk::wkt(c("POLYGON ((1.5 1.5, 6.5 2.5, 5.5 6.5, 2 5, 1.5 1.5))",
                 "LINESTRING (0.5 7.5, 7.5 0.5)",
                 "POLYGON ((5 5, 7 5, 7 7, 5 7, 5 5), (5.5 5.5, 6 5.5, 6 6, 5.5 6, 5.5 5.5))",
                 "POINT (3.5 3.5)"))
  args <- list(g, extent = c(0, 8, 0, 8), dimension = c(8, 8))
  for (extra in list(list(), list(id_offset = 10), list(coverage = FALSE), list(threads = 2))) {
    onecall <- do.call(cb_burn, c(args, extra, relate = TRUE))
    posthoc <- cb_relate(do.call(cb_burn, c(args, extra)), g)
    expect_equal(onecall, posthoc)
    ## only the two polygons have relations, ids offset like the other tables
    off <- if (is.null(extra$id_offset)) 0L else extra$id_offset
    expect_setequal(unique(onecall$relations$id), off + c(1L, 3L))
    expect_identical(order(onecall$relations$id, onecall$relations$row, onecall$relations$col),
                     seq_len(nrow(onecall$relations)))
  }
  plain <- do.call(cb_burn, args)
  expect_null(plain$relations)
  expect_named(plain, c("runs", "edges", "lines", "points", "notes"))
})

test_that("empty, unparseable and non-polygon input has no relations", {
  good <- unclass(wk::as_wkb(wk::wkt("POLYGON ((1 1, 4 1, 4 4, 1 4, 1 1))")))[[1L]]
  blobs <- structure(list(good, as.raw(c(1, 2, 3)), raw()), class = "wk_wkb")
  b <- cb_burn(blobs, extent = c(0, 5, 0, 5), dimension = c(5, 5), relate = TRUE)
  expect_setequal(unique(b$relations$id), 1L)
  expect_equal(b$notes$geom_index, 2L)
  none <- cb_burn(wk::wkt("POINT (2.5 2.5)"), extent = c(0, 5, 0, 5), dimension = c(5, 5),
                  relate = TRUE)
  expect_equal(nrow(none$relations), 0L)
  expect_named(none$relations, c("row", "col", "id", "matrix"))
})

test_that("cropping crops the relations along with the other tables", {
  g <- wk::wkt(ring <- "POLYGON ((2.5 2.5, 9.5 3.5, 8.5 9.5, 3 8, 2.5 2.5))")
  b <- cb_burn(g, extent = c(0, 12, 0, 12), dimension = c(12, 12), relate = TRUE)
  cr <- cb_crop(b, c(2, 9, 3, 10))
  off <- attr(cr, "parent")$offset
  dm <- attr(cr, "dimension")
  expect_true(nrow(cr$relations) > 0L)
  expect_true(nrow(cr$relations) < nrow(b$relations))
  key <- function(t, dr = 0L, dc = 0L) paste(t$id, t$row + dr, t$col + dc)
  expect_true(all(key(cr$relations, off[2], off[1]) %in% key(b$relations)))
  expect_equal(cr$relations$matrix,
               b$relations$matrix[match(key(cr$relations, off[2], off[1]), key(b$relations))])
  inside <- b$relations$row > off[2] & b$relations$row <= off[2] + dm[2] &
    b$relations$col > off[1] & b$relations$col <= off[1] + dm[1]
  expect_equal(nrow(cr$relations), sum(inside))
  ## a burn without relations stays without
  expect_null(cb_crop(cb_burn(g, extent = c(0, 12, 0, 12), dimension = c(12, 12)),
                      c(2, 9, 3, 10))$relations)
})

test_that("combining chunks with relations equals relating the whole", {
  g <- wk::wkt(c("POLYGON ((1.5 1.5, 6.5 2.5, 5.5 6.5, 2 5, 1.5 1.5))",
                 "POLYGON ((5 5, 7 5, 7 7, 5 7, 5 5))",
                 "MULTIPOLYGON (((0.5 0.5, 1.5 0.5, 1.5 1.5, 0.5 1.5, 0.5 0.5)), ((6 0.5, 7.5 1, 7 2, 6 0.5)))"))
  args <- list(extent = c(0, 8, 0, 8), dimension = c(8, 8), relate = TRUE)
  whole <- do.call(cb_burn, c(list(g), args))
  a <- do.call(cb_burn, c(list(g[1]), args))
  b <- do.call(cb_burn, c(list(g[2:3]), args, id_offset = 1))
  expect_equal(cb_combine(b, a), whole)

  without <- cb_burn(g[2:3], extent = c(0, 8, 0, 8), dimension = c(8, 8), id_offset = 1)
  expect_error(cb_combine(a, without), "relations")
})

test_that("cb_relate refuses what it cannot relate", {
  g <- wk::wkt("POLYGON ((1 1, 4 1, 4 4, 1 4, 1 1))")
  b <- cb_burn(g, extent = c(0, 5, 0, 5), dimension = c(5, 5))
  expect_error(cb_relate(1, g), "cb_burn")
  expect_error(cb_relate(b, wk::wkt(c("POINT (1 1)", "POINT (2 2)"))), "geometries")
  expect_error(cb_relate(cb_crop(b, c(1, 4, 1, 4)), g), "cropped")
})

test_that("print and summary cope with the relations table", {
  g <- wk::wkt("POLYGON ((1.5 1.5, 4.5 1.5, 4.5 4.5, 1.5 1.5))")
  b <- cb_burn(g, extent = c(0, 6, 0, 6), dimension = c(6, 6), relate = TRUE)
  expect_output(print(b), "relations=")
  expect_s3_class(summary(b), "summary.cb_burn")
  expect_output(print(summary(b)), "relations=")
  expect_equal(sum(cb_materialize(b, background = 0)),
               sum(cb_materialize(cb_burn(g, extent = c(0, 6, 0, 6), dimension = c(6, 6)),
                                  background = 0)))
})
