## Integration tests against the compiled crate: these pin the R-side index
## convention (1-based row/col, inclusive col_end, id = 1-based input
## position) regardless of what the crate emits internally.

sq <- wk::wkt("POLYGON ((2.5 4.5, 6.5 4.5, 6.5 8.5, 2.5 8.5, 2.5 4.5))")
grid <- list(extent = c(0, 10, 0, 10), dimension = c(10L, 10L))

test_that("polygon runs and edges are 1-based, row 1 on top, col_end inclusive", {
  b <- cb_burn(sq, extent = grid$extent, dimension = grid$dimension)
  expect_identical(b$runs$row, 3:5)
  expect_true(all(b$runs$col_start == 4L))
  expect_true(all(b$runs$col_end == 6L))
  expect_true(all(b$runs$id == 1L))
  expect_equal(range(b$edges$row), c(2L, 6L))
  expect_equal(range(b$edges$col), c(3L, 7L))
  expect_equal(sum(b$edges$fraction), 7, tolerance = 1e-6)
  expect_equal(attr(b, "n_geom"), 1L)
  ## agrees cell for cell with the hand-built fixture used by the other tests
  expect_equal(cb_materialize(b), cb_materialize(square_burn(crs = NA_character_)))
})

test_that("a run reaching the last column ends at ncol", {
  full <- wk::wkt("POLYGON ((0 0, 10 0, 10 10, 0 10, 0 0))")
  b <- cb_burn(full, extent = grid$extent, dimension = grid$dimension)
  expect_identical(b$runs$row, 1:10)
  expect_true(all(b$runs$col_start == 1L))
  expect_true(all(b$runs$col_end == 10L))
  expect_equal(nrow(b$edges), 0L)
  expect_equal(sum(cb_materialize(b)), 100)
})

test_that("points, lines, ids and notes are 1-based", {
  pt <- wk::wkt("POINT (0.5 9.5)")
  ln <- wk::wkt("LINESTRING (0 0.5, 10 0.5)")
  blobs <- c(unclass(wk::as_wkb(sq)), unclass(wk::as_wkb(pt)),
             unclass(wk::as_wkb(ln)), list(as.raw(c(1, 2, 3))))
  b <- cb_burn(wk::new_wk_wkb(blobs), extent = grid$extent, dimension = grid$dimension)
  expect_equal(attr(b, "n_geom"), 4L)
  expect_equal(b$points, data.frame(row = 1L, col = 1L, id = 2L))
  expect_true(all(b$lines$row == 10L))
  expect_identical(sort(b$lines$col), 1:10)
  expect_true(all(b$lines$id == 3L))
  expect_equal(sum(b$lines$length), 10, tolerance = 1e-5)
  expect_equal(b$notes$geom_index, 4L)
  s <- summary(b)
  expect_equal(s$area$id, 1:4)
  expect_equal(s$area$interior, c(9, 0, 0, 0))
  expect_equal(s$area$points, c(0, 1, 0, 0))
})

test_that("approx mode emits runs only", {
  b <- cb_burn(sq, extent = grid$extent, dimension = grid$dimension, coverage = FALSE)
  expect_equal(nrow(b$edges), 0L)
  expect_true(all(b$runs$col_start == 3L))
  expect_true(all(b$runs$col_end == 6L))
  expect_equal(nrow(b$runs), 4L)
  expect_false(attr(b, "coverage"))
})

test_that("defaults burn onto the geometry bbox", {
  b <- cb_burn(sq)
  expect_equal(attr(b, "extent"), c(2.5, 6.5, 4.5, 8.5))
  expect_equal(attr(b, "dimension"), c(256L, 256L))
  expect_equal(sum(cb_materialize(b, background = 0)), 256^2)
})
