test_that("cb_materialize produces read_ds-shaped output", {
  b <- square_burn()
  v <- cb_materialize(b)
  expect_type(v, "double")
  expect_length(v, 100L)
  gis <- attr(v, "gis")
  expect_identical(gis$type, "raster")
  expect_equal(gis$bbox, c(0, 0, 10, 10))
  expect_equal(gis$dim, c(10, 10, 1))
  expect_identical(gis$srs, "EPSG:3857")
  expect_identical(gis$datatype, "Float32")
})

test_that("coverage values are laid out row-major from the top-left", {
  b <- square_burn()
  m <- cb_as_matrix(cb_materialize(b))
  expect_equal(dim(m), c(10L, 10L))
  expect_equal(sum(m, na.rm = TRUE), 16)
  expect_equal(m[4, 5], 1)
  expect_equal(m[2, 3], 0.25)
  expect_equal(m[2, 5], 0.5)
  expect_true(all(is.na(m[1, ])))
  expect_true(all(is.na(m[, 1:2])))
})

test_that("background and integer modes work", {
  b <- square_burn()
  v0 <- cb_materialize(b, background = 0)
  expect_equal(sum(v0), 16)

  vid <- cb_materialize(b, "id")
  expect_type(vid, "integer")
  expect_identical(attr(vid, "gis")$datatype, "Int32")
  expect_equal(sum(vid == 1L, na.rm = TRUE), 25L)

  vcount <- cb_materialize(b, "count", background = 0L)
  expect_equal(sum(vcount), 25L)
})

test_that("overlapping geometries sum, and highest id wins", {
  b <- square_burn()
  b$runs <- rbind(b$runs, data.frame(row = 4L, col_start = 5L, col_end = 5L, id = 2L))
  v <- cb_materialize(b)
  expect_equal(cb_as_matrix(v)[4, 5], 2)
  vid <- cb_as_matrix(cb_materialize(b, "id"))
  expect_equal(vid[4, 5], 2L)
  expect_equal(vid[4, 4], 1L)
  vcount <- cb_as_matrix(cb_materialize(b, "count"))
  expect_equal(vcount[4, 5], 2L)
  expect_equal(vcount[4, 4], 1L)
})

test_that("empty burn materializes to background", {
  b <- square_burn()
  b$runs <- b$runs[0, ]
  b$edges <- b$edges[0, ]
  v <- cb_materialize(b)
  expect_true(all(is.na(v)))
  expect_identical(attr(v, "gis")$srs, "EPSG:3857")
  attr(b, "crs") <- NA_character_
  expect_identical(attr(cb_materialize(b), "gis")$srs, "")
})

test_that("cb_as_matrix requires a gis attribute", {
  expect_error(cb_as_matrix(1:4), "gis")
})
