test_that("cb_crop snaps out and re-bases indices", {
  b <- square_burn()
  cr <- cb_crop(b, c(2.2, 5.1, 4.9, 9))
  expect_s3_class(cr, "cb_burn")
  expect_equal(attr(cr, "extent"), c(2, 6, 4, 9))
  expect_equal(attr(cr, "dimension"), c(4L, 5L))
  expect_equal(attr(cr, "parent")$offset, c(2L, 1L))
  expect_identical(attr(cr, "crs"), "EPSG:3857")
  m <- cb_as_matrix(cb_materialize(cr))
  full <- cb_as_matrix(cb_materialize(b))
  expect_equal(m, full[2:6, 3:6])
})

test_that("cb_crop clips runs at the window edge", {
  b <- square_burn()
  cr <- cb_crop(b, c(4, 6, 5, 7), snap = "in")
  expect_equal(attr(cr, "dimension"), c(2L, 2L))
  expect_equal(cr$runs$col_start, c(1L, 1L))
  expect_equal(cr$runs$col_end, c(2L, 2L))
  expect_equal(cr$runs$row, c(1L, 2L))
  expect_equal(nrow(cr$edges), 0L)
})

test_that("cb_crop intersects with the parent grid", {
  b <- square_burn()
  cr <- cb_crop(b, c(-5, 3, 5, 20))
  expect_equal(attr(cr, "extent"), c(0, 3, 5, 10))
  expect_equal(attr(cr, "dimension"), c(3L, 5L))
  expect_equal(attr(cr, "parent")$offset, c(0L, 0L))
})

test_that("cb_crop snap modes", {
  b <- square_burn()
  expect_equal(attr(cb_crop(b, c(2.4, 5.6, 2.4, 5.6), snap = "near"), "extent"),
               c(2, 6, 2, 6))
  expect_equal(attr(cb_crop(b, c(2.4, 5.6, 2.4, 5.6), snap = "in"), "extent"),
               c(3, 5, 3, 5))
  ## already-aligned extents are not pushed out by floating point noise
  expect_equal(attr(cb_crop(b, c(3, 5, 3, 5) + 1e-12), "extent"), c(3, 5, 3, 5))
  ## the full extent is a no-op crop
  cr <- cb_crop(b, attr(b, "extent"))
  expect_equal(cr$runs, b$runs)
  expect_equal(cr$edges, b$edges)
})

test_that("cb_crop errors on empty windows and bad input", {
  b <- square_burn()
  expect_error(cb_crop(b, c(20, 30, 0, 1)), "does not intersect")
  expect_error(cb_crop(b, c(3.2, 3.8, 3.2, 3.8), snap = "in"), "does not intersect")
  expect_error(cb_crop(b, c(0, 1)), "length-4")
  expect_error(cb_crop(list(), c(0, 1, 0, 1)), "cb_burn")
})

test_that("print reports the crop", {
  b <- square_burn()
  expect_output(print(cb_crop(b, c(2, 6, 4, 9))), "crop:")
  expect_output(print(b), "<cb_burn>")
})
