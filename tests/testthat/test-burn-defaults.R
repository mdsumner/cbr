test_that("default_nrow gives square cells and at least one row", {
  expect_equal(default_nrow(c(0, 10, 0, 10), 256L), 256L)
  expect_equal(default_nrow(c(0, 10, 0, 5), 256L), 128L)
  expect_equal(default_nrow(c(0, 10, 0, 5.03), 256L), 129L)
  expect_equal(default_nrow(c(0, 1000, 0, 1), 256L), 1L)
  expect_type(default_nrow(c(0, 10, 0, 10), 256), "integer")
})

test_that("default extent follows wk_bbox in terra ordering", {
  sq <- wk::wkt("POLYGON ((2.5 4.5, 6.5 4.5, 6.5 8.5, 2.5 8.5, 2.5 4.5))")
  ext <- as.numeric(wk::wk_bbox(wk::as_wkb(sq)))[c(1L, 3L, 2L, 4L)]
  expect_equal(ext, c(2.5, 6.5, 4.5, 8.5))
  expect_equal(default_nrow(ext, 256L), 256L)
})
