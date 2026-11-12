two_halves <- wk::wkt(
  "MULTIPOLYGON (((0.5 0.5, 1.5 0.5, 1.5 2.5, 0.5 2.5, 0.5 0.5)),
                 ((1.5 0.5, 2.5 0.5, 2.5 2.5, 1.5 2.5, 1.5 0.5)))"
)

covered <- function(b) sum(cb_materialize(b, background = 0))

test_that("multi-part geometries get one record per cell by default", {
  merged <- cb_burn(two_halves, extent = c(0, 3, 0, 3), dimension = c(3, 3))
  parity <- cb_burn(two_halves, extent = c(0, 3, 0, 3), dimension = c(3, 3),
                    parity = TRUE)

  expect_lt(nrow(merged$edges), nrow(parity$edges))
  expect_false(anyDuplicated(merged$edges[c("row", "col", "id")]) > 0)
  expect_equal(covered(merged), covered(parity), tolerance = 1e-6)
  expect_equal(covered(merged), 4, tolerance = 1e-6)
})

test_that("threads give the same burn as the serial path", {
  n <- 200
  set.seed(4217)
  x <- runif(n, 0, 90)
  y <- runif(n, 0, 90)
  sq <- wk::wkt(sprintf(
    "POLYGON ((%1$f %2$f, %3$f %2$f, %3$f %4$f, %1$f %4$f, %1$f %2$f))",
    x, y, x + 5.3, y + 4.7
  ))
  args <- list(sq, extent = c(0, 100, 0, 100), dimension = c(500, 500))
  serial <- do.call(cb_burn, args)
  expect_equal(do.call(cb_burn, c(args, threads = 4)), serial, ignore_attr = "threads")
  expect_equal(do.call(cb_burn, c(args, threads = 0)), serial, ignore_attr = "threads")
})

test_that("unparseable input is noted with threads, as serially", {
  good <- unclass(wk::as_wkb(wk::wkt("POLYGON ((1 1, 4 1, 4 4, 1 4, 1 1))")))[[1L]]
  blobs <- structure(list(good, as.raw(c(1, 2, 3)), raw(), good), class = "wk_wkb")
  serial <- cb_burn(blobs, extent = c(0, 5, 0, 5), dimension = c(5, 5))
  par <- cb_burn(blobs, extent = c(0, 5, 0, 5), dimension = c(5, 5), threads = 2)
  expect_equal(serial$notes$geom_index, 2L)
  expect_equal(par$notes, serial$notes)
  expect_equal(par$runs, serial$runs)
})

test_that("`threads` is validated", {
  expect_error(cb_burn(two_halves, threads = -1), "threads")
  expect_error(cb_burn(two_halves, threads = NA), "threads")
})
