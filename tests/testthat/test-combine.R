chunk_geoms <- function() {
  wk::wkt(c(
    "POLYGON ((1 1, 4.5 1, 4.5 4.5, 1 4.5, 1 1))",
    "MULTIPOLYGON (((5.5 5.5, 7.5 5.5, 7.5 8.5, 5.5 8.5, 5.5 5.5)),
                   ((7.5 5.5, 9.5 5.5, 9.5 8.5, 7.5 8.5, 7.5 5.5)))",
    "LINESTRING (0.5 9.5, 9.5 0.7)",
    "POINT (2.5 8.5)",
    "POLYGON ((3 3, 6 3, 6 6, 3 6, 3 3))"
  ))
}
grid10 <- list(extent = c(0, 10, 0, 10), dimension = c(10L, 10L))
burn10 <- function(g, ...) {
  cb_burn(g, extent = grid10$extent, dimension = grid10$dimension, ...)
}

test_that("id_offset shifts ids and note indices", {
  g <- wk::as_wkb(chunk_geoms())
  g <- structure(c(unclass(g), list(as.raw(c(1, 2, 3)))), class = class(g))
  b0 <- burn10(g)
  b7 <- burn10(g, id_offset = 7)
  expect_equal(attr(b7, "id_offset"), 7L)
  expect_equal(attr(b7, "n_geom"), 6L)
  expect_equal(b7$runs$id, b0$runs$id + 7L)
  expect_equal(b7$edges$id, b0$edges$id + 7L)
  expect_equal(b7$notes$geom_index, b0$notes$geom_index + 7L)
  expect_equal(b0$notes$geom_index, 6L)
})

test_that("`id_offset` is validated", {
  g <- chunk_geoms()
  expect_error(burn10(g, id_offset = -1), "id_offset")
  expect_error(burn10(g, id_offset = 1.5), "id_offset")
  expect_error(burn10(g, id_offset = NA), "id_offset")
  expect_error(burn10(g, id_offset = .Machine$integer.max), "id_offset")
})

test_that("combining chunks equals the burn of the whole collection", {
  g <- chunk_geoms()
  for (args in list(list(), list(coverage = FALSE), list(parity = TRUE))) {
    whole <- do.call(burn10, c(list(g), args))
    a <- do.call(burn10, c(list(g[1:2]), args))
    b <- do.call(burn10, c(list(g[3L]), args, id_offset = 2))
    d <- do.call(burn10, c(list(g[4:5]), args, id_offset = 3))
    expect_equal(cb_combine(a, b, d), whole)
    expect_equal(cb_combine(list(d, a, b)), whole)
  }
})

test_that("combine reproduces the summary and materialized ids", {
  g <- chunk_geoms()
  whole <- burn10(g)
  comb <- cb_combine(burn10(g[1:2]), burn10(g[3:5], id_offset = 2))
  expect_equal(cb_materialize(comb, "id"), cb_materialize(whole, "id"))
  expect_equal(summary(comb)$area, summary(whole)$area)
})

test_that("summary and plot respect an offset", {
  g <- chunk_geoms()
  b0 <- burn10(g)
  b5 <- burn10(g, id_offset = 5)
  s0 <- summary(b0)
  s5 <- summary(b5)
  expect_equal(s5$area$id, s0$area$id + 5L)
  expect_equal(s5$area[-1L], s0$area[-1L])
  expect_equal(s5$n_geom, 5L)
  grDevices::png(tempfile(fileext = ".png"))
  on.exit(grDevices::dev.off())
  expect_invisible(plot(b5, by_id = TRUE))
})

test_that("crop keeps the offset", {
  b <- burn10(chunk_geoms(), id_offset = 3)
  expect_equal(attr(cb_crop(b, c(0, 5, 0, 5)), "id_offset"), 3L)
})

test_that("combine refuses mismatched grids and overlapping ids", {
  g <- chunk_geoms()
  a <- burn10(g[1:2])
  expect_error(cb_combine(a, burn10(g[3], id_offset = 2, coverage = FALSE)), "coverage")
  expect_error(cb_combine(a, cb_burn(g[3], extent = c(0, 20, 0, 10),
                                     dimension = c(10, 10), id_offset = 2)), "extent")
  expect_error(cb_combine(a, burn10(g[3], id_offset = 1)), "overlap")
  expect_error(cb_combine(a, 1), "cb_burn")
  expect_error(cb_combine(), "nothing")
})

test_that("a gap in the ids counts as empty geometries", {
  g <- chunk_geoms()
  comb <- cb_combine(burn10(g[1]), burn10(g[2], id_offset = 3))
  expect_equal(attr(comb, "n_geom"), 4L)
  expect_equal(summary(comb)$area$id, 1:4)
  expect_equal(summary(comb)$area$interior[2:3], c(0, 0))
})
