test_that("summary is consistent with materialize on a small burn", {
  b <- square_burn()
  s <- summary(b)
  expect_s3_class(s, "summary.cb_burn")
  m <- cb_as_matrix(cb_materialize(b, background = 0))
  ## block 1 overview is the exact coverage
  ov <- s$overview
  expect_equal(attr(ov, "block"), 1L)
  expect_equal(as.numeric(ov), as.numeric(cb_materialize(b)))
  ## marginals are exact row / col sums
  expect_equal(s$row_marginal, rowSums(m))
  expect_equal(s$col_marginal, colSums(m))
  ## content is the touched ring
  expect_equal(s$content$rows, c(2L, 6L))
  expect_equal(s$content$cols, c(3L, 7L))
  expect_equal(s$content$extent, c(2, 7, 4, 9))
  ## per-geometry area
  expect_equal(s$area$interior, 9)
  expect_equal(s$area$boundary, 7)
  expect_equal(s$n_geom, 1L)
  expect_output(print(s), "content:  cols 3..7, rows 2..6")
  expect_output(print(s, ascii = TRUE), "@")
})

test_that("overview tiles aggregate coverage, including run-spanning tiles", {
  b <- square_burn()
  m <- cb_as_matrix(cb_materialize(b, background = 0))
  for (blk in c(2L, 3L, 4L, 7L)) {
    ov <- cb_overview(b, block = blk)
    ri <- (seq_len(10) - 1L) %/% blk
    agg <- t(rowsum(t(rowsum(m, ri)), ri)) / outer(tabulate(ri + 1L), tabulate(ri + 1L))
    agg[agg == 0] <- NA
    expect_equal(as.numeric(ov), as.numeric(t(agg)), info = paste("block", blk))
    expect_equal(attr(ov, "gis")$dim[1:2], rep(ceiling(10 / blk), 2))
  }
  ## a long run spanning several tiles: interior tiles are full
  b$runs <- data.frame(row = 1L, col_start = 2L, col_end = 9L, id = 1L)
  b$edges <- b$edges[0, ]
  ov <- cb_overview(b, block = 2L)
  expect_equal(cb_as_matrix(ov)[1, ], c(0.25, 0.5, 0.5, 0.5, 0.25))
  expect_equal(summary(b)$col_marginal, c(0, rep(1, 8), 0))
})

test_that("summary of an empty burn", {
  b <- square_burn()
  b$runs <- b$runs[0, ]
  b$edges <- b$edges[0, ]
  s <- summary(b)
  expect_null(s$content)
  expect_true(all(is.na(s$overview)))
  expect_equal(sum(s$row_marginal), 0)
  expect_output(print(s), "empty")
  expect_error(cb_crop(b), "empty")
})

test_that("cb_crop defaults to the content extent", {
  b <- square_burn()
  cr <- cb_crop(b)
  expect_equal(attr(cr, "extent"), c(2, 7, 4, 9))
  expect_equal(attr(cr, "dimension"), c(5L, 5L))
  expect_equal(sum(cb_materialize(cr), na.rm = TRUE), 16)
})

test_that("plot runs in both modes", {
  b <- square_burn()
  b$lines <- data.frame(row = 9L, col = 1:3, length = c(0.5, 1, 1.2), id = 2L)
  b$points <- data.frame(row = 10L, col = 10L, id = 3L)
  attr(b, "n_geom") <- 3L
  tf <- tempfile(fileext = ".png")
  grDevices::png(tf)
  expect_invisible(plot(b, main = "exact"))
  expect_invisible(plot(b, what = "overview", marginals = FALSE))
  expect_invisible(plot(b, by_id = TRUE, content = TRUE))
  grDevices::dev.off()
  expect_true(file.exists(tf))
  out <- plot(b)
  expect_s3_class(attr(out, "summary"), "summary.cb_burn")
})
