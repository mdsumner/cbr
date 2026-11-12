px2w <- function(gt, col, row) {
  c(gt[1] + col * gt[2] + row * gt[3], gt[4] + col * gt[5] + row * gt[6])
}
shoelace <- function(m) {
  n <- nrow(m)
  abs(sum(m[-n, 1] * m[-1, 2] - m[-1, 1] * m[-n, 2])) / 2
}
covered_cells <- function(b) {
  sum(b$runs$col_end - b$runs$col_start + 1L) + sum(b$edges$fraction)
}
cell_area <- function(gt) abs(gt[2] * gt[6] - gt[3] * gt[5])
ring_wkt <- function(m) {
  sprintf("POLYGON ((%s))", paste(sprintf("%.15g %.15g", m[, 1], m[, 2]), collapse = ", "))
}
## a star about world point `c0` with alternating radii
star <- function(c0, r1, r2, n = 5) {
  a <- seq(0, 2 * pi, length.out = 2 * n + 1)
  r <- rep(c(r1, r2), length.out = 2 * n + 1)
  r[2 * n + 1] <- r[1]
  cbind(c0[1] + r * cos(a), c0[2] + r * sin(a))
}

test_that("cb_geotransform builds GDAL-order geotransforms", {
  expect_equal(cb_geotransform(), c(0, 1, 0, 0, 0, -1))
  expect_equal(cb_geotransform(c(5, 9), res = c(2, 3)), c(5, 2, 0, 9, 0, -3))
  ## columns run up the page, rows run to the right
  expect_equal(cb_geotransform(c(10, 20), res = 2, rotation = 90),
               c(10, 0, 2, 20, 2, 0))
  ## shear turns the row direction away from perpendicular
  gt <- cb_geotransform(res = 1, shear = 45)
  expect_equal(gt[c(3, 6)], c(sin(pi / 4), -cos(pi / 4)))
  expect_error(cb_geotransform(res = -1), "res")
  expect_error(cb_geotransform(origin = 1), "origin")
  expect_error(cb_geotransform(rotation = NA), "rotation")
})

test_that("a north-up geotransform gives the regular burn", {
  gt <- c(1e6, 10, 0, 6e6, 0, -10)
  g <- wk::wkt("POLYGON ((1000013.3 5999871.7, 1000241.9 5999905.2,
                          1000180.4 5999790.6, 1000013.3 5999871.7))")
  a <- cb_burn(g, extent = c(1e6, 1e6 + 300, 6e6 - 200, 6e6), dimension = c(30, 20))
  b <- cb_burn_index(g, gt, c(30, 20))
  for (nm in names(a)) {
    expect_equal(b[[nm]], a[[nm]], tolerance = 1e-5)
  }
})

test_that("affine burns carry the world definition", {
  gt <- cb_geotransform(c(100, 200), res = c(2, 1.5), rotation = 30, shear = 10)
  b <- cb_burn_index(wk::wkt("POINT (101 195)"), gt, c(8, 6), crs = "EPSG:32755")
  expect_s3_class(b, "cb_burn")
  expect_equal(attr(b, "geotransform"), gt)
  expect_equal(attr(b, "extent"), c(0, 8, 0, 6))
  expect_equal(attr(b, "dimension"), c(8L, 6L))
  expect_equal(attr(b, "crs"), "EPSG:32755")
  corners <- rbind(px2w(gt, 0, 0), px2w(gt, 8, 0), px2w(gt, 0, 6), px2w(gt, 8, 6))
  expect_equal(attr(b, "world_extent"),
               c(range(corners[, 1]), range(corners[, 2])))
  expect_output(print(b), "affine")
})

test_that("a point lands in the cell that contains it", {
  gt <- cb_geotransform(c(100, 200), res = c(2, 1.5), rotation = 30, shear = 10)
  for (rc in list(c(1, 1), c(4, 3), c(8, 6), c(2, 5))) {
    p <- px2w(gt, rc[1] - 0.5, rc[2] - 0.5)
    b <- cb_burn_index(wk::wkt(sprintf("POINT (%.15g %.15g)", p[1], p[2])), gt, c(8, 6))
    expect_equal(c(b$points$col, b$points$row), rc)
  }
  p <- px2w(gt, 9.5, 2.5)
  b <- cb_burn_index(wk::wkt(sprintf("POINT (%.15g %.15g)", p[1], p[2])), gt, c(8, 6))
  expect_equal(nrow(b$points), 0L)
})

test_that("rotated and sheared grids conserve area", {
  for (gt in list(
    cb_geotransform(c(2, 13), res = 1, rotation = 25),
    cb_geotransform(c(0, 9), res = 1, rotation = 0, shear = 20),
    cb_geotransform(c(5, 20), res = c(2, 0.75), rotation = -15, shear = -10)
  )) {
    dm <- c(14L, 10L)
    centre <- px2w(gt, 7, 5)
    ring <- star(centre, 1.8 * min(abs(gt[c(2, 6)])), 0.8 * min(abs(gt[c(2, 6)])))
    b <- cb_burn_index(wk::wkt(ring_wkt(ring)), gt, dm)
    expect_equal(covered_cells(b) * cell_area(gt), shoelace(ring), tolerance = 1e-5)
  }
})

test_that("fractions match an independent clip of each cell", {
  skip_if_not_installed("sf")
  gt <- cb_geotransform(c(100, 200), res = c(2, 1.5), rotation = 30, shear = 10)
  dm <- c(8L, 6L)
  centre <- px2w(gt, 4, 3)
  star_ring <- star(centre, 8, 3.5, n = 6)
  sq <- function(c0, h) {
    rbind(c0 + c(-h, -h), c0 + c(h, -h), c0 + c(h, h), c0 + c(-h, h), c0 + c(-h, -h))
  }
  holey <- sprintf("POLYGON ((%s), (%s))",
                   paste(sprintf("%.15g %.15g", sq(centre, 6)[, 1], sq(centre, 6)[, 2]), collapse = ", "),
                   paste(sprintf("%.15g %.15g", sq(centre, 2)[, 1], sq(centre, 2)[, 2]), collapse = ", "))
  geoms <- list(ring_wkt(star_ring), holey)

  cells <- lapply(seq_len(dm[2]) - 1L, function(r) {
    lapply(seq_len(dm[1]) - 1L, function(c) {
      m <- rbind(px2w(gt, c, r), px2w(gt, c + 1, r), px2w(gt, c + 1, r + 1),
                 px2w(gt, c, r + 1), px2w(gt, c, r))
      sf::st_polygon(list(m))
    })
  })
  cells <- unlist(cells, recursive = FALSE)

  for (w in geoms) {
    target <- sf::st_as_sfc(w)[[1]]
    ref <- vapply(cells, function(p) {
      i <- sf::st_intersection(p, target)
      if (length(i)) sf::st_area(i) / sf::st_area(p) else 0
    }, numeric(1))
    b <- cb_burn_index(wk::wkt(w), gt, dm)
    got <- as.vector(cb_materialize(b, background = 0))
    expect_equal(got, ref, tolerance = 1e-5)
    ## a case worth having: it touches partial cells and full ones
    expect_true(any(ref > 0 & ref < 1) && any(ref == 1))
  }
})

test_that("the cell-centre rule works on a rotated grid", {
  gt <- cb_geotransform(c(0, 12), res = 1, rotation = 20)
  b <- cb_burn_index(wk::wkt(ring_wkt(star(px2w(gt, 7, 5), 3.5, 3.5, n = 3))),
                     gt, c(14, 10), coverage = FALSE)
  expect_equal(nrow(b$edges), 0L)
  expect_gt(nrow(b$runs), 0L)
})

test_that("line lengths are in world units on similar cells", {
  gt <- cb_geotransform(c(0, 20), res = 2, rotation = 20)
  p <- px2w(gt, 1.3, 2.2); q <- px2w(gt, 9.1, 5.7)
  b <- cb_burn_index(wk::wkt(sprintf("LINESTRING (%.15g %.15g, %.15g %.15g)", p[1], p[2], q[1], q[2])),
                     gt, c(12, 8))
  expect_equal(sum(b$lines$length), sqrt(sum((p - q)^2)), tolerance = 1e-5)

  skewed <- cb_geotransform(c(0, 20), res = 2, rotation = 20, shear = 15)
  expect_warning(
    cb_burn_index(wk::wkt(sprintf("LINESTRING (%.15g %.15g, %.15g %.15g)", p[1], p[2], q[1], q[2])),
                  skewed, c(12, 8)),
    "cell units"
  )
})

test_that("cb_materialize records where an affine grid is", {
  gt <- cb_geotransform(c(100, 200), res = c(2, 1.5), rotation = 30)
  b <- cb_burn_index(wk::wkt(ring_wkt(star(px2w(gt, 4, 3), 5, 2.5))), gt, c(8, 6))
  v <- cb_materialize(b)
  gis <- attr(v, "gis")
  we <- attr(b, "world_extent")
  expect_equal(gis$bbox, c(we[1], we[3], we[2], we[4]))
  expect_equal(gis$dim, c(8, 6, 1))
  expect_equal(gis$geotransform, gt)
  ## regular burns are unchanged
  r <- cb_materialize(cb_burn(wk::wkt("POINT (1.5 1.5)"), extent = c(0, 3, 0, 3), dimension = c(3, 3)))
  expect_null(attr(r, "gis")$geotransform)
})

test_that("cropping an affine burn moves its origin, not its cells", {
  gt <- cb_geotransform(c(100, 200), res = c(2, 1.5), rotation = 30, shear = 10)
  b <- cb_burn_index(wk::wkt(ring_wkt(star(px2w(gt, 5, 4), 4, 2, n = 4))), gt, c(12, 9))
  cr <- cb_crop(b)
  off <- attr(cr, "parent")$offset
  expect_true(all(off >= 0) && any(off > 0))
  expect_equal(attr(cr, "extent"), c(0, attr(cr, "dimension")[1], 0, attr(cr, "dimension")[2]))
  expect_equal(c(attr(cr, "geotransform")[1], attr(cr, "geotransform")[4]),
               px2w(gt, off[1], off[2]))
  expect_equal(attr(cr, "geotransform")[c(2, 3, 5, 6)], gt[c(2, 3, 5, 6)])

  centroids <- function(x) {
    cc <- cb_cells(x)
    m <- rowsum(cbind(x = cc$x, y = cc$y), cc$cell) / 4
    d <- cbind(m, value = cc$value[!duplicated(cc$cell)])
    d[order(round(d[, "x"], 6), round(d[, "y"], 6)), ]
  }
  expect_equal(unname(centroids(cr)), unname(centroids(b)), tolerance = 1e-8)
  expect_equal(sum(cb_materialize(cr, background = 0)),
               sum(cb_materialize(b, background = 0)), tolerance = 1e-6)
})

test_that("cb_combine needs the same affine grid", {
  gt <- cb_geotransform(c(0, 10), res = 1, rotation = 10)
  gt2 <- cb_geotransform(c(0, 10), res = 1, rotation = 11)
  g <- wk::wkt(c("POLYGON ((3 5, 6 5, 6 8, 3 8, 3 5))", "POINT (4 6)"))
  a <- cb_burn_index(g[1], gt, c(10, 10))
  expect_error(cb_combine(a, cb_burn_index(g[2], gt2, c(10, 10), id_offset = 1)),
               "geotransform")
  whole <- cb_burn_index(g, gt, c(10, 10))
  comb <- cb_combine(a, cb_burn_index(g[2], gt, c(10, 10), id_offset = 1))
  expect_equal(comb, whole)
})

test_that("cb_cells returns world-space cell polygons", {
  gt <- cb_geotransform(c(100, 200), res = c(2, 1.5), rotation = 30, shear = 10)
  b <- cb_burn_index(wk::wkt(ring_wkt(star(px2w(gt, 4, 3), 5, 2.5))), gt, c(8, 6))
  cc <- cb_cells(b)
  ncell <- sum(b$runs$col_end - b$runs$col_start + 1L) + nrow(b$edges)
  expect_equal(nrow(cc), 4L * ncell)
  expect_named(cc, c("cell", "x", "y", "row", "col", "id", "kind", "value"))
  one <- cc[cc$cell == 1L, ]
  corners <- rbind(px2w(gt, one$col[1] - 1, one$row[1] - 1), px2w(gt, one$col[1], one$row[1] - 1),
                   px2w(gt, one$col[1], one$row[1]), px2w(gt, one$col[1] - 1, one$row[1]))
  expect_equal(unname(as.matrix(one[c("x", "y")])), unname(corners))
  ## cell areas are the geotransform's
  areas <- vapply(split(cc, cc$cell), function(p) shoelace(rbind(as.matrix(p[c("x", "y")]),
                                                               as.matrix(p[1, c("x", "y")]))),
                  numeric(1))
  expect_equal(unname(areas), rep(cell_area(gt), ncell), tolerance = 1e-8)
  expect_error(cb_cells(b, max_cells = 3), "max_cells")
  empty <- cb_burn_index(wk::wkt("POINT (1000 1000)"), gt, c(8, 6))
  expect_equal(nrow(cb_cells(empty)), 0L)
  regular <- cb_burn(wk::wkt("POINT (1.5 2.5)"), extent = c(0, 4, 0, 4), dimension = c(4, 4))
  expect_equal(unname(as.matrix(cb_cells(regular)[c("x", "y")])),
               cbind(c(1, 2, 2, 1), c(3, 3, 2, 2)))
})

test_that("bad geotransforms and dimensions are refused", {
  g <- wk::wkt("POINT (1 1)")
  expect_error(cb_burn_index(g, 1:5, c(4, 4)), "six")
  expect_error(cb_burn_index(g, c(0, 1, 0, 0, NA, -1), c(4, 4)), "six")
  expect_error(cb_burn_index(g, c(0, 1, 2, 0, 2, 4), c(4, 4)), "singular")
  expect_error(cb_burn_index(g, c(0, 0, 0, 0, 0, 0), c(4, 4)), "singular")
  expect_error(cb_burn_index(g, c(0, 1, 0, 0, 0, -1)), "dimension")
  expect_error(cb_burn_index(g, c(0, 1, 0, 0, 0, -1), c(4, 0)), "dimension")
})
