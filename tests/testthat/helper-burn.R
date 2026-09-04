## A hand-built <cb_burn> so that the pure-R helpers can be tested without
## running the Rust rasterizer: a 4x4-cell square on a 10x10 grid over
## [0, 10] x [0, 10], with half-covered edge cells and quarter-covered corners.
## Interior full cells are rows 3:5, cols 4:6; the boundary ring is rows 2:6,
## cols 3:7.
square_burn <- function(coverage = TRUE, crs = "EPSG:3857") {
  ring <- expand.grid(row = 2:6, col = 3:7)
  ring <- ring[!(ring$row %in% 3:5 & ring$col %in% 4:6), ]
  corner <- ring$row %in% c(2L, 6L) & ring$col %in% c(3L, 7L)
  structure(
    list(
      runs = data.frame(row = 3:5, col_start = 4L, col_end = 6L, id = 1L),
      edges = data.frame(row = ring$row, col = ring$col,
                         fraction = ifelse(corner, 0.25, 0.5), id = 1L),
      lines = data.frame(row = integer(), col = integer(),
                         length = numeric(), id = integer()),
      points = data.frame(row = integer(), col = integer(), id = integer()),
      notes = data.frame(geom_index = integer(), message = character())
    ),
    extent = c(0, 10, 0, 10),
    dimension = c(10L, 10L),
    crs = crs,
    coverage = coverage,
    class = "cb_burn"
  )
}
