# Plot a burn

Draw a
[`cb_burn()`](https://mdsumner.github.io/cbr/reference/cb_burn.md) in
map coordinates (`asp = 1`) with the grid extent as the frame. Two
renderings are available:

## Usage

``` r
# S3 method for class 'cb_burn'
plot(
  x,
  what = c("auto", "exact", "overview"),
  marginals = TRUE,
  col = grDevices::hcl.colors(64, "YlOrRd", rev = TRUE),
  by_id = FALSE,
  tiles = 256L,
  threshold = 50000,
  content = FALSE,
  main = NULL,
  detail = 64L,
  ...
)
```

## Arguments

- x:

  A `cb_burn` object.

- what:

  `"auto"`, `"exact"`, or `"overview"`.

- marginals:

  Draw the row and column marginal strips. The layout is sized from the
  extent's aspect ratio so the strips span exactly the map; leftover
  device space is left blank.

- col:

  Colour ramp (a vector of colours, low to high coverage) for
  `"overview"` and for edge shading in `"exact"`.

- by_id:

  In `"exact"` mode, colour runs and edges by geometry `id` instead of
  by coverage.

- tiles:

  Passed to
  [`cb_overview()`](https://mdsumner.github.io/cbr/reference/cb_overview.md)
  for the `"overview"` rendering.

- threshold:

  Record count above which `"auto"` switches to `"overview"`.

- content:

  Draw the content extent (bounding box of touched cells) as a dashed
  rectangle.

- main:

  Optional title, drawn in the outer margin above the strips.

- detail:

  Largest grid side (cells) at which `"exact"` mode draws the run lines
  and edge numbers. `0` never, `Inf` always.

- ...:

  Passed to
  [`graphics::plot()`](https://rdrr.io/r/graphics/plot.default.html) for
  the main panel.

## Value

`x`, invisibly. The summary computed for the marginals is attached as
attribute `summary`.

## Details

- `"exact"`:

  every record is drawn: runs as rectangles, edge cells as rectangles
  shaded by coverage fraction, line cells shaded by length, points as
  cells. Cost is one rectangle per record, so this is the picture of the
  sparse representation itself.

- `"overview"`:

  the coarse mean-coverage raster from
  [`cb_overview()`](https://mdsumner.github.io/cbr/reference/cb_overview.md)
  drawn with
  [`graphics::image()`](https://rdrr.io/r/graphics/image.html). Cost
  does not depend on grid size.

`"auto"` (default) picks `"exact"` when the total record count is at
most `threshold`, else `"overview"`. On small grids (no side longer than
`detail` cells) the exact rendering also shows the record structure: a
thin white line through the cells of each run, so runs can be seen and
counted, and each edge cell numbered by its position in the `edges`
table. Either way the optional marginal strips are exact: covered cells
per grid row (right) and per grid column (top), from
[`summary.cb_burn()`](https://mdsumner.github.io/cbr/reference/summary.cb_burn.md).
