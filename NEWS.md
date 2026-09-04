# cbr 0.0.0.9000

Initial development version. R bindings to the `controlledburn` Rust crate for
sparse scanline rasterization of polygons, lines, and points.

## Burning

* `cb_burn()` rasterizes any geometry that `wk::as_wkb()` understands onto a
  regular grid given `extent = c(xmin, xmax, ymin, ymax)` and
  `dimension = c(ncol, nrow)`. The result is a `cb_burn` object: five data
  frames (`runs`, `edges`, `lines`, `points`, `notes`) with the grid stored in
  attributes. Indices are 1-based and row 1 is the top row; geometry `k` gets
  `id = k`.

* `extent` defaults to the bounding box of `geom` (`wk::wk_bbox()`, reordered
  to `c(xmin, xmax, ymin, ymax)`), and `dimension` defaults to 256 columns with
  the number of rows computed so that cells are square in the units of
  `extent`. A single `dimension` value is taken as `ncol`.

* `coverage = TRUE` (default) returns polygon interiors as horizontal runs of
  fully covered cells and boundary cells with their exact coverage fraction in
  `edges`. `coverage = FALSE` uses the faster cell-centre (fasterize-style)
  rule and emits runs only.

* `cb_burn_wkb()` is the low-level worker taking a list of raw WKB blobs and
  the grid parameters directly.

* Non-fatal per-geometry problems (for example unclosed rings, which are closed
  silently) are reported in the `notes` table rather than raised as errors.

## Cropping and materializing

* `cb_crop()` snaps an extent to the cell boundaries of a burn grid
  (`snap = "out"`, `"in"`, or `"near"`) and subsets the tables to that window.
  No geometry is re-rasterized: runs that straddle the window edge are clipped
  and all indices are re-based, so the result is a valid `cb_burn` on the
  smaller grid. The parent grid and the window's `c(col, row)` offset within it
  are recorded in the `parent` attribute.

* `cb_materialize()` expands a burn into a dense cell vector in the layout of
  `gdalraster::read_ds()`: one element per cell, row-major from the top-left,
  with a `gis` attribute (`type`, `bbox`, `dim`, `srs`, `datatype`). This
  format is produced without any dependency on gdalraster. `value` selects the
  cell contents: `"coverage"` (summed coverage, Float32), `"id"` (highest
  geometry id wins, Int32), or `"count"` (distinct geometries per cell,
  Int32). Untouched cells take `background` (default `NA`).

* `cb_as_matrix()` reshapes a materialized vector (from `cb_materialize()` or
  `gdalraster::read_ds()`) into an `nrow x ncol` matrix oriented as on a map.

## Summarizing and plotting

* `summary()` for `cb_burn` describes a burn in time proportional to the
  number of records, never to the number of cells: exact per-row and
  per-column covered-cell marginals (the column marginal via a difference
  array over run endpoints), the content extent of touched cells, per-geometry
  area (interior cells, summed boundary fractions, line length, point count),
  edge-fraction quantiles, and a coarse `overview` raster. `print()` reports
  these, with `ascii = TRUE` for a text rendering of the overview.

* `cb_overview()` aggregates a burn onto square tiles of `block` cells and
  returns mean coverage per tile in the `read_ds` layout, so it is the same
  shape as `cb_materialize()` and equal to it when `block = 1`. Runs spanning
  several tiles are handled with a difference array along the tile row, so
  the cost is records plus tiles regardless of block size.

* `plot()` for `cb_burn` draws the grid extent with `asp = 1` and either every
  record as a rectangle (`what = "exact"`: runs, edge cells shaded by
  fraction, optionally coloured `by_id`) or the overview raster
  (`what = "overview"`), switching automatically on a record-count
  `threshold`. Marginal strips along the top and right show the exact column
  and row totals even when the main panel is an overview.

* `cb_crop()` with no `extent` crops to the content extent.

* `cb_burn()` records the number of input geometries in attribute `n_geom`.

* `print()` for `cb_burn` summarizes the grid, mode, table sizes, and (for a
  cropped burn) the offset within the parent grid.
