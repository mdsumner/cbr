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

* `print()` for `cb_burn` summarizes the grid, mode, table sizes, and (for a
  cropped burn) the offset within the parent grid.
