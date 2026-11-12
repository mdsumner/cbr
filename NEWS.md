# cbr 0.1.0

Initial development version. R bindings to the `controlledburn` Rust crate for
sparse scanline rasterization of polygons, lines, and points.

## Crate

* Built on controlledburn 0.2.0. The crate moved to 0-based indices,
  half-open run ranges (`col_end` exclusive) and `id = k`; `cb_burn_wkb()`
  converts at the boundary so R sees 1-based `row`/`col`, inclusive `col_end`
  and `id`/`geom_index` equal to the 1-based input position. Nothing changes
  for R code, and the `cb_burn` tables are identical to those from 0.1.0.

* `cb_burn_wkb()` writes each output column straight into its R vector and
  frees the crate's table as soon as its columns exist, and borrows the input
  WKB rather than copying it. Peak memory for a large burn is about twice the
  result size instead of three times.

* Now built on controlledburn 0.3.1, from crates.io. Notable changes that
  reach R: multi-part geometries now get one record per cell (parts' coverage
  is summed, and a cell the parts together fill becomes a run) instead of
  per-part records, so tables for multipolygons can differ from earlier
  builds, with the same total coverage; and rings too short to burn, and
  one-coordinate lines, now get a note instead of being skipped silently.

* `cb_burn(parity = TRUE)` gives the C++ core's records exactly (no run
  coalescing, per-part records for multi-part geometries).

* `cb_burn(threads = )` burns on several threads via the crate's `burn_par`,
  with the same result as the serial burn. `0` uses all available threads;
  the default `1` is serial. Small burns stay serial regardless.

* `cb_burn(id_offset = )` adds a non-negative offset to every geometry id (and
  to `geom_index` in `notes`), recorded in attribute `id_offset`; `n_geom` stays
  the number of input geometries. `cb_combine()` joins burns of consecutive
  slices of a collection that share a grid (identical `extent`, `dimension`,
  `crs` and `coverage`) into the burn of the whole, record for record,
  whatever order they are given in. It refuses overlapping id ranges; a gap in
  the ids counts as empty geometries. `summary()`, `plot()` and `cb_crop()`
  respect the offset.

* `cb_burn_index()` burns onto affine grids, rotated and sheared, given a GDAL
  geotransform and `dimension`. Geometry is moved into index space by the
  inverse of the geotransform and burned on a unit-cell lattice, so edge
  fractions are exact fractions of each cell's own area, whatever its shape.
  The result is a `cb_burn` whose `extent` is the index lattice, with the world
  definition in attributes `geotransform`, `world_extent` and `crs`.
  `cb_geotransform()` builds a geotransform from an origin, cell size,
  rotation and shear. `cb_cells()` returns the cells as world-space polygons
  for plotting. `cb_crop()` moves an affine burn's origin to the window,
  `cb_combine()` requires matching geotransforms, and `cb_materialize()` adds
  `geotransform` to the `gis` list (with `bbox` the world bounding box of the
  grid's outline). `lines$length` is in world units on grids whose cells are
  similar to squares, and in cell units, with a warning, otherwise. See the
  "Rotated and sheared grids" vignette for a gallery.

* New vignette, "Burning onto a model's native ocean grid": a polygon burned
  exactly onto the rectilinear part of the ACCESS-ESM1.6 ocean (MOM tripolar)
  grid in its native array order, with edges cut at the cell walls taken from
  the model's corner coordinates, the 80°E seam handled by burning on two
  periods and folding, figures with `ximage`, and an area-weighted mean of
  temperature over a Southern Ocean sector. A 317 KB extract of the model's
  grid and one temperature field ships in `inst/extdata`. Suggests `ximage`.

* New article, "Cell predicates: exact all-touched and slivers": the DE-9IM
  relations on Natural Earth countries at one degree. cbr's exact all-touched
  (`cb_materialize(value = "intersects")`) against GDAL's all-touched via
  `terra::rasterize(touches = TRUE)`: 82 cells cbr finds that GDAL misses
  (borders on cell walls) and 4 GDAL burns that do not touch the country, all
  settled by GEOS. And the slivers the burn's tolerance rounds away, with a
  check for invalid input: one of the crate book's 18 slivers comes from
  Sudan's self-intersecting polygon in Natural Earth. The countries ship in
  `inst/extdata` (public domain). Suggests `terra`.

* The pkgdown site is built and deployed by GitHub Actions
  (`.github/workflows/pkgdown.yaml`), with a hex logo drawn by cbr itself
  (`data-raw/logo.R`).

* `cb_relate()` and `cb_burn(relate = TRUE)` (also `cb_burn_index(relate = )`)
  add a sixth table, `relations`: for each polygon, the DE-9IM intersection
  matrix (nine characters, e.g. `"212F11FF2"`) against every cell whose closed
  box its boundary touches, including cells it only touches along a wall or at
  a corner, which have no record in the burn. Every other cell is `"212FF1FF2"`
  if it is in the polygon's runs and `"FF2FF1212"` otherwise. The matrices come
  from the geometry, not the burn's records, so they see the slivers the burn's
  1e-6 coverage tolerance rounds away. Computed by the crate's `relate_cells`
  (controlledburn 0.3.1); matrices match `sf::st_relate` cell by cell in the
  tests, including on rotated and sheared grids. `cb_crop()` crops the table
  with the others and `cb_combine()` binds it (all burns or none). `cb_relate()`
  needs the geometry again, since a burn does not keep it, and refuses cropped
  burns.

* Spatial predicates per cell, from a burn's runs and `relations`:
  `cb_matrix()` looks up the DE-9IM matrix of any cell (relations where the
  boundary touches, otherwise covered for run cells and disjoint for the
  rest; `NA` off the grid and for lines and points). `cb_cells_where()`
  selects cells by predicate (`"intersects"`, `"interiors_meet"`, `"covers"`,
  `"touches"`, `"covered_by"`) or by DE-9IM pattern. `cb_materialize()` gains
  those predicates as values, counting the polygons satisfying each per cell;
  `value = "intersects"` is the exact all-touched mask, including wall and
  corner contacts and the slivers the burn rounds away. `cb_slivers()` lists
  the cells where the burn's records and the exact topology disagree. Tested
  against GEOS (`sf::st_relate` and its named predicates) on every cell; on
  Natural Earth countries at one degree it reproduces the crate's figures
  (12 interior and 6 exterior slivers, 78 touching cells).

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

* `plot(what = "exact")` on small grids (no side longer than `detail`, default
  64 cells) also shows the record structure: a thin white line through the
  cells of each run, and each edge cell numbered by its position in the
  `edges` table, so the records can be counted off the picture.

* `cb_crop()` with no `extent` crops to the content extent.

* `cb_burn()` records the number of input geometries in attribute `n_geom`.

* `print()` for `cb_burn` summarizes the grid, mode, table sizes, and (for a
  cropped burn) the offset within the parent grid.
