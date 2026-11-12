# Package index

## Burning

Rasterize geometry onto a grid. The result is sparse: runs of fully
covered cells, boundary cells with their exact coverage fraction, and
tables for lines, points and notes.

- [`cb_burn()`](https://mdsumner.github.io/cbr/reference/cb_burn.md) :
  Burn geometries onto a regular grid
- [`cb_combine()`](https://mdsumner.github.io/cbr/reference/cb_combine.md)
  : Combine burns of consecutive slices of a collection

## Rotated and sheared grids

Burn onto affine grids given a GDAL geotransform, in index space.

- [`cb_burn_index()`](https://mdsumner.github.io/cbr/reference/cb_burn_index.md)
  : Burn geometries onto an affine (rotated or sheared) grid
- [`cb_geotransform()`](https://mdsumner.github.io/cbr/reference/cb_geotransform.md)
  : Build a GDAL-style geotransform
- [`cb_cells()`](https://mdsumner.github.io/cbr/reference/cb_cells.md) :
  Cells of a burn as polygons in world coordinates

## Cropping and materializing

Cut a burn to a window of its grid, and expand it to a dense raster in
the layout of `gdalraster::read_ds()`.

- [`cb_crop()`](https://mdsumner.github.io/cbr/reference/cb_crop.md) :
  Crop a burn to a sub-window of its grid
- [`cb_materialize()`](https://mdsumner.github.io/cbr/reference/cb_materialize.md)
  : Materialize a burn as a dense raster
- [`cb_as_matrix()`](https://mdsumner.github.io/cbr/reference/cb_as_matrix.md)
  : Convert a materialized raster to a matrix

## Summaries and plots

- [`summary(`*`<cb_burn>`*`)`](https://mdsumner.github.io/cbr/reference/summary.cb_burn.md)
  : Summarize a burn without materializing it
- [`cb_overview()`](https://mdsumner.github.io/cbr/reference/cb_overview.md)
  : Coarse overview raster of a burn
- [`plot(`*`<cb_burn>`*`)`](https://mdsumner.github.io/cbr/reference/plot.cb_burn.md)
  : Plot a burn

## Cell predicates (DE-9IM)

The DE-9IM relation of each polygon to the grid cells, exact where the
boundary touches a cell, and the predicates and masks that follow.

- [`cb_relate()`](https://mdsumner.github.io/cbr/reference/cb_relate.md)
  : DE-9IM relations of polygons to the cells their boundaries touch
- [`cb_matrix()`](https://mdsumner.github.io/cbr/reference/cb_matrix.md)
  : DE-9IM matrix of any cell
- [`cb_cells_where()`](https://mdsumner.github.io/cbr/reference/cb_cells_where.md)
  : Select the cells that satisfy a spatial predicate
- [`cb_slivers()`](https://mdsumner.github.io/cbr/reference/cb_slivers.md)
  : Cells where the burn's records and the exact topology disagree

## Low-level bindings

The functions that call into Rust, for building other front ends.

- [`cb_burn_wkb()`](https://mdsumner.github.io/cbr/reference/cb_burn_wkb.md)
  : Burn WKB geometries onto a regular grid.
- [`cb_relate_wkb()`](https://mdsumner.github.io/cbr/reference/cb_relate_wkb.md)
  : DE-9IM of polygons against grid cells, for an existing burn.
