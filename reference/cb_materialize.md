# Materialize a burn as a dense raster

Expand the sparse tables of a
[`cb_burn()`](https://mdsumner.github.io/cbr/reference/cb_burn.md) into
a dense cell vector in the layout returned by `gdalraster::read_ds()`: a
plain vector, one element per cell, in row-major order starting at the
top-left cell, carrying a `"gis"` attribute that describes the grid
(`type`, `bbox`, `dim`, `srs`, `datatype`). This is a stable format
definition, so no dependency on 'gdalraster' is needed to produce or
consume it, and the result can be handed straight to
`gdalraster::plot_raster()` or wrapped as a matrix with
[`cb_as_matrix()`](https://mdsumner.github.io/cbr/reference/cb_as_matrix.md).

## Usage

``` r
cb_materialize(
  x,
  value = c("coverage", "id", "count", "intersects", "interiors_meet", "covers",
    "touches", "covered_by"),
  background = NA
)
```

## Arguments

- x:

  A `cb_burn` object.

- value:

  What each cell should hold:

  `"coverage"`

  :   (default) summed coverage. Interior run cells contribute 1,
      boundary cells their `fraction`, line cells their `length` (CRS
      units), points 1 each. Overlapping geometries add, so values above
      1 are possible; use `pmin(., 1)` if a clamped mask is wanted.
      Output `datatype` is `"Float32"`.

  `"id"`

  :   the 1-based input position of the geometry covering the cell;
      where geometries overlap the highest `id` wins (later inputs are
      burned over earlier ones). Boundary, line, and point cells count
      as covered. Output `datatype` is `"Int32"`.

  `"count"`

  :   the number of distinct geometries with a record in the cell.
      Output `datatype` is `"Int32"`.

  `"intersects"`, `"interiors_meet"`, `"covers"`, `"touches"`, `"covered_by"`

  :   the number of polygons satisfying that predicate against the cell,
      as in
      [`cb_cells_where()`](https://mdsumner.github.io/cbr/reference/cb_cells_where.md).
      Needs a `relations` table (`cb_burn(relate = TRUE)` or
      [`cb_relate()`](https://mdsumner.github.io/cbr/reference/cb_relate.md)).
      `"intersects"` is the exact all-touched mask: unlike `"count"` it
      includes cells a polygon only touches along a wall or at a corner,
      and cells holding slivers the burn rounds away. Output `datatype`
      is `"Int32"`.

- background:

  Value for cells not touched by any geometry. Default `NA`.

## Value

A numeric (`"coverage"`) or integer (every other `value`) vector of
length `ncol * nrow` with attribute `gis`, a list: `type = "raster"`,
`bbox = c(xmin, ymin, xmax, ymax)`, `dim = c(ncol, nrow, 1)`, `srs` (the
burn's `crs`, or `""`), and `datatype`.

## Details

Consider
[`cb_crop()`](https://mdsumner.github.io/cbr/reference/cb_crop.md) first
when the parent grid is large and the geometry occupies only part of it.

## Affine grids

For a burn from
[`cb_burn_index()`](https://mdsumner.github.io/cbr/reference/cb_burn_index.md)
the cells are rotated or sheared, which the `gis` format cannot describe
with a bounding box alone. The vector is still in grid order (row-major
from the top-left cell, one element per cell), `bbox` is the bounding
box of the grid's outline in world coordinates, and an extra element
`geotransform` holds the six GDAL numbers that define where each cell
is. Consumers that only read `bbox` and `dim` will place the raster
approximately but with square, unrotated cells.

## Examples

``` r
if (FALSE) { # \dontrun{
library(wk)
sq <- wkt("POLYGON ((2.5 4.5, 6.5 4.5, 6.5 8.5, 2.5 8.5, 2.5 4.5))")
b <- cb_burn(sq, extent = c(0, 10, 0, 10), dimension = c(10, 10))
v <- cb_materialize(b)
cb_as_matrix(v)
} # }
```
