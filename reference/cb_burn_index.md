# Burn geometries onto an affine (rotated or sheared) grid

Rasterize onto a grid defined by a geotransform rather than by an
axis-aligned extent. The geometries are moved into the grid's index
space (fractional column and row) by the inverse of the geotransform and
burned on a lattice of unit cells; the records are then already in the
grid's own `(row, col)`. Because the map is affine, straight edges stay
straight and the coverage fractions are exact: each is the fraction of
the (possibly rotated, sheared or non-square) cell's area covered by the
geometry, up to the single precision of the `edges` table.

## Usage

``` r
cb_burn_index(
  geom,
  geotransform,
  dimension,
  coverage = TRUE,
  crs = NA_character_,
  parity = FALSE,
  threads = 1L,
  id_offset = 0L,
  relate = FALSE
)
```

## Arguments

- geom:

  Geometry input, anything
  [`wk::as_wkb()`](https://paleolimbot.github.io/wk/reference/wkb.html)
  understands, in the world coordinates the geotransform maps to.

- geotransform:

  Numeric vector of length 6, `c(x0, dx, rx, y0, ry, dy)`, as in GDAL;
  see
  [`cb_geotransform()`](https://mdsumner.github.io/cbr/reference/cb_geotransform.md)
  to build one. The grid's top-left corner is `(x0, y0)`.

- dimension:

  Grid dimensions as `c(ncol, nrow)`.

- coverage:

  Logical. If `TRUE` (default) use exact coverage fractions (boundary
  cells appear in `edges`). If `FALSE` use the faster cell-centre
  ("approx", fasterize-style) rule, which emits runs only.

- crs:

  Optional coordinate reference system as a WKT string (or anything
  coercible to character). Carried through to
  [`cb_materialize()`](https://mdsumner.github.io/cbr/reference/cb_materialize.md)
  output.

- parity:

  Logical. By default a multi-part geometry gets one record per cell
  (parts' coverage is summed, and a cell the parts together fill becomes
  a run), and adjacent runs on a row are coalesced. If `TRUE`, give the
  records of the C++ 'controlledburn' core instead: one run per fully
  covered boundary cell and separate records for each part.

- threads:

  Number of threads, `1` (default) for none, `0` for all available. The
  result is identical to the serial burn. Threads only pay for large
  inputs (the crate uses at most one per 8192 units of work), and with
  more than one thread the input is parsed up front, which uses more
  memory.

- id_offset:

  Non-negative whole number added to each geometry's position to give
  its `id` (and `geom_index` in `notes`), so `id = id_offset + k` for
  the `k`th input geometry. Burn consecutive slices of a collection with
  offsets `0`, `n1`, `n1 + n2`, ... and join them with
  [`cb_combine()`](https://mdsumner.github.io/cbr/reference/cb_combine.md)
  to get exactly the burn of the whole collection. Recorded in attribute
  `id_offset`.

- relate:

  Logical. If `TRUE`, also compute the DE-9IM relation of each polygon
  to the cells its boundary touches and return it as a sixth table,
  `relations`; see
  [`cb_relate()`](https://mdsumner.github.io/cbr/reference/cb_relate.md),
  which does the same for an existing burn.

## Value

A `cb_burn` object; see
[`cb_burn()`](https://mdsumner.github.io/cbr/reference/cb_burn.md).

## Details

The result is an ordinary
[`cb_burn()`](https://mdsumner.github.io/cbr/reference/cb_burn.md)
object, so
[`cb_crop()`](https://mdsumner.github.io/cbr/reference/cb_crop.md),
[`cb_combine()`](https://mdsumner.github.io/cbr/reference/cb_combine.md),
[`cb_materialize()`](https://mdsumner.github.io/cbr/reference/cb_materialize.md)
and [`summary()`](https://rdrr.io/r/base/summary.html) work on it. Its
`extent` attribute is the index lattice `c(0, ncol, 0, nrow)` (so
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) draws in index
space), and the world definition is kept in attributes:

- `geotransform`:

  the six numbers, as given.

- `world_extent`:

  `c(xmin, xmax, ymin, ymax)` of the grid's outline in world coordinates
  (its bounding box, for a rotated grid).

- `crs`:

  as given.

Use [`cb_cells()`](https://mdsumner.github.io/cbr/reference/cb_cells.md)
to get the cells as polygons in world coordinates, for plotting or
export.

## Line lengths

The `length` column of `lines` is in world units when the cells are
similar to squares (equal sides at right angles, whatever the rotation).
For other affine grids the length of a piece depends on its direction,
and the crate reports lengths in cell units, which is what is returned,
with a warning.

## Examples

``` r
library(wk)
gt <- cb_geotransform(c(0, 6), res = 1, rotation = 20)
sq <- wkt("POLYGON ((3 3, 7 3, 7 7, 3 7, 3 3))")
b <- cb_burn_index(sq, gt, dimension = c(12, 8))
b
#> <cb_burn>
#>   grid:   12 col x 8 row
#>   extent: x[0, 12] y[0, 8]
#>   affine: geotransform c( 0.000000,  0.939693,  0.342020,  6.000000,  0.342020, -0.939693)
#>   world:  x[0, 14.0125] y[-1.51754, 10.1042] (extent above is index space)
#>   mode:   coverage
#>   tables: runs=3, edges=19, lines=0, points=0, notes=0 
attr(b, "world_extent")
#> [1]  0.000000 14.012473 -1.517541 10.104242
```
