# Burn geometries onto a regular grid

Rasterize polygons, lines, and points onto a regular grid using the
`controlledburn` Rust crate. Polygon interiors are returned as
horizontal runs of fully covered cells, and (in coverage mode) boundary
cells are returned with their exact coverage fraction.

## Usage

``` r
cb_burn(
  geom,
  extent,
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

  Geometry input. Anything that
  [`wk::as_wkb()`](https://paleolimbot.github.io/wk/reference/wkb.html)
  understands, e.g. a `wk_wkb` vector, an `sf`/`sfc` column, or a list
  of raw WKB blobs.

- extent:

  Grid extent as `c(xmin, xmax, ymin, ymax)` (the ordering used by
  [`terra::ext()`](https://rspatial.github.io/terra/reference/ext.html)).
  If missing, the bounding box of `geom` is used (via
  [`wk::wk_bbox()`](https://paleolimbot.github.io/wk/reference/wk_bbox.html)).

- dimension:

  Grid dimensions as `c(ncol, nrow)`. If missing, `ncol` is 256 and
  `nrow` is computed from it so that cells are square in the units of
  `extent` (at least 1 row). A single value is taken as `ncol`, with
  `nrow` computed the same way.

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

An object of class `cb_burn`: a named list of five data frames with the
grid definition stored in attributes (`extent`, `dimension`, `crs`,
`coverage`) and the number of input geometries in `n_geom`; with
`relate = TRUE` it also has a sixth table, `relations`. Indices are
1-based and row 1 is the top row; geometry `k` (input position) gets
`id = k`.

- `runs`:

  `row`, `col_start`, `col_end`, `id`: fully covered cells.

- `edges`:

  `row`, `col`, `fraction`, `id`: boundary cells, fraction in (0, 1).

- `lines`:

  `row`, `col`, `length`, `id`: line length within the cell, CRS units.

- `points`:

  `row`, `col`, `id`: one row per point inside the grid.

- `notes`:

  `geom_index`, `message`: non-fatal per-geometry problems.

## Examples

``` r
if (FALSE) { # \dontrun{
library(wk)
sq <- wkt("POLYGON ((2.5 4.5, 6.5 4.5, 6.5 8.5, 2.5 8.5, 2.5 4.5))")
cb_burn(sq, extent = c(0, 10, 0, 10), dimension = c(10, 10))

## defaults: bbox of the geometry, 256 columns, square cells
cb_burn(sq)
} # }
```
