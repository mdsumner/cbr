# Summarize a burn without materializing it

Compute a compact description of a
[`cb_burn()`](https://mdsumner.github.io/cbr/reference/cb_burn.md)
object in time proportional to the number of records in its tables,
never to the number of grid cells. The pieces are:

## Usage

``` r
# S3 method for class 'cb_burn'
summary(object, tiles = 256L, ...)
```

## Arguments

- object:

  A `cb_burn` object.

- tiles:

  Target number of tiles along the longer grid axis for the overview.
  The tile size in cells is `ceiling(max(dimension) / tiles)`.

- ...:

  Ignored.

## Value

A list of class `summary.cb_burn`.

## Details

- `overview`:

  a coarse raster of mean coverage per tile, in the same `read_ds`
  layout as
  [`cb_materialize()`](https://mdsumner.github.io/cbr/reference/cb_materialize.md)
  (see
  [`cb_overview()`](https://mdsumner.github.io/cbr/reference/cb_overview.md)).
  When the grid is small enough that one tile is one cell, this is the
  exact coverage.

- `row_marginal`, `col_marginal`:

  exact covered-cell totals per grid row and per grid column, at full
  resolution.

- `content`:

  the extent and 1-based row/col range of touched cells, or `NULL` for
  an empty burn.
  [`cb_crop()`](https://mdsumner.github.io/cbr/reference/cb_crop.md)
  uses this as its default window.

- `area`:

  per-geometry totals: `interior` (run cells), `boundary` (summed edge
  fractions), `lines` (summed length), `points` (count).

- `fractions`:

  quantiles of the edge coverage fractions.

- `records`:

  row counts of the five tables.

Coverage here means: run cells 1 each, edge cells their fraction, line
and point cells 1 each (presence).
