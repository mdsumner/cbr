# Cells where the burn's records and the exact topology disagree

The burn treats coverage within 1e-6 of 0 or 1 as exactly 0 or 1. That
is right for area, but it means a cell holding a sliver of a polygon can
have no record, and a cell holding a sliver of the polygon's exterior
can be in a run. The matrices in the `relations` table come from the
geometry, not the records, so comparing the two finds those cells.

## Usage

``` r
cb_slivers(x)
```

## Arguments

- x:

  A `cb_burn` object with a `relations` table (see
  [`cb_relate()`](https://mdsumner.github.io/cbr/reference/cb_relate.md)).

## Value

A data frame with one row per geometry and disagreeing cell, sorted by
`id`, `row` and `col`: `row`, `col`, `id`, `matrix` and `kind`, one of

- `"interior"`:

  the polygon shares area with the cell, but the burn has no record
  there (a sliver of the polygon, rounded to 0).

- `"exterior"`:

  the cell is in a run, but the polygon does not cover it (a sliver of
  the exterior, rounded to 1).

No rows when the records agree with the topology everywhere, which is
the usual case.

## Examples

``` r
## a square with a spike of area 5e-8 into the next cell: no record there
g <- wk::wkt("POLYGON ((1 1, 3 1, 3 2, 3.0000001 2.5, 3 3, 1 3, 1 1))")
b <- cb_burn(g, extent = c(0, 4, 0, 4), dimension = c(4, 4), relate = TRUE)
cb_slivers(b)
#>   row col id    matrix     kind
#> 1   2   4  1 212101212 interior
```
