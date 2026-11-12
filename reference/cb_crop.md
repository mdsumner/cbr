# Crop a burn to a sub-window of its grid

Snap an arbitrary extent to the cell boundaries of a
[`cb_burn()`](https://mdsumner.github.io/cbr/reference/cb_burn.md) grid
and subset the burn tables to that window. This is a pure index
operation: no geometry is re-rasterized, runs that straddle the window
edge are clipped, and row/column indices are re-based so that the result
is a valid `cb_burn` on the smaller grid. Use it before
[`cb_materialize()`](https://mdsumner.github.io/cbr/reference/cb_materialize.md)
to avoid allocating a dense array for the whole parent grid.

## Usage

``` r
cb_crop(x, extent, snap = c("out", "in", "near"))
```

## Arguments

- x:

  A `cb_burn` object.

- extent:

  Crop extent as `c(xmin, xmax, ymin, ymax)`, in the CRS of the grid
  (for a burn from
  [`cb_burn_index()`](https://mdsumner.github.io/cbr/reference/cb_burn_index.md),
  in its index space: x is the column position and y counts rows up from
  the bottom). It is snapped to the parent grid and intersected with it.
  If missing, the burn is cropped to its content: the cells touched by
  any geometry (see `content` in
  [`summary.cb_burn()`](https://mdsumner.github.io/cbr/reference/summary.cb_burn.md)).

- snap:

  How to align `extent` to cell boundaries. `"out"` (default) expands to
  the enclosing cells, `"in"` shrinks to fully contained cells, `"near"`
  moves each edge to the nearest cell boundary.

## Value

A `cb_burn` object on the cropped grid. The `extent` and `dimension`
attributes describe the new grid; `crs` and `coverage` are carried
through. The parent grid is recorded in attribute `parent` as a list
with `extent`, `dimension`, and the 0-based cell `offset`
(`c(col, row)`) of the window within it. An error is raised if the
snapped window is empty.

## Examples

``` r
if (FALSE) { # \dontrun{
library(wk)
sq <- wkt("POLYGON ((2.5 4.5, 6.5 4.5, 6.5 8.5, 2.5 8.5, 2.5 4.5))")
b <- cb_burn(sq, extent = c(0, 10, 0, 10), dimension = c(10, 10))
cb_crop(b, c(2.2, 5.1, 4.9, 9))   # snapped to x[2, 6] y[4, 9]
} # }
```
