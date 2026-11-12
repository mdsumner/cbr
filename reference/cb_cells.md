# Cells of a burn as polygons in world coordinates

One quadrilateral per record cell, in the grid's world coordinates: the
cells of every run expanded, then the edge, line and point cells. For a
burn from
[`cb_burn_index()`](https://mdsumner.github.io/cbr/reference/cb_burn_index.md)
the quadrilaterals are rotated or sheared as the grid is; for a regular
[`cb_burn()`](https://mdsumner.github.io/cbr/reference/cb_burn.md) they
are rectangles. Intended for plotting
([`ggplot2::geom_polygon()`](https://ggplot2.tidyverse.org/reference/geom_polygon.html))
and for export.

## Usage

``` r
cb_cells(x, max_cells = 1e+05)
```

## Arguments

- x:

  A `cb_burn` object.

- max_cells:

  Refuse to build more cells than this (runs can cover very many cells),
  to avoid exhausting memory by accident.

## Value

A data frame with four rows per cell and columns `cell` (cell number,
the polygon group), `x`, `y` (the corner, in order around the cell),
`row`, `col`, `id`, `kind` (`run`, `edge`, `line` or `point`) and
`value` (1 for runs and points, the coverage fraction for edges, the
length for lines).

## Examples

``` r
b <- cb_burn(wk::wkt("POLYGON ((1.5 1.5, 4.5 1.5, 4.5 3.5, 1.5 1.5))"),
             extent = c(0, 6, 0, 5), dimension = c(6, 5))
head(cb_cells(b), 8)
#>   cell x y row col id kind      value
#> 1    1 3 4   2   4  1 edge 0.02083333
#> 2    1 4 4   2   4  1 edge 0.02083333
#> 3    1 4 3   2   4  1 edge 0.02083333
#> 4    1 3 3   2   4  1 edge 0.02083333
#> 5    2 4 4   2   5  1 edge 0.16666667
#> 6    2 5 4   2   5  1 edge 0.16666667
#> 7    2 5 3   2   5  1 edge 0.16666667
#> 8    2 4 3   2   5  1 edge 0.16666667
```
