# DE-9IM matrix of any cell

Look up the DE-9IM intersection matrix of a polygon against grid cells,
for a burn with a `relations` table (see
[`cb_relate()`](https://mdsumner.github.io/cbr/reference/cb_relate.md)).
A cell its boundary touches has its matrix in `relations`; any other
cell is `"212FF1FF2"` (covered: wholly inside the polygon) if it is in
one of the polygon's runs, and `"FF2FF1212"` (disjoint) if not.

## Usage

``` r
cb_matrix(x, row, col, id)
```

## Arguments

- x:

  A `cb_burn` object with a `relations` table.

- row, col:

  Cell indices (1-based, row 1 at the top), recycled to a common length.

- id:

  Geometry ids, recycled with `row` and `col`. May be omitted when the
  burn has a single geometry.

## Value

A character vector of nine-character DE-9IM strings, the polygon's
interior, boundary and exterior (rows) against the cell's (columns).
`NA` for cells outside the grid, ids not in the burn, and ids of lines
and points (which have no polygon matrix).

## See also

[`cb_cells_where()`](https://mdsumner.github.io/cbr/reference/cb_cells_where.md)
to select cells by predicate, and `cb_materialize(value = "intersects")`
and friends for dense masks.

## Examples

``` r
g <- wk::wkt("POLYGON ((1 1, 4 1, 4 4, 1 4, 1 1))")
b <- cb_burn(g, extent = c(0, 6, 0, 6), dimension = c(6, 6), relate = TRUE)
## a cell deep inside, one below the bottom edge (sharing a wall), one
## touching a corner, and one far away
cb_matrix(b, row = c(4, 6, 6, 1), col = c(3, 3, 1, 6))
#> [1] "212FF1FF2" "FF2F11212" "FF2F01212" "FF2FF1212"
```
