# DE-9IM relations of polygons to the cells their boundaries touch

A burn says how much of each cell a polygon covers. The DE-9IM matrix
says how the polygon and a cell (a closed box) meet: interiors,
boundaries and exteriors, each against each. For almost every cell the
burn already answers that: a cell the polygon's boundary does not reach
is either wholly inside it (a run) or wholly outside (no record).
`cb_relate()` computes the matrix for the remaining cells, those whose
closed box the boundary touches, and stores them in a sixth table,
`relations`. That includes cells the polygon only touches along a wall
or at a corner, which have no record in a burn at all.

## Usage

``` r
cb_relate(x, geom)
```

## Arguments

- x:

  A `cb_burn` object that has not been cropped.

- geom:

  The geometries that were burned into `x`, in the same coordinates and
  order: anything
  [`wk::as_wkb()`](https://paleolimbot.github.io/wk/reference/wkb.html)
  understands. There must be `n_geom` of them. The burn does not keep
  its geometry, so it is passed again; passing different geometry gives
  meaningless matrices.

## Value

`x`, with a `relations` table. Use `cb_burn(relate = TRUE)` to get the
same from the geometry in one step.

## Details

The matrix is read from the geometry, not from the burn's records. The
burn treats a coverage within 1e-6 of 0 or 1 as exactly 0 or 1, so a
cell holding a sliver of a polygon has no record and a cell holding a
sliver of its exterior is a run; the relations see both.

## The relations table

One row for each polygon and in-grid cell whose closed box the polygon's
boundary meets, sorted by `id`, `row` and `col`:

- `row`, `col`, `id`:

  as in the other tables: 1-based, row 1 at the top, `id` the position
  of the geometry plus `id_offset`.

- `matrix`:

  the DE-9IM matrix of the polygon against the cell as nine characters,
  row-major over the polygon's interior, boundary and exterior (rows)
  against the cell's (columns), each `F` or the dimension `0`, `1` or
  `2` of the intersection, for example `"212F11FF2"`.

Every other cell of a polygon is `"212FF1FF2"` if it is in one of that
polygon's runs and `"FF2FF1212"` otherwise, so the runs and the
relations together give the matrix of every cell. Lines, points and
empty geometries have no relations.

## Assumptions

Polygons are assumed valid (rings that do not cross, holes inside their
shell, parts of a multipolygon meeting only at points); the matrices of
invalid polygons are not defined. Contacts are decided with exact
predicates on the coordinates as given, with no snapping, so a polygon
offset from a cell wall by `1e-14` does not touch it. For a burn from
[`cb_burn_index()`](https://mdsumner.github.io/cbr/reference/cb_burn_index.md)
the geometry is moved into index space exactly as it was burned, and the
matrix is exact for the rotated or sheared cell; on grids burned in
index space by other front ends it is the index-space matrix.

## Examples

``` r
library(wk)
g <- wkt("POLYGON ((1 1, 3 1, 3 3, 1 3, 1 1))")
b <- cb_relate(cb_burn(g, extent = c(0, 4, 0, 4), dimension = c(4, 4)), g)
b$relations
#>    row col id    matrix
#> 1    1   1  1 FF2F01212
#> 2    1   2  1 FF2F11212
#> 3    1   3  1 FF2F11212
#> 4    1   4  1 FF2F01212
#> 5    2   1  1 FF2F11212
#> 6    2   2  1 212F11FF2
#> 7    2   3  1 212F11FF2
#> 8    2   4  1 FF2F11212
#> 9    3   1  1 FF2F11212
#> 10   3   2  1 212F11FF2
#> 11   3   3  1 212F11FF2
#> 12   3   4  1 FF2F11212
#> 13   4   1  1 FF2F01212
#> 14   4   2  1 FF2F11212
#> 15   4   3  1 FF2F11212
#> 16   4   4  1 FF2F01212
table(b$relations$matrix)
#> 
#> 212F11FF2 FF2F01212 FF2F11212 
#>         4         4         8 
```
