# Coarse overview raster of a burn

Aggregate a
[`cb_burn()`](https://mdsumner.github.io/cbr/reference/cb_burn.md) onto
a grid of square tiles of `block` cells and return the mean coverage per
tile (in `[0, 1]` for polygons) in the `read_ds` layout used by
[`cb_materialize()`](https://mdsumner.github.io/cbr/reference/cb_materialize.md).
Edge tiles that are cut off by the grid boundary are normalized by their
actual cell count. Cost is proportional to the number of records plus
the number of tiles: a run that spans several tiles contributes partial
lengths to its two end tiles and the fully covered tiles between them
are filled with a difference array, so nothing is ever expanded to cells
whatever the block size.

## Usage

``` r
cb_overview(x, tiles = 256L, block = NULL)
```

## Arguments

- x:

  A `cb_burn` object.

- tiles:

  Target number of tiles along the longer axis; ignored if `block` is
  given.

- block:

  Tile size in cells. `block = 1` gives the exact coverage.

## Value

A numeric vector with a `gis` attribute (`datatype = "Float32"`), and
additionally attribute `block` giving the tile size used. Untouched
tiles are `NA`.
