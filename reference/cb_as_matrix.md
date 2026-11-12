# Convert a materialized raster to a matrix

Reshape the vector returned by
[`cb_materialize()`](https://mdsumner.github.io/cbr/reference/cb_materialize.md)
(or by `gdalraster::read_ds()`) into a matrix with `nrow` rows and
`ncol` columns, oriented as the grid appears on a map (row 1 at the
top). The `gis` attribute is dropped.

## Usage

``` r
cb_as_matrix(x)
```

## Arguments

- x:

  A vector with a `gis` attribute holding `dim = c(ncol, nrow, ...)`.

## Value

A matrix.
