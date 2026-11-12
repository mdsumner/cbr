# Combine burns of consecutive slices of a collection

Concatenate
[`cb_burn()`](https://mdsumner.github.io/cbr/reference/cb_burn.md)
results that share a grid and whose geometry ids do not overlap, as when
a collection is burned in chunks (one per worker, say) with `id_offset`
set to the number of geometries before each chunk. The result is the
same, record for record, as burning the whole collection at once: each
geometry's records are produced independently, so joining the chunks in
id order changes nothing.

## Usage

``` r
cb_combine(...)
```

## Arguments

- ...:

  `cb_burn` objects, or a single list of them. They may be given in any
  order; the tables are joined in order of `id_offset`.

## Value

A `cb_burn` object on the shared grid. Its `id_offset` is the lowest of
the inputs' and `n_geom` spans to the highest id, so if the chunks leave
a gap in the ids those geometries count as empty. `parent` is kept when
all inputs agree on it.

## Details

All inputs must have identical `extent`, `dimension`, `crs`, `coverage`
and (for
[`cb_burn_index()`](https://mdsumner.github.io/cbr/reference/cb_burn_index.md)
burns) `geotransform`, and their id ranges
`(id_offset, id_offset + n_geom]` must not overlap, and either all or
none must have a `relations` table; otherwise an error is raised rather
than silently merging different geometries under one id.

## Examples

``` r
if (FALSE) { # \dontrun{
library(wk)
g <- wkt(c("POLYGON ((1 1, 4 1, 4 4, 1 4, 1 1))",
           "POLYGON ((5 5, 9 5, 9 9, 5 9, 5 5))"))
grid <- list(extent = c(0, 10, 0, 10), dimension = c(10, 10))
a <- cb_burn(g[1], extent = grid$extent, dimension = grid$dimension)
b <- cb_burn(g[2], extent = grid$extent, dimension = grid$dimension,
             id_offset = 1)
cb_combine(a, b)
} # }
```
