# DE-9IM of polygons against grid cells, for an existing burn.

Low-level worker behind
[`cb_relate()`](https://mdsumner.github.io/cbr/reference/cb_relate.md).
`runs` is the burn's runs table (R's conventions: 1-based, inclusive
`col_end`); `wkb` must be the geometries that were burned, in the same
grid coordinates, and `coverage` and `id_offset` the burn's. Returns the
relations table: `row`, `col`, `id`, `matrix`.

## Usage

``` r
cb_relate_wkb(
  wkb,
  runs,
  xmin,
  ymin,
  xmax,
  ymax,
  ncol,
  nrow,
  coverage,
  id_offset
)
```

## Arguments

- wkb:

  List of raw vectors, one WKB blob per geometry.

- runs:

  Data frame with integer columns `row`, `col_start`, `col_end`, `id`.

- xmin, ymin, xmax, ymax:

  Grid extent.

- ncol, nrow:

  Grid dimensions.

- coverage:

  Logical; the burn's mode.

- id_offset:

  Integer; the burn's `id_offset`.

## Value

A list of equal-length columns `row`, `col`, `id`, `matrix`.
