# Burn WKB geometries onto a regular grid.

Low-level worker behind
[`cb_burn()`](https://mdsumner.github.io/cbr/reference/cb_burn.md).
`wkb` is a list of raw vectors (one WKB blob per geometry, as produced
by [`wk::wkb()`](https://paleolimbot.github.io/wk/reference/wkb.html)).
The grid is defined by its extent (`xmin`, `ymin`, `xmax`, `ymax`) and
dimensions (`ncol`, `nrow`), with row 1 at the top (`ymax` edge).
`coverage = TRUE` uses exact coverage fractions; `FALSE` uses the faster
cell-centre (approx) rule.

## Usage

``` r
cb_burn_wkb(
  wkb,
  xmin,
  ymin,
  xmax,
  ymax,
  ncol,
  nrow,
  coverage,
  parity,
  threads,
  id_offset,
  relate
)
```

## Arguments

- wkb:

  List of raw vectors, one WKB blob per geometry.

- xmin, ymin, xmax, ymax:

  Grid extent.

- ncol, nrow:

  Grid dimensions.

- coverage:

  Logical; exact coverage fractions if `TRUE`, cell-centre rule if
  `FALSE`.

- parity:

  Logical; if `TRUE` give the C++ core's records exactly: one run per
  fully covered boundary cell (no coalescing) and per-part records for
  multi-part geometries (no merging).

- threads:

  Integer; number of threads, `0` for all available. With more than one
  the input is parsed up front, and the crate uses at most one thread
  per 8192 units of work, so small burns stay serial.

- id_offset:

  Integer, non-negative; added to each geometry's position, so the `id`
  (and `geom_index`) columns are `id_offset` plus the 1-based input
  position.

- relate:

  Logical; if `TRUE` also compute the DE-9IM relations of each polygon
  against the cells its boundary touches (see
  [`cb_relate()`](https://mdsumner.github.io/cbr/reference/cb_relate.md)).

## Value

A list of five lists (`runs`, `edges`, `lines`, `points`, `notes`), each
a set of equal-length columns, and `relations` (the relations table if
`relate`, otherwise `NULL`).

## Details

Returns a list of five column-oriented tables: `runs`, `edges`, `lines`,
`points`, and `notes`. The crate emits 0-based indices, half-open run
ranges and `id = k`; this function converts to the R convention of
1-based `row`/`col`, inclusive `col_end`, and `id`/`geom_index` equal to
the 1-based input position.
