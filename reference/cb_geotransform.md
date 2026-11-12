# Build a GDAL-style geotransform

A geotransform is the six numbers `c(x0, dx, rx, y0, ry, dy)` that map a
pixel position (`col`, `row`, counted from the top-left corner of the
grid, 0-based) to world coordinates:

## Usage

``` r
cb_geotransform(origin = c(0, 0), res = 1, rotation = 0, shear = 0)
```

## Arguments

- origin:

  Position of the top-left corner of the grid, `c(x, y)`.

- res:

  Cell size: one value, or `c(width, height)`, in world units
  (positive).

- rotation:

  Angle in degrees, counter-clockwise, of the direction in which columns
  increase, measured from the +x axis.

- shear:

  Extra angle in degrees by which the direction in which rows increase
  is turned from being perpendicular to the columns (0 leaves the cells
  rectangular).

## Value

A numeric vector of length 6 in GDAL order.

## Details

    x = x0 + col * dx + row * rx
    y = y0 + col * ry + row * dy

This helper builds one from an origin, cell sizes, a rotation and a
shear, so that rotated and sheared grids can be written down without
doing the trigonometry by hand. A north-up grid is
`rotation = 0, shear = 0`, which gives `c(x0, res, 0, y0, 0, -res)`.

## Examples

``` r
cb_geotransform(c(100, 200), res = 10)
#> [1] 100  10   0 200   0 -10
cb_geotransform(c(0, 0), res = c(2, 1), rotation = 30)
#> [1]  0.0000000  1.7320508  0.5000000  0.0000000  1.0000000 -0.8660254
```
