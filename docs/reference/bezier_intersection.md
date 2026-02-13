# Find the intersection of a 2D Bezier curve B(t) and a line passing through the origin with slope m

Find the intersection of a 2D Bezier curve B(t) and a line passing
through the origin with slope m

## Usage

``` r
bezier_intersection(p0, p1, p2, m)
```

## Arguments

- p0:

  Start control point

- p1:

  Intermediate control point

- p2:

  End control point

- m:

  Slope of straight line passing through the origin

## Value

Value of t where the line intersects the curve B(t)

## Examples

``` r
bezier_intersection(p0 = c(0.4,0.1),p1 = c(0.9,0.9),p2 = c(0.1,0.4),m = 1)
#> [1] 0.5
```
