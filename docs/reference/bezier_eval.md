# Find a point on a 2D Bezier curve B(t) at a particular value of t

Find a point on a 2D Bezier curve B(t) at a particular value of t

## Usage

``` r
bezier_eval(p0, p1, p2, t)
```

## Arguments

- p0:

  Start control point

- p1:

  Intermediate control point

- p2:

  End control point

- t::

  Location on the curve

## Value

The x-coordinate of B(t)

## Examples

``` r
bezier_eval(p0 = c(0.4,0.1),p1 = c(0.9,0.9),p2 = c(0.1,0.4),t = 0)
#> [1] 0.4
```
