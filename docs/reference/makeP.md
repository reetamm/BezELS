# Make a 7x2 matrix that defines a Bezier spline, based on parameters and constraints

Make a 7x2 matrix that defines a Bezier spline, based on parameters and
constraints

## Usage

``` r
makeP(theta, pmix0, pmix1, pmix2, pmix3)
```

## Arguments

- theta:

  Vector of 9 control points that define the shape of the spline

- pmix0:

  Unit mass probability of p0y and p6x

- pmix1:

  Unit mass probability for p1x and p5y

- pmix2:

  Unit mass probability for p2x and p4y

- pmix3:

  Unit mass probability for p3

## Value

P 7x2 matrix that defines a Bezier spline comprised of 3 Bezier curves
