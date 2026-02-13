# Evaluate the rate parameter of the truncated Gamma distribution

Evaluate the rate parameter of the truncated Gamma distribution

## Usage

``` r
gx(N = N, p = P, m = x[, 2]/x[, 1], p_l = x[, 1])
```

## Arguments

- N:

  Number of data points that the truncated Gamma distribution is to be
  fitted to

- p:

  7x2 matrix of points that define the Bezier spline

- m:

  Vector of slopes of length N; entries are slopes of lines connecting
  the origin to each point

- p_l:

  x-coordinate of each point

## Value

Vector of rate parameters of the truncated Gamma distribution
