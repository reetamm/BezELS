# Compute \\\lambda(\omega)\\ based on Bezier spline control points, usually a posterior draw obtained from `fit_mcmc_bezier`

Compute \\\lambda(\omega)\\ based on Bezier spline control points,
usually a posterior draw obtained from `fit_mcmc_bezier`

## Usage

``` r
lambda_empirical(control_points, omega)
```

## Arguments

- control_points:

  A 9x1 vector of Bezier spline control points

- omega:

  scalar between (0,1)

## Value

scalar value of \\\lambda(\omega)\\
