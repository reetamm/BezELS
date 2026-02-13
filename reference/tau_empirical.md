# Compute \\\tau_1(\delta)\\ based on Bezier spline control points, usually a posterior draw obtained from `fit_mcmc_bezier`

Compute \\\tau_1(\delta)\\ based on Bezier spline control points,
usually a posterior draw obtained from `fit_mcmc_bezier`

## Usage

``` r
tau_empirical(control_points, delta)
```

## Arguments

- control_points:

  A 9x1 vector of Bezier spline control points

- delta:

  scalar between (0,1)

## Value

scalar value of \\\tau_1(\delta)\\
