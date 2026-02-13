# Evaluate the analytical value of \\\tau_1(\delta)\\ for bivariate Gaussian, logistic, asymmetric logistic, and inverted logistic copulas

Evaluate the analytical value of \\\tau_1(\delta)\\ for bivariate
Gaussian, logistic, asymmetric logistic, and inverted logistic copulas

## Usage

``` r
tau_analytical(theta, delta, dep = NULL)
```

## Arguments

- theta:

  Strength of dependence, between (0,1)

- delta:

  scalar between (0,1)

- dep:

  Copula; `g` for Gaussian, `l` for logistic, `il` for inverted
  logistic, `al` for asymmetric logistic

## Value

scalar value of \\\tau_1(\delta)\\
