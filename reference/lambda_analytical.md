# Evaluate the analytical value of \\\lambda(\omega)\\ for bivariate Gaussian, logistic, asymmetric logistic, and inverted logistic copulas

Evaluate the analytical value of \\\lambda(\omega)\\ for bivariate
Gaussian, logistic, asymmetric logistic, and inverted logistic copulas

## Usage

``` r
lambda_analytical(theta, omega, dep = NULL)
```

## Arguments

- theta:

  Strength of dependence, between (0,1)

- omega:

  scalar between (0,1)

- dep:

  Copula; `g` for Gaussian, `l` for logistic, `il` for inverted
  logistic, `al` for asymmetric logistic

## Value

scalar value of \\\lambda(\omega)\\
