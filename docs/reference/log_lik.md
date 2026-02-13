# Evaluate the truncated Gamma likelihood for a vector of radii

Evaluate the truncated Gamma likelihood for a vector of radii

## Usage

``` r
log_lik(r, r_0, shape, rates)
```

## Arguments

- r:

  Vector of radii

- r_0:

  Vector of truncation points for each radii

- shape:

  Common shape parameter for the truncated Gamma distribution

- rates:

  Vector of rate parameters for each truncated Gamma distribution

## Value

Sum of the truncated Gamma log-likelihoods

## Examples

``` r
log_lik(r = 10, r_0 = 3, shape = 2, rates = 1)
#> [1] -2.584957
```
