# Generate data with unit exponential margins from 4 different bivariate copulas

Generate data with unit exponential margins from 4 different bivariate
copulas

## Usage

``` r
gen_data_exp(n, theta, copula = NULL, tau = 0.75)
```

## Arguments

- n:

  Number of observations

- theta:

  Dependence parameter (between 0 and 1)

- copula:

  `l` for logistic, `il` for inverted-logist, `al` for asymmetric
  logistic, `g` for Gaussian

- tau:

  Quantile threshold using marginal thresholding

## Value

x n x 2 matrix of data with unit exponential marginal distributions

r Vector of radii

w Vector of angles

r_0_marg truncation point for truncated Gamma distribution for all n
data points

above_thresh_marg TRUE/FALSE vector of whether points are above the
marginal threshold

data_marg_r_0 List with n, r, w, and r_0 for data points where
above_thresh_marg=T

## Examples

``` r
simdata <- gen_data_exp(n = 600, theta = 0.8, tau=0.75, copula = 'g')
summary(simdata$x)
#>        V1                 V2         
#>  Min.   :0.005722   Min.   :0.00232  
#>  1st Qu.:0.263048   1st Qu.:0.26039  
#>  Median :0.643098   Median :0.68423  
#>  Mean   :0.925769   Mean   :0.93596  
#>  3rd Qu.:1.240948   3rd Qu.:1.29205  
#>  Max.   :6.279349   Max.   :6.81087  
```
