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
#>  Min.   :0.005073   Min.   :0.000534  
#>  1st Qu.:0.239599   1st Qu.:0.263604  
#>  Median :0.652507   Median :0.662748  
#>  Mean   :0.931075   Mean   :0.913779  
#>  3rd Qu.:1.275914   3rd Qu.:1.206342  
#>  Max.   :6.457210   Max.   :4.989393  
```
