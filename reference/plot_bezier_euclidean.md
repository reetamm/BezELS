# Plot the limit set based on MCMC output in Euclidean coordinates and export median curve

Plot the limit set based on MCMC output in Euclidean coordinates and
export median curve

## Usage

``` r
plot_bezier_euclidean(
  mcmc_samples,
  x,
  above_thresh_marg,
  theta = NULL,
  copula = NULL,
  thin.by = 10,
  medcurve = NULL,
  plottitle = "Euclidean coordinates",
  plot.truth = T
)
```

## Arguments

- mcmc_samples:

  Output object from mcmc_bezier function

- x:

  Original data in exponential margins

- above_thresh_marg:

  Vector of TRUE/FALSE values, where TRUE indicates points selected for
  fitting the truncated Gamma model

- theta:

  True value of copula dependence parameter (if known)

- copula:

  True parametric copula form (if known)

- thin.by:

  Thin the posterior for plotting - should ideally correspond to the
  thin.by parameter of the plot_bezier_polar function

- medcurve:

  Median curve output from the plot.bezier.polar function

- plottitle:

  Title of the plot

- plot.truth:

  TRUE/FALSE for whether the true limit set should be plotted (if known)
