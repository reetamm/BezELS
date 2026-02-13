# Plot the limit set based on MCMC output in polar coordinates and export median curve

Plot the limit set based on MCMC output in polar coordinates and export
median curve

## Usage

``` r
plot_bezier_polar(
  mcmc_samples,
  x,
  r_0_marg,
  above_thresh_marg,
  theta = NULL,
  copula = NULL,
  thin.by = 10,
  plot.fig = T,
  plottitle = "Polar coordinates",
  plot.truth = T
)
```

## Arguments

- mcmc_samples:

  Output object from mcmc_bezier function

- x:

  Original data in exponential margins

- r_0_marg:

  Vector of truncation thresholds for all points

- above_thresh_marg:

  Vector of TRUE/FALSE values, where TRUE indicates points selected for
  fitting the truncated Gamma model

- theta:

  True value of copula dependence parameter (if known)

- copula:

  True parametric copula form (if known)

- thin.by:

  Thin the posterior for plotting (and median curve calculation)

- plot.fig:

  TRUE/FALSE for whether the limit set should be plotted

- plottitle:

  Title of plot (only used if plot.fig = T)

- plot.truth:

  TRUE/FALSE for whether the true limit set should be plotted (if known)

## Value

Row of thinned mcmc_samples which corresponds to the median curve
