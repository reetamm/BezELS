% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_functions.R
\name{plot_bezier_euclidean}
\alias{plot_bezier_euclidean}
\title{Plot the limit set based on MCMC output in Euclidean coordinates and export median curve}
\usage{
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
}
\arguments{
\item{mcmc_samples}{Output object from mcmc_bezier function}

\item{x}{Original data in exponential margins}

\item{above_thresh_marg}{Vector of TRUE/FALSE values, where TRUE indicates points selected for fitting the truncated Gamma model}

\item{theta}{True value of copula dependence parameter (if known)}

\item{copula}{True parametric copula form (if known)}

\item{thin.by}{Thin the posterior for plotting - should ideally correspond to the thin.by parameter of the plot_bezier_polar function}

\item{medcurve}{Median curve output from the plot.bezier.polar function}

\item{plottitle}{Title of the plot}

\item{plot.truth}{TRUE/FALSE for whether the true limit set should be plotted (if known)}
}
\description{
Plot the limit set based on MCMC output in Euclidean coordinates and export median curve
}
