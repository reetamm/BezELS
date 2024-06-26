% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bezier_functions.R
\name{bezier_intersection}
\alias{bezier_intersection}
\title{Find the intersection of a 2D Bezier curve B(t) and a line passing through the origin with slope m}
\usage{
bezier_intersection(p0, p1, p2, m)
}
\arguments{
\item{p0}{Start control point}

\item{p1}{Intermediate control point}

\item{p2}{End control point}

\item{m}{Slope of straight line passing through the origin}
}
\value{
Value of t where the line intersects the curve B(t)
}
\description{
Find the intersection of a 2D Bezier curve B(t) and a line passing through the origin with slope m
}
\examples{
bezier_intersection(p0 = c(0.4,0.1),p1 = c(0.9,0.9),p2 = c(0.1,0.4),m = 1)
}
