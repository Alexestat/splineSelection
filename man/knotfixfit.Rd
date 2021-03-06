% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/knotfixfit.R
\name{knotfixfit}
\alias{knotfixfit}
\title{Knots locations selection}
\usage{
knotfixfit(
  x,
  y,
  nknots,
  delta,
  degree,
  b.knots = c(min(x), max(x)),
  bspline = TRUE
)
}
\arguments{
\item{x}{Vector of function arguments.}

\item{y}{Vector of function values.}

\item{nknots}{Number of internal knots.}

\item{delta}{Minimal distance value between two consecutive knots.}

\item{degree}{Degree of the spline basis.}

\item{b.knots}{External knots vector. Default is the (min(x),max(x)).}

\item{bspline}{Logical value. If TRUE, B-spline is fitted. Else, natural splines.}
}
\value{
A list with selected internal knots, mean square error and fitted model.
}
\description{
Provides internal knots locations in spline-based models for a fixed number of knots.
}
