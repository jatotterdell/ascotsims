% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/util_functions.R
\name{utility_ent}
\alias{utility_ent}
\title{Calculate change in utility using entropy}
\usage{
utility_ent(n, y, a = 1, b = 1, ...)
}
\arguments{
\item{n}{Sample size in each arm}

\item{y}{Response in each arm}

\item{a}{Prior shape}

\item{b}{Prior shape}

\item{...}{Other function arguments}
}
\description{
Calculate change in utility using entropy
}
\examples{
ascotsims:::utility_ent(rep(0, 4), rep(0, 4))
ascotsims:::utility_ent(c(20, 19, 18, 17), c(8,9,10,11))
}
