% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/util_functions.R
\name{diff_beta_var}
\alias{diff_beta_var}
\title{Variance of a added beta random variables}
\usage{
diff_beta_var(a, b)
}
\arguments{
\item{a}{Collection of shape 1 par}

\item{b}{Collection of shape 2 par}
}
\description{
Always taken in reference to control, e.g. `(a[1], b[1])`.
}
\examples{
diff_beta_var(c(5, 5, 4), c(5, 5, 4))
}
