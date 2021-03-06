% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/util_functions.R
\name{beta_prob_supr_approx}
\alias{beta_prob_supr_approx}
\title{Probability that each arm superior to reference.}
\usage{
beta_prob_supr_approx(a, b, reverse = FALSE)
}
\arguments{
\item{a}{First value is reference}

\item{b}{First value is reference}

\item{reverse}{Reverse direction of superiority}
}
\description{
Calculates event that one Beta RV larger (or smaller if reverse) to another Beta RV by normal approximation.
}
\examples{
beta_prob_supr_approx(c(5, 10, 15), c(15, 10, 5))
beta_prob_supr_approx(c(5, 10, 15), c(15, 10, 5), reverse = TRUE)
}
