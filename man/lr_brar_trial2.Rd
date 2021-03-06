% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lr_brar_trial.R
\name{lr_brar_trial2}
\alias{lr_brar_trial2}
\title{Run a logistic regression BRAR trial with two domains.}
\usage{
lr_brar_trial2(
  n_seq,
  full_dat,
  XA,
  XB,
  M0,
  S0,
  impute = F,
  impute_sim = 100,
  early_t = 7,
  trans_prob = 0.75,
  perm_drop = F,
  rar = F,
  rar_control = F,
  hpar1 = 0.5,
  hpar2 = 0,
  h = function(m) hpar1 * (m/length(n_seq))^hpar2,
  epsilon = 0.05,
  best_ref = 0.975,
  best_ref_a = 0.975^(ncol(XA) - 1),
  best_ref_b = 0.975^(ncol(XB)),
  supr_ref = 1,
  mc_draw = 20000,
  ...
)
}
\arguments{
\item{n_seq}{Sequence of interim analyses}

\item{full_dat}{Potential outcomes data from `gen_potential_outcomes`}

\item{XA}{Design matrix for domain A (include reference group as intercept)}

\item{XB}{Design matrix for domain B (exclude reference group as intercept)}

\item{M0}{Prior mean}

\item{S0}{Prior variance}

\item{impute}{Impute missing outcomes using early event information}

\item{impute_sim}{Number of imputation draws}

\item{early_t}{Time to early outcome (should be less than value of `delay` in `gen_potential_outcomes`)}

\item{trans_prob}{Probability of transitioning from no early event to a final event}

\item{perm_drop}{Permanently drop arms}

\item{rar}{Turn on RAR}

\item{rar_control}{Include control arm in RAR or fix}

\item{hpar1}{Parameter to control amount of RAR}

\item{hpar2}{Parameter to control scaling of RAR}

\item{h}{Function for RAR}

\item{supr_ref}{If dropping control, what threshold must an intervention arm achieve in terms of superiority}

\item{mc_draw}{Number of Monte Carlo draws form approximating posterior}

\item{...}{Other function arguments}

\item{fpar1}{Parameter to control lack of effect boundary (futility)}

\item{fpar2}{Parameter to control lack of effect boundary (futility)}

\item{f}{Function for utility boundary}

\item{gpar1}{Parameter to control most effective boundary (success)}

\item{gpar2}{Parameter to control most effective boundary (success)}

\item{g}{Function for most effective boundary}
}
\value{
A list of the trial quantities
}
\description{
Run a logistic regression BRAR trial with two domains.
}
