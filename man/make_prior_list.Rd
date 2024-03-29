% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/prior_specification.R
\name{make_prior_list}
\alias{make_prior_list}
\title{Make a list of priors to be passed to the SIR function.}
\usage{
make_prior_list(
  r_max = make_prior(runif, 0, 0.118),
  K = make_prior(use = FALSE),
  var_N = make_prior(0),
  N_obs = make_prior(runif, 500, 40000),
  add_CV = make_prior(use = FALSE),
  add_VAR_IA = make_prior(use = FALSE),
  catch_sample = make_prior(runif, 0, 1),
  z = make_prior(2.39),
  Pmsy = make_prior(use = FALSE),
  P50 = make_prior(0),
  q_IA1 = make_prior(use = FALSE),
  q_IA2 = make_prior(0),
  q_count = make_prior(use = FALSE)
)
}
\arguments{
\item{r_max}{Population growth rate; defaults to Uniform(0, 0.106).}

\item{K}{Carrying capacity, defaults to unused, and solution found using
recent observation and sampled \code{r_max}.}

\item{var_N}{Variance of lognormal (logmean = 0) population process error: defaults to 0 where no process error occurs.}

\item{N_obs}{Prior on a recent abundance estimate. Defaults to Uniform(500,
20,000).}

\item{add_CV}{Defaults to unused. Additional variability.}

\item{add_VAR_IA}{Defaults to unused. Additional variability on IA}

\item{catch_sample}{Defaults to unused. Samples between the minimum and maximum catch values.}

\item{z}{Defaults to constant 2.39. Shape parameter for generalized logistic
population dynamics function. Both z and Pmsy are
confounded and only one can be used. If \code{use = FALSE},
 Pmsy is used.}

\item{Pmsy}{Parameter that determines the level of depletion at MSY.
Both z and Pmsy are confounded and only one can be used.
 If \code{use = FALSE}, z is used.}

\item{P50}{Parameter that determines the level of depensation.}

\item{q_IA1}{Defaults to unused. Prior on annual slope for q for indices of abundance. If
\code{use = FALSE} for q_IA1, an analytic solution for q is used.}

\item{q_count}{Defaults to unused. Prior for q on counts.}
}
\value{
A named list containing each of the specified priors in a form that
  can be used by the SIR function.
}
\description{
Make a list of priors to be passed to the SIR function.
}
