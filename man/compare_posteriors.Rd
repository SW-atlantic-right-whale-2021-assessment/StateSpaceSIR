% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_sir.R
\name{compare_posteriors}
\alias{compare_posteriors}
\title{Function to compare posteriors}
\usage{
compare_posteriors(
  SIR,
  model_names = NULL,
  file_name = NULL,
  bayes_factor = NULL,
  reference_sir = TRUE,
  model_average = TRUE,
  years = NULL
)
}
\arguments{
\item{SIR}{list of SIR objects}

\item{model_names}{names of sir objects}

\item{file_name}{name of a file to identified the files exported by the
function. If NULL, does not save.}

\item{bayes_factor}{Optional. Vector of bayesfactors of length SIR}

\item{reference_sir}{Default = TRUE, reference SIR is in \code{SIR}, should be first if so.}

\item{model_average}{Default = TRUE, model average SIR is in \code{SIR}, should be last if so.}

\item{years}{Optional vector of years to compare.}
}
\description{
Function to compare posteriors
}
