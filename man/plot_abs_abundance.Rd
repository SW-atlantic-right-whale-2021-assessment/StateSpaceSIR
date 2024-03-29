% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_sir.R
\name{plot_abs_abundance}
\alias{plot_abs_abundance}
\title{OUTPUT FU#' OUTPUT FUNCTION}
\usage{
plot_abs_abundance(
  SIR,
  Reference = NULL,
  file_name = NULL,
  posterior_pred = TRUE,
  coolors = "#941B0C",
  coolors2 = "#104F55"
)
}
\arguments{
\item{SIR}{A fitted SIR model}

\item{Reference}{A fitted SIR model as the reference model}

\item{file_name}{name of a file to identified the files exported by the
function. If NULL, does not save.}

\item{posterior_pred}{Logical. If true, includes a posterior predictive distribution of the estimated IOA}

\item{coolors}{Colors for abundance.}

\item{coolors2}{Colors for catch}
}
\value{
Returns and saves a figure with the abundance and catch trajectories.
}
\description{
Function that provides a plot of the estimated abundance trends and catch history from a SIR  model including: median, 95%
credible interval, 50% credible interval, catch, and absolute indices of abundance.
}
