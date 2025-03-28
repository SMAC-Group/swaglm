% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{sawglm}
\alias{sawglm}
\alias{swaglm}
\title{sawglm}
\usage{
swaglm(
  X,
  y,
  p_max = 2L,
  family = NULL,
  method = 0L,
  alpha = 0.3,
  verbose = FALSE,
  seed = 123L
)
}
\arguments{
\item{X}{A numeric matrix of predictors.}

\item{y}{A numeric vector of responses.}

\item{p_max}{maximum dimension to explore}

\item{family}{A glm family object (default is binomial).}

\item{method}{an integer scalar with value 0 for the column-pivoted QR decomposition, 1 for the unpivoted QR decomposition, 2 for the LLT Cholesky, or 3 for the LDLT Cholesky}

\item{alpha}{lower quantile to of criterion}

\item{verbose}{boolean for verbose}

\item{seed}{A \code{int} that is the random seed used when creating the set of model to explore for the next dimension}
}
\value{
A list of list of matrices. Each list correspond to the estimated coefficients from the estimated models, xxx.
}
\description{
This function runs the SWAG algorithm.
}
