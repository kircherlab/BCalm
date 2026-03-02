\name{fit_elements}
\alias{fit_elements}
\title{Fit Elements with Linear Models}
\description{
  Fit linear models to MPRA data using \code{mpralm}, with optional
  endomorphic output.
}
\usage{
fit_elements(object, normalize = TRUE, block = NULL, endomorphic = FALSE,
  normalizeSize = 1e9, plot = FALSE, ...)
}
\arguments{
  \item{object}{An MPRASet object containing DNA and RNA counts.}
  \item{normalize}{Logical indicating whether to normalize counts.}
  \item{block}{Optional vector or factor specifying replicate blocks.}
  \item{endomorphic}{Logical indicating whether to return an MPRASet result.}
  \item{normalizeSize}{Library size for normalization.}
  \item{plot}{Logical indicating whether to plot the voom precision weights.}
  \item{...}{Additional arguments passed to \code{mpralm}.}
}
\value{
  An \code{MArrayLM} object, or an \code{MPRASet} if \code{endomorphic} is TRUE.
}
\examples{
\dontrun{
  data(BcSetExample)
  block_vector <- rep(1:2, length.out = ncol(BcSetExample))
  fit <- fit_elements(BcSetExample, normalize = TRUE, block = block_vector)
}
}
