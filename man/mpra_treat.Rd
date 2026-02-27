\name{mpra_treat}
\alias{mpra_treat}
\title{Apply logFC Threshold for MPRA Results}
\description{
  Apply a moderated t-statistic with a logFC threshold derived from
  negative controls.
}
\usage{
mpra_treat(fit, percentile = 0.975, neg_label, trend = FALSE,
  robust = FALSE, winsor.tail.p = c(0.05, 0.1))
}
\arguments{
  \item{fit}{An \code{MArrayLM} or \code{MPRASet} object.}
  \item{percentile}{Percentile of negative controls used as threshold.}
  \item{neg_label}{Label for negative control elements.}
  \item{trend}{Logical indicating if a trend is fitted.}
  \item{robust}{Logical indicating robust variance estimation.}
  \item{winsor.tail.p}{Tail probabilities for winsorization.}
}
\value{
  A data frame containing moderated t-statistics and adjusted p-values.
}
\examples{
\dontrun{
  data(BcSetExample)
  block_vector <- rep(1:2, length.out = ncol(BcSetExample))
  fit <- fit_elements(BcSetExample, normalize = TRUE, block = block_vector)
  mpra_treat(fit, percentile = 0.95, neg_label = "control_name")
}
}
