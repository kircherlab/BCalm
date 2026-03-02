\name{plot_groups}
\alias{plot_groups}
\title{Plot logFC Distributions by Group}
\description{
  Plot the distribution of logFC values for MPRA results, optionally
  highlighting negative control percentiles.
}
\usage{
plot_groups(object, percentile = NULL, neg_label = NULL, test_label = NULL)
}
\arguments{
  \item{object}{An MPRASet object or a data frame with \code{logFC} values.}
  \item{percentile}{Percentile for negative control thresholds.}
  \item{neg_label}{Label for negative control elements.}
  \item{test_label}{Label for test elements to subset.}
}
\value{
  A \code{ggplot} object.
}
\examples{
\dontrun{
  data(BcSetExample)
  block_vector <- rep(1:2, length.out = ncol(BcSetExample))
  fit <- fit_elements(BcSetExample, normalize = TRUE, block = block_vector)
  plot_groups(fit, percentile = 0.95, neg_label = "control_name")
}
}
