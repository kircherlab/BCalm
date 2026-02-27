\name{getLabel}
\alias{getLabel}
\title{Get Labels from an MPRASet}
\description{
  Retrieve the label vector stored in the row data of an MPRASet object.
}
\usage{
getLabel(object)
}
\arguments{
  \item{object}{An MPRASet object.}
}
\value{
  A character vector of labels.
}
\seealso{
  \code{\link{MPRASet}}
}
\examples{
\dontrun{
  data(BcSetExample)
  getLabel(BcSetExample)
}
}
