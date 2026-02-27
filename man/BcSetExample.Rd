\name{BcSetExample}
\alias{BcSetExample}
\alias{LabelExample}
\alias{MapExample}
\alias{MapExampleOriginalIds}
\docType{data}
\title{Example Data for BCalm}
\description{
  Example datasets for the BCalm package demonstrating typical MPRA data
  structures and preprocessing workflows.
}
\usage{
data("BcSetExample")
data("LabelExample")
data("MapExample")
data("MapExampleOriginalIds")
}
\format{
  \describe{
    \item{\code{BcSetExample}}{An MPRASet object with barcode-level DNA and RNA counts.}
    \item{\code{LabelExample}}{A data frame with labels for elements in BcSetExample.}
    \item{\code{MapExample}}{A data frame mapping barcodes to variant IDs.}
    \item{\code{MapExampleOriginalIds}}{A data frame mapping barcodes to original variant IDs.}
  }
}
\details{
  These datasets are provided for demonstration and testing of the
  barcode-counting linear model workflow.
}
\examples{
data(BcSetExample)
}
