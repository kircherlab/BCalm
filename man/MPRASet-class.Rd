\name{MPRASet-class}
\Rdversion{1.1}
\docType{class}
\alias{MPRASet-class}
\alias{show,MPRASet-method}
\alias{MPRASet}
\alias{getRNA}
\alias{getDNA}
\alias{getBarcode}
\alias{getEid}
\alias{getEseq}
\title{Class \code{"MPRASet"}}
\description{
  A container for data from massively parallel reporter assays
  (MPRA). Builds on the SummarizedExperiment class.
}
\usage{
## Constructor
MPRASet(label = new("character"), ...)

## Accessors
getRNA(object, aggregate = FALSE)
getDNA(object, aggregate = FALSE)
getBarcode(object)
getEid(object)
getEseq(object)
}
\arguments{
  \item{object}{A \code{MPRASet} object.}
  \item{aggregate}{A \code{logical} indicating if data should be
    aggregated to the element level (by summing across barcodes).}
  \item{label}{A character vector of labels to store with the elements.
    When provided, the labels are aligned to the element IDs.}
  \item{...}{Further arguments passed to \code{mpra::MPRASet}.}
}
\section{Slots}{
  Slots are as described in a \code{SummarizedExperiment}.  Of
  particular interest are \code{colData} which describes the phenotype
  data.  The \code{assay} slot holds the assay data, with specific assay
  names \code{RNA} and \code{DNA} (accessed by \code{getRNA} and
  \code{getDNA}). Element and barcode data are accessible in the
  \code{rowData} slot.  We have chosen to store barcode and element as
  \code{character} to be flexible, although they are often DNA sequences
  (and could therefore be considered \code{DNAStringSet} (from package
  Biostrings)).
}
\section{Extends}{
Class \code{"\linkS4class{SummarizedExperiment}"}, directly.
}
\value{
  The constructor function \code{MPRASet} returns an object of class
  \code{"MPRASet"}.
}
\section{Accessors}{
\describe{
\item{\code{getDNA}:}{Gets the DNA channel data.}
\item{\code{getRNA}:}{Gets the RNA channel data.}
\item{\code{getBarcode}:}{Gets the barcode, if present.}
\item{\code{getEid}:}{Gets the element ID}
\item{\code{getEseq}:}{Gets the element sequence, if present.}
}
}
\seealso{
  \code{\linkS4class{SummarizedExperiment}} for the basic class that is
  used as a building block.
}
\examples{
showClass("MPRASet")
}
\keyword{classes}
