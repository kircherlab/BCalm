\name{downsample_barcodes}
\alias{downsample_barcodes}
\title{Downsample Barcodes}
\description{
  Downsample barcodes to a specified percentile of their distribution to
  balance barcode counts across groups.
}
\usage{
downsample_barcodes(df, id_column_name = "name", percentile = 0.95)
}
\arguments{
  \item{df}{A data frame containing barcode information.}
  \item{id_column_name}{Column name for the element or sequence ID.}
  \item{percentile}{Percentile threshold used to downsample barcodes.}
}
\value{
  A data frame with downsampled barcodes.
}
\examples{
df <- data.frame(
  name = c("A", "B", "C", "A", "B", "C"),
  allele = c("ref", "ref", "alt", "ref", "alt", "alt"),
  Barcode = c("BC1", "BC2", "BC3", "BC1", "BC2", "BC3"),
  count = c(10, 50, 30, 20, 40, 35)
)
downsample_barcodes(df)
}
