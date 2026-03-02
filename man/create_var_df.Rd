\name{create_var_df}
\alias{create_var_df}
\title{Create Variant Mapping Data}
\description{
  Map barcode-level data to variant IDs using a mapping data frame that
  contains reference and alternate allele sequences.
}
\usage{
create_var_df(df, map_df)
}
\arguments{
  \item{df}{A data frame containing barcode data with a "name" column.}
  \item{map_df}{A data frame with columns "ID", "REF", and "ALT".}
}
\value{
  A data frame with variant IDs and allele assignments.
}
\examples{
df <- data.frame(
  name = c("A", "B", "C"),
  Barcode = c("BC1", "BC2", "BC3"),
  count = c(10, 20, 30)
)
map_df <- data.frame(
  ID = c("var1", "var2", "var3"),
  REF = c("A", "B", "D"),
  ALT = c("E", "F", "C")
)
create_var_df(df, map_df)
}
