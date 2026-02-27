\name{create_dna_df}
\alias{create_dna_df}
\title{Create DNA Count Matrix}
\description{
  Convert a long-format DNA count data frame into a wide-format matrix
  suitable for MPRA analysis.
}
\usage{
create_dna_df(df, id_column_name = "variant_id", allele_column_name = NULL)
}
\arguments{
  \item{df}{A data frame containing DNA counts in long format.}
  \item{id_column_name}{Column name for variant identifiers.}
  \item{allele_column_name}{Column name for allele values. If NULL, uses "allele" when present.}
}
\value{
  A wide-format data frame with DNA counts per allele.
}
\examples{
df <- data.frame(
  variant_id = c("var1", "var2", "var1", "var2"),
  allele = c("ref", "ref", "alt", "alt"),
  DNA_A = c(100, 200, 150, 250)
)
create_dna_df(df)
}
