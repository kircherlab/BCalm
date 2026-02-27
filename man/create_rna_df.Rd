\name{create_rna_df}
\alias{create_rna_df}
\title{Create RNA Count Matrix}
\description{
  Convert a long-format RNA count data frame into a wide-format matrix
  suitable for MPRA analysis.
}
\usage{
create_rna_df(df, id_column_name = "variant_id", allele_column_name = NULL)
}
\arguments{
  \item{df}{A data frame containing RNA counts in long format.}
  \item{id_column_name}{Column name for variant identifiers.}
  \item{allele_column_name}{Column name for allele values. If NULL, uses "allele" when present.}
}
\value{
  A wide-format data frame with RNA counts per allele.
}
\examples{
df <- data.frame(
  variant_id = c("var1", "var2", "var1", "var2"),
  allele = c("ref", "ref", "alt", "alt"),
  RNA_A = c(500, 600, 450, 700)
)
create_rna_df(df)
}
