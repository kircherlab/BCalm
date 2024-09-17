library(mpra)
library(testthat)
library(usethis)
library(devtools)
library(tidyr)
library(dplyr)
##Data preparation for tests

#Downsample_barcodes
empty_df <- data.frame(
  name = character(),
  allele = character(),
  Barcode = character(),
  count = numeric()
)

nan_df <- data.frame(
  name = c("A", "B", NA, "D"),
  allele = c("ref", NA, "alt", "ref"),
  Barcode = c("BC1", "BC2", NA, "BC4"),
  count = c(10, NA, 5, NA)
)

single_column_df <- data.frame(
  name = c("A", "B", "C", "D")
)

#create_var_df
df_var_1 <- data.frame(
  name = c("A", "B", "C"),
  Barcode = c("BC1", "BC2", "BC3"),
  count = c(10, 20, 30)
)

map_df_var_1 <- data.frame(
  ID = c("var1", "var2", "var3", "var4"),
  REF = c("A", "A", "B", "D"),
  ALT = c("E", "F", "G", "C")
)

res_var_1 <- data.frame(
  variant_id = c("var1", "var2", "var3", "var4"),
  allele = c("ref", "ref", "ref", "alt"),
  Barcode = c("BC1", "BC1", "BC2", "BC3"),
  count = c(10, 10, 20, 30)
)

map_df_var_2 <- data.frame(
  ID = c("var1", "var2", "var3", "var4"),
  REF = c("U", "V", "W", "X"),
  ALT = c("Y", "Z", "Q", "R")
)

res_var_2 <- data.frame(
  variant_id = character(),
  allele = character(),
  Barcode = character(),
  count = numeric()
)

#create_dna_df
df_dna_1 <- data.frame(
  variant_id = c("var1", "var2"),
  allele = c("ref", "alt"),
  DNA_A = c(100, 200),
  DNA_B = c(300, 400))

res_dna_1 <- data.frame(
  sample_A_bc1_ref = c(100, NA),
  sample_A_bc1_alt = c(NA, 200),
  sample_B_bc1_ref = c(300, NA),
  sample_B_bc1_alt = c(NA, 400),
  row.names = c("var1", "var2"))

df_dna_2 <- data.frame(
  variant_id = c("var1", "var2"),
  custom_allele = c("ref", "alt"),
  DNA_A = c(100, 200),
  DNA_B = c(300, 400))

res_dna_2 <- data.frame(
  sample_A_bc1_ref = c(100, NA),
  sample_A_bc1_alt = c(NA, 200),
  sample_B_bc1_ref = c(300, NA),
  sample_B_bc1_alt = c(NA, 400),
  row.names = c("var1", "var2"))

df_dna_3 <- data.frame(
  DNA_A = c(100, 200),
  DNA_B = c(300, 400))

df_dna_4 <- data.frame(
  variant_id = character(),
  DNA_A = numeric(),
  DNA_B = numeric())

df_dna_5 <- data.frame(
  name = c("var1", "var2"),
  allele = c("ref", "alt"),
  DNA_A = c(100, 200))

##testthat calls
test_that("downsample_barcodes", {

  expect_equal(as.data.frame(downsample_barcodes(empty_df)), empty_df)
  expect_equal(as.data.frame(downsample_barcodes(nan_df)), nan_df)
  expect_true(is.data.frame(downsample_barcodes(empty_df)))

})

test_that("create_var_df", {

  expect_equal(create_var_df(df_var_1,map_df_var_1), res_var_1)
  expect_equal(create_var_df(df_var_1,map_df_var_2), res_var_2)

})

test_that("create_dna_df", {
  
  expect_equal(create_dna_df(df_dna_1), res_dna_1)
  expect_equal(create_dna_df(df_dna_2, allele_column_name = "custom_allele"), res_dna_2)
  expect_error(create_dna_df(df_dna_3))
  expect_equal(create_dna_df(df_dna_4), data.frame())
  expect_error(create_dna_df(df_dna_5, id_column_name = "variant_id"))

})
