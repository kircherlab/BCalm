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
  name = c("var1", "var2", "var3", "var4"),
  Barcode = c("BC1", "BC2", "BC3", "BC4"),
  count = c(10, 20, 30, 40)
)

map_df_var_1 <- data.frame(
  ID = c("var1", "var2", "var3", "var4"),
  REF = c("var1", "var3", "var5", "var7"),
  ALT = c("var2", "var4", "var6", "var8")
)

res_var_1 <- data.frame(
  var_id = c("var1", "var2", "var3", "var4"),
  allele = c("ref", "ref", "alt", "alt"),
  Barcode = c("BC1", "BC2", "BC3", "BC4"),
  code = c(10, 20, 30, 40)
)

##testthat functions
test_that("downsample_barcodes", {

  expect_equal(as.data.frame(downsample_barcodes(empty_df)), empty_df)
  expect_equal(as.data.frame(downsample_barcodes(nan_df)), nan_df)
  expect_true(is.data.frame(downsample_barcodes(empty_df)))

})

test_that("create_var_df", {

  expect_equal(as.data.frame(create_var_df(df_var_1,map_df_var_1)), res_var_1)

})