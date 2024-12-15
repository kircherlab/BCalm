##Data preparation
# Controlls without effect and tests with half positive half negative effect
dna <- as.data.frame(matrix(rnorm(20 * 10, mean = 10, sd = sqrt(0.5)), nrow = 20, ncol = 10))
rownames(dna) <- paste0("label_", sprintf("%06d", 1:20))
colnames(dna) <- paste0("sample_count_", rep(1:2, each = 5), "_bc_", rep(1:5, 2))

rna_controlls <- as.data.frame(matrix(rnorm(10 * 10, mean = 10, sd = sqrt(0.5)), nrow = 10, ncol = 10))
rna_pos <- matrix(rnorm(5 * 10, mean = 14, sd = 1), nrow = 5, ncol = 10)
rna_neg <- matrix(rnorm(5 * 10, mean = 6, sd = 1), nrow = 5, ncol = 10)
rna <- round(rbind(rna_neg, rna_pos, rna_controlls), digits = 0)

rownames(rna) <- paste0("label_", sprintf("%06d", 1:20))
colnames(rna) <- colnames(dna)

labels_vec <- c(rep("test_name", 10), rep("control_name", 10))
names(labels_vec) <- paste0("label_", sprintf("%06d", 1:20))

mpra <- MPRASet(DNA = dna, RNA = rna, eid = rownames(dna), barcode = NULL, label=labels_vec)
nr_reps <- 2
bcs <- ncol(dna)/ nr_reps
block_vector <- rep(1:nr_reps, each = bcs)

mpralm_fit <- fit_elements(object = mpra, normalize=TRUE, block = block_vector)

# forcing the output to be MPRASEt as well
mpralm_fit_endo <- fit_elements(object = mpra, normalize = TRUE, block = block_vector, endomorphic = TRUE)

# with different percentiles
result_95 <- mpra_treat(mpralm_fit, percentile = 0.95, neg_label="control_name", test_label="test_name", side="both")
result_50 <- mpra_treat(mpralm_fit, percentile = 0.50, neg_label="control_name", test_label="test_name", side="both")

# with different side options
result_right <- mpra_treat(mpralm_fit, percentile = 0.95, neg_label="control_name", test_label="test_name", side = "right")
result_left <- mpra_treat(mpralm_fit, percentile = 0.95, neg_label="control_name", test_label="test_name", side = "left")

#simple result to check output structure
result <- mpra_treat(mpralm_fit, percentile = 0.95, neg_label="control_name", test_label="test_name", side="both")
result_endo <- mpra_treat(mpralm_fit_endo, percentile = 0.95, neg_label = "control_name", test_label = "test_name", side = "both")

## testthat calls
test_that("mpra_treat", {
  expect_error(mpra_treat(mpralm_fit, percentile = 0.95, neg_label="neg_name", test_label="test_name", side="both"))
  expect_error(mpra_treat(mpralm_fit, percentile = 0.95, neg_label="control_name", test_label="name", side="both"))

  expect_true(nrow(result_95) > 0)
  expect_true(nrow(result_50) > 0)

  expect_true(all(result_right$logFC > 0))
  expect_true(all(result_left$logFC < 0))

  expect_true(nrow(result) > 0)
  expect_true("logFC" %in% colnames(result))
  expect_true("AveExpr" %in% colnames(result))

  expect_equal(result, result_endo)
  expect_error(mpra_treat(mpra))
})