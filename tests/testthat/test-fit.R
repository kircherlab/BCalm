##Data preparation for tests

# compute_logratio
dna <- as.data.frame(matrix(10, nrow = 10, ncol = 10))
rownames(dna) <- paste0("name_", sprintf("%06d", 1:10))
colnames(dna) <- paste0("sample_count_", rep(1:2, each = 5), "_bc_", rep(1:5, 2))

rna <- as.data.frame(matrix(10, nrow = 10, ncol = 10))
rownames(rna) <- paste0("name_", sprintf("%06d", 1:10))
colnames(rna) <- paste0("sample_count_", rep(1:2, each = 5), "_bc_", rep(1:5, 2))

# first set with no effect
labels_vec <- c(rep("test_label_1", 5), rep("test_label_2",5))
names(labels_vec) <- rownames(dna)  
mpra_log_1 <- MPRASet(DNA = dna, RNA = dna, eid = row.names(dna), barcode = NULL, label=labels_vec)

# second set with missing value in first cell
dna_log_2 <- dna
dna_log_2["name_000001", "sample_count_1_bc_1"] <- NA
mpra_log_2 <- MPRASet(DNA = dna_log_2, RNA = dna, eid = row.names(dna), barcode = NULL, label=labels_vec)

# to test outcome type endometric == TRUE
fit_MPRA_Set <- fit_elements(object = mpra_log_2, normalize=TRUE, block = block_vector, endomorphic = TRUE)
fit_MArrayLM <- fit_elements(object = mpra_log_2, normalize=TRUE, block = block_vector)

#expected results
res_log_1 <- as.data.frame(matrix(0, nrow = length(dna), ncol = length(dna)))
rownames(res_log_1) <- paste0("name_", sprintf("%06d", 1:10))
colnames(res_log_1) <- paste0("sample_count_", rep(1:2, each = 5), "_bc_", rep(1:5, 2))
res_log_2 <- res_log_1
res_log_2[1, 1] <- NA

## testthat calls

test_that("compute_logratio", {
  expect_equal(compute_logratio(mpra_log_1, aggregate="none"), res_log_1)
  expect_equal(compute_logratio(mpra_log_2, aggregate="none"), res_log_2)
})

test_that("fit_elements", {
  expect_equal(class(fit_MPRA_Set)[1], "MPRASet")
  expect_equal(class(fit_MArrayLM)[1], "MArrayLM")
})