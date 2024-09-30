##Data preparation for tests

#fit_elements
dna <- as.data.frame(matrix(10, nrow = 10, ncol = 10))
rownames(dna) <- paste0("label_", sprintf("%06d", 1:10))
colnames(dna) <- paste0("sample_count_", rep(1:2, each = 5), "_bc_", rep(1:5, 2))

rna <- as.data.frame(matrix(20, nrow = 10, ncol = 10))
rownames(rna) <- paste0("label_", sprintf("%06d", 1:10))
colnames(rna) <- paste0("sample_count_", rep(1:2, each = 5), "_bc_", rep(1:5, 2))

labels_vec <- rep("test_name", 10)  
names(labels_vec) <- paste0("label_", sprintf("%06d", 1:10))  
mpra_1 <- MPRASet(DNA = dna, RNA = rna, eid = row.names(dna), barcode = NULL, label=labels_vec)

#compute_logratio
dna_log_2 <- dna
dna_log_2["label_000001", "sample_count_1_bc_1"] <- NA
rownames(dna) <- paste0("label_", sprintf("%06d", 1:10))
colnames(dna) <- paste0("sample_count_", rep(1:2, each = 5), "_bc_", rep(1:5, 2))
mpra_log_2 <- MPRASet(DNA = dna_log_2, RNA = rna, eid = row.names(dna), barcode = NULL, label=labels_vec)

res_log_1 <- as.data.frame(matrix(0, nrow = length(dna), ncol = length(dna)))
res_log_2 <- res_log_1
res_log_2[1, 1] <- NA

## testthat calls

test_that("compute_logratio", {
  expect_equal(compute_logratio(mpra_1, aggregate=none), res_log_1)
  expect_equal(compute_logratio(mpra_log_2, aggregate=none), res_log_2)
})