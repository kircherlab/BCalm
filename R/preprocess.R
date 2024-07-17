downsample_barcodes <- function(df) {
  	df <- df %>% group_by(variant_id, allele) %>%
    	mutate(rank = sample(row_number())) %>%
    	filter(rank <= quantile(rank, 0.95)) %>%
    	ungroup()

	df$rank <- NULL

	return(df)
}

create_dna_df <- function(df, id_column_name="variant_id", allele_column_name="allele") {
	df_dna <- .pivot_df(df, id_column_name, allele_column_name, "DNA")
	return(df_dna)
}

create_rna_df <- function(df, id_column_name="variant_id", allele_column_name="allele") {
	df_rna <- .pivot_df(df, id_column_name, allele_column_name, "RNA")
	return(df_rna)
}

.pivot_df <- function(df, id_column_name="variant_id", allele_column_name="allele", type) {
	df <- df %>% group_by(!!sym(id_column_name), !!sym(allele_column_name)) %>%
		mutate(bc = row_number()) %>%
		unite("new_idx", c("bc", "allele"), sep = "_", remove = FALSE) %>%
		ungroup()

	df_pivot <- df %>%
	  pivot_wider(names_from = new_idx, values_from = matches(type, ignore.case=TRUE), names_prefix = "bc", id_cols = id_column_name) %>%
	  rename_with(~ gsub(paste0("(?i)", type), "sample", .)) %>%
	  arrange(!!sym(id_column_name)) %>%
	  as.data.frame()

	row.names(df_pivot) <- df_pivot[,id_column_name]
	df_pivot[,id_column_name] <- NULL
	return(df_pivot)
}
