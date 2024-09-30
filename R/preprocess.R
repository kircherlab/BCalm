downsample_barcodes <- function(df, id_column_name="name", percentile=0.95) {
	if (!id_column_name %in% names(df)) {
		warning(paste("Column", id_column_name, "does not exist in the DataFrame.
			Provide an existing column name to the variable id_column_name. Returning the original DataFrame."))
		return(df)  # Return the original DataFrame if the column does not exist
	}
	if (any(names(df) == "allele")) {

		# Calculate the 0.95th quantile of the number of barcodes across all groups
		max_bc <- df %>%
			group_by(!!sym(id_column_name), allele) %>%
			summarise(n = n(), .groups = 'drop') %>%
			summarise(max_bc = quantile(n, percentile)) %>%
			pull(max_bc)
		
		# Downsample barcodes
		df <- df %>%
			group_by(!!sym(id_column_name), allele) %>%
			mutate(row_num = sample(row_number())) %>%
			filter(row_num <= max_bc) %>%
			ungroup()
	} else {
		# Calculate the 0.95th quantile of the number of barcodes across all groups
		max_bc <- df %>%
			group_by(!!sym(id_column_name)) %>%
			summarise(n = n(), .groups = 'drop') %>%
			summarise(max_bc = quantile(n, percentile)) %>%
			pull(max_bc)
		
		# Downsample barcodes
		df <- df %>%
			group_by(!!sym(id_column_name)) %>%
			mutate(row_num = sample(row_number())) %>%
			filter(row_num <= max_bc) %>%
			ungroup()
	}

	df$row_num <- NULL

	return(df)
}

create_dna_df <- function(df, id_column_name="variant_id", allele_column_name=NULL) {
	suppressWarnings({
		if (is.null(allele_column_name) && !is.null(df$allele)) {
		allele_column_name <- "allele"
	}})
		
	df_dna <- .pivot_df(df, id_column_name, allele_column_name, "DNA")
	return(df_dna)
}

create_rna_df <- function(df, id_column_name="variant_id", allele_column_name=NULL) {
	suppressWarnings({
		if (is.null(allele_column_name) && !is.null(df$allele)) {
			allele_column_name <- "allele"
	}})

	df_rna <- .pivot_df(df, id_column_name, allele_column_name, "RNA")
	return(df_rna)
}

create_var_df <- function(df, map_df) {
	if (!all(c("ID", "REF", "ALT") %in% colnames(map_df))) {
		stop("map_df must contain columns 'ID', 'REF', and 'ALT'")
	}

	if (!all(c("name") %in% colnames(df))) {
		stop("df must contain column 'name'")
	}

	if (!any(df$name %in% map_df$REF) & !any(df$name %in% map_df$ALT)) {
		stop("No matches found between the 'name' column in 'df' and the 'REF'/'ALT' columns in 'map_df'. Please ensure that these columns have matching values.")
	}

	map_df <- map_df %>% select(ID, REF, ALT)
	
	# Merge on REF
	df_ref <- merge(df, map_df, by.x = "name", by.y = "REF", all.x = FALSE)
	df_ref$allele <- "ref"
	df_ref$ALT <- NULL

	# Merge on ALT
	df_alt <- merge(df, map_df, by.x = "name", by.y = "ALT", all.x = FALSE)
	df_alt$allele <- "alt"
	df_alt$REF <- NULL

	# Combine the results
	df_combined <- rbind(df_ref, df_alt)

	# Select and rename columns as necessary
	var_df <- df_combined %>% select(variant_id = ID, allele, Barcode, matches("count"))

	return(var_df)
}

.pivot_df <- function(df, id_column_name="variant_id", allele_column_name=NULL, type) {
	if (is.null(allele_column_name)) {
		df <- df %>% group_by(!!sym(id_column_name)) %>%
			mutate(new_idx = row_number()) %>%
			ungroup()
	} else {
		df <- df %>% group_by(!!sym(id_column_name), !!sym(allele_column_name)) %>%
			mutate(bc = row_number()) %>%
			unite("new_idx", c("bc", "allele"), sep = "_", remove = FALSE) %>%
			ungroup()
	}
	df_pivot <- df %>%
	  pivot_wider(names_from = new_idx, values_from = matches(type, ignore.case=TRUE), names_prefix = "bc", id_cols = id_column_name) %>%
	  rename_with(~ gsub(paste0("(?i)", type), "sample", .)) %>%
	  arrange(!!sym(id_column_name)) %>%
	  as.data.frame()

	row.names(df_pivot) <- df_pivot[,id_column_name]
	df_pivot[,id_column_name] <- NULL
	return(df_pivot)
}
