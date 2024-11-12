#' Downsample barcodes in a data frame
#'
#' This function downsamples the number of barcodes for each group in a data frame
#' to a specified percentile of the distribution of barcode counts across groups.
#'
#' @param df A data frame containing barcode data.
#' @param id_column_name String indicating the name of the column
#'   used for grouping. Default is "name".
#' @param percentile Numeric value between 0 and 1 specifying the percentile
#'   to use for downsampling. Default is 0.95.
#'
#' @return A data frame with downsampled barcodes. If the specified id_column_name
#'   doesn't exist, returns the original data frame with a warning.
#'
#' @details
#' The function first checks if an "allele" column exists in the data frame.
#' If it does, downsampling is performed within each combination of id and allele.
#' If not, downsampling is performed within each id group.
#'
#' The number of barcodes kept for each group is determined by the specified
#' percentile of the distribution of barcode counts across all groups.
#'
#' @importFrom dplyr %>% group_by summarise ungroup mutate filter pull
#' @importFrom rlang sym
#' @importFrom stats quantile
#'
#' @examples
#' df <- data.frame(
#'     name = rep(c(rep("A_seq", 2), rep("B_seq", 2), rep("C_seq", 6)), each = 2),
#'     barcode = paste0("barcode_", 1:20),
#'     allele = rep(c("ref", "alt"), 10),
#'     count = runif(20)
#' )
#'
#' downsampled_df <- downsample_barcodes(df, id_column_name = "name", percentile = 0.8)
#' @export
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


#' Generate DNA count specific dataframe
#'
#' @param df Data frame with 'name' column and the barcode and count information for each sequence
#' @param id_column_name String indicating the name of the ID column. Default is "variant_id"
#' @param allele_column_name String indicating the name of the allele column. Default is NULL
#'
#' @return Data frame with variant_id, allele, Barcode, and count columns
#'
#' @export
create_dna_df <- function(df, id_column_name="variant_id", allele_column_name=NULL) {
	suppressWarnings({
		if (is.null(allele_column_name) && !is.null(df$allele)) {
		allele_column_name <- "allele"
	}})

	df_dna <- .pivot_df(df, id_column_name, allele_column_name, "DNA")
	return(df_dna)
}


#' Generate RNA count specific dataframe
#'
#' @param df Data frame with 'name' column and the barcode and count information for each sequence
#' @param id_column_name String indicating the name of the ID column. Default is "variant_id"
#' @param allele_column_name String indicating the name of the allele column. Default is NULL
#'
#' @return Data frame with variant_id, allele, Barcode, and count columns
#'
#' @export
create_rna_df <- function(df, id_column_name="variant_id", allele_column_name=NULL) {
	suppressWarnings({
		if (is.null(allele_column_name) && !is.null(df$allele)) {
			allele_column_name <- "allele"
	}})

	df_rna <- .pivot_df(df, id_column_name, allele_column_name, "RNA")
	return(df_rna)
}


#' Match count information and variant alleles to create a variant data frame
#'
#' @param df Data frame with 'name' column and the barcode and count information for each sequence
#' @param map_df Data frame with 'ID', 'REF', and 'ALT' columns
#'
#' @return Data frame with variant_id, allele, Barcode, and count columns
#'
#' @importFrom dplyr %>% select matches
#' @examples
#' df <- data.frame(name = c("ref1", "ref2", "alt1", "alt2"), Barcode = 1:4, dna_count= 2:5, rna_count = 10:13)
#' map_df <- data.frame(ID = c("rs1", "rs2"), REF = c("ref1", "ref2"), ALT = c("alt1", "alt2"))
#' result <- create_var_df(df, map_df)
#' head(result)
#' @export
create_var_df <- function(df, map_df) {
	if (!all(c("ID", "REF", "ALT") %in% colnames(map_df))) {
		stop("map_df must contain columns 'ID', 'REF', and 'ALT'")
	}

	if (!all(c("name") %in% colnames(df))) {
		stop("df must contain column 'name'")
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

#' Pivot a data frame
#'
#' This function pivots a data frame, creating a wider format based on specified columns.
#'
#' @param df A data frame to be pivoted
#' @param id_column_name String indicating the name of the ID column. Default is "variant_id"
#' @param allele_column_name String indicating the name of the allele column. Default is NULL
#' @param type String indicating the type of values to pivot
#'
#' @return A pivoted data frame with rows named by the ID column and columns prefixed with "bc"
#'
#' @importFrom dplyr %>% group_by mutate ungroup rename_with arrange
#' @importFrom tidyr pivot_wider unite
#' @importFrom rlang sym
#'
#' @examples
#' df <- data.frame(
#'     variant_id = c("rs1", "rs1", "rs2", "rs2"),
#'     allele = c("ref", "alt", "ref", "alt"),
#'     count = 1:4
#' )
#' pivoted_df <- .pivot_df(df,
#'     id_column_name = "variant_id",
#'     allele_column_name = "allele", type = "count"
#' )
#'
#' @keywords internal
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
