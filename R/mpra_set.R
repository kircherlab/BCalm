#' Create an MPRASet Object
#'
#' This function initializes a new MPRASet object from the `mpra` package, optionally adding sequence labels to the object. The `label` parameter allows you to assign specific labels to each sequence in the dataset.
#'
#' @param label A character vector containing the labels for each sequence. The length of the `label` vector must match the number of elements in the object, and it will be assigned to the `label` column in the row data of the MPRASet. Defaults to an empty character vector, meaning no labels will be applied unless specified.
#' @param ... Additional arguments passed to the `mpra::MPRASet` constructor, allowing for further customization of the object (e.g., data matrices, metadata).
#'
#' @details
#' The function creates a new `MPRASet` object, which is a container for MPRA data. If the `label` vector is provided, it will be assigned to the `label` field of the row data. The `label` will be matched to the element identifiers (EIDs) of the object.
#'
#' If no `label` is provided, the object will be returned without labels, and the user can assign them later using the `rowData` function.
#'
#' @return An `MPRASet` object with or without labels, depending on the input.
#'
#' @examples
#' \dontrun{
#' # Create an MPRASet with sequence labels (other_args need to be defined before running this example)
#' mpra_data <- MPRASet(label = c("seq1", "seq2", "seq3"), other_args)
#' }
#'
#' @importFrom SummarizedExperiment rowData
#' @export
MPRASet <- function(label=new("character"), ...) {
	# Create a new MPRASet object
	# label: A character vector with the labels of the sequences
	object <- mpra::MPRASet(...)
	if (length(label) != 0) {
		eid <- getEid(object)
		label <- label[eid]
		rowData(object)$label <- label
	}
	object
}


#' Extract Sequence Labels from an MPRA Object
#'
#' This function retrieves the sequence labels stored in the `label` field of the row data from an MPRA object. It is a simple accessor function designed to extract the label information for sequences in the dataset.
#'
#' @param object An MPRA object, typically created using the `MPRASet` function, containing sequence data and labels. The object must contain a `label` column in its row data.
#'
#' @details
#' The function checks that the input object is a valid MPRA object using `.is_mpra_or_stop`. It then extracts the sequence labels from the `label` field of the row data. This field is expected to be a vector of labels for each sequence in the object, which is typically assigned during the creation of the object.
#'
#' If the object does not contain a `label` field, this function will raise an error.
#'
#' @return A character vector of sequence labels from the row data of the MPRA object.
#'
#' @examples
#' \dontrun{
#' # Assuming `mpra_object` is an MPRA object with sequence labels
#' labels <- getLabel(mpra_object)
#' }
#'
#' @importFrom SummarizedExperiment rowData
#' @export
getLabel <- function(object) {
	.is_mpra_or_stop(object)
	rowData(object)$label
}