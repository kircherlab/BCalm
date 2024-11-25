#' @importFrom mpra MPRASet getRNA getDNA getBarcode getEid getEseq
#' @export MPRASet getRNA getDNA getBarcode getEid getEseq

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

getLabel <- function(object) {
	.is_mpra_or_stop(object)
	rowData(object)$label
}

setLabel <- function(object, label) {
	.is_mpra_or_stop(object)
	eid <- getEid(object)
	label <- label[eid]
	rowData(object)$label <- label
	object
}