#' @importFrom mpra MPRASet getRNA getDNA getBarcode getEid getEseq
#' @export MPRASet getRNA getDNA getBarcode getEid getEseq

MPRASet <- function(label = new("character"), ...) {
    # Create a new MPRASet object
    # label: A character vector with the labels of the sequences
    object <- mpra::MPRASet(...)
    .check_element_order(object)
    if (length(label) != 0) {
        eid <- getEid(object)
        label <- label[eid]
        rowData(object)$label <- label
    }
    object
}

# mpra joins measurements and weights by position; the weights follow the order of
# getDNA(aggregate = TRUE), so the rows must already be in that (locale-dependent) order.
.check_element_order <- function(object) {
    if (length(getBarcode(object)) != 0) {
        return(invisible(TRUE))
    }
    if (!identical(rownames(getDNA(object, aggregate = TRUE)), rownames(getDNA(object, aggregate = FALSE)))) {
        stop(
            "Row order of the DNA/RNA input differs from the order mpra uses for the weights ",
            "(base R sort() in the current locale). Sort the rows of DNA and RNA with ",
            "sort(rownames(...)) before creating the MPRASet."
        )
    }
    invisible(TRUE)
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
