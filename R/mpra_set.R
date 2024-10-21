setClass("MPRASet", contains = "SummarizedExperiment")

setValidity("MPRASet", function(object) {
    ## Required information: DNA, RNA, and element ID
    msg <- validMsg(NULL, .check_assay_names(object, c("DNA", "RNA")))
    msg <- validMsg(msg, !is.null(getEid(object)))
    if("barcode" %in% names(rowData(object))) {
        if(!is.character(rowData(object)$barcode) ||
           anyDuplicated(rowData(object)$barcode))
            msg <- validMsg(msg, "`barcode` should be a character vector without duplicate values.")
    }
    if("eseq" %in% names(rowData(object)) &&
       !is.character(rowData(object)$eseq))
        msg <- validMsg(msg, "`eseq` should be a character vector")
    if(! "eid" %in% names(rowData(object)) ||
       !is.character(rowData(object)$eid))
        msg <- validMsg(msg, "`eid` should be present and be a character vector")
    if (is.null(msg)) TRUE else msg
})

MPRASet <- function(DNA = new("matrix"), RNA = new("matrix"),
                    barcode = new("character"), eid = new("character"),
                    eseq = new("character"), label = new("character"), ...) {
	#eid <- eid[order(eid)]
	DNA <- DNA[eid, ]
	RNA <- RNA[eid, ]
    assays <- SimpleList(DNA = DNA, RNA = RNA)
    rowData <- DataFrame(eid = eid)
	if (length(barcode) != 0) {
		rowData$barcode <- barcode
	}
	if (length(eseq) != 0) {
		rowData$eseq <- eseq
	}
	if (length(label) != 0) {
		label <- label[eid]
		rowData$label <- label
	}
    new("MPRASet",
        SummarizedExperiment(assays = assays, rowData = rowData, ...)
    )
}

setMethod("show", signature(object = "MPRASet"),
          function(object) {
    callNextMethod()
    .show.barcodePresence(object)
})

getDNA <- function(object, aggregate = FALSE) {
    .is_mpra_or_stop(object)
    raw <- assay(object, "DNA")
    if (aggregate) {
		eid <- getEid(object)
        by_out <- by(raw, eid, colSums, na.rm = FALSE)
        agg <- do.call("rbind", by_out)
        rownames(agg) <- names(by_out)
		#agg <- agg[order(rownames(agg)), ]
        return(agg)
    } else {
		#raw <- raw[order(rownames(raw)), ]
        return(raw)
    }
}

getRNA <- function(object, aggregate = FALSE) {
    .is_mpra_or_stop(object)
    raw <- assay(object, "RNA")
    if (aggregate) {
    	eid <- getEid(object)
        by_out <- by(raw, eid, colSums, na.rm = FALSE)
        agg <- do.call("rbind", by_out)
        rownames(agg) <- names(by_out)
        #agg <- agg[order(rownames(agg)), ]
        return(agg)
    } else {
		#raw <- raw[order(rownames(raw)), ]
        return(raw)
    }
}

getBarcode <- function(object) {
    .is_mpra_or_stop(object)
    rowData(object)$barcode
}

getEid <- function(object) {
    .is_mpra_or_stop(object)
    rowData(object)$eid
}

getEseq <- function(object) {
    .is_mpra_or_stop(object)
    rowData(object)$eseq
}

getLabel <- function(object) {
	.is_mpra_or_stop(object)
	rowData(object)$label
}