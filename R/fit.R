#' Fit Elements Using MPRA Linear Modeling
#'
#' This function fits an MPRA (Massively Parallel Reporter Assay) object using a linear modeling approach, optionally normalizing the data. It returns the model fit, including calculated log ratios.
#'
#' @param object An MPRA data object or expression matrix where each column represents a sample, and each row represents an element or sequence to be analyzed.
#' @param normalize Logical; if \code{TRUE} (default), normalizes the data before fitting the model.
#' @param block Optional factor specifying blocking structure, if applicable. Use when observations are not independent (e.g., paired or grouped samples).
#' @param ... Additional arguments passed to \code{mpralm}, allowing customization of the model fitting process.
#'
#' @details
#' This function uses the \code{mpralm} function to perform linear modeling on each element or sequence in the MPRA data, using an intercept-only design matrix. The model type is set to "indep_groups" by default, assuming independent groups.
#' If \code{normalize = TRUE}, normalization is applied to the data before fitting.
#' The resulting fit object includes log ratios as \code{logratio}, calculated from model coefficients, and labels retrieved from \code{getLabel(object)}.
#'
#' @return A list-like object (of class \code{mpralm_fit}) containing the model fit for each element, with the following components:
#' \itemize{
#'   \item \code{coefficients}: The model coefficients representing log ratios.
#'   \item \code{logratio}: Calculated log ratios for each element.
#'   \item \code{label}: Labels extracted from the \code{object}, indicating the element or sequence grouping.
#' }
#'
#' @examples
#' \dontrun{
#' # Fit MPRA elements with default normalization
#' fit <- fit_elements(mpra_object, normalize = TRUE)
#' }
#'
#' @importFrom mpralm mpralm
#' @export
fit_elements <- function(object, normalize=TRUE, block = NULL, ...) {
	design <- data.frame(rep(1, ncol(object)))
	mpralm_fit <- mpralm(object = object, design = design, aggregate = "none",
						normalize = normalize, model_type = "indep_groups",
						block = block, ...)
	mpralm_fit$logratio <- mpralm_fit$coefficients
	mpralm_fit$label <- getLabel(object)
	return(mpralm_fit)
}

#' Compute Log Ratios for MPRA Data
#'
#' This function computes log ratios between RNA and DNA counts in an MPRA (Massively Parallel Reporter Assay) object. The log ratio calculation can be aggregated by "mean" or "sum," or performed without aggregation.
#'
#' @param object An MPRA data object containing RNA and DNA counts to be used for log ratio computation.
#' @param aggregate Character string specifying the aggregation method to be used. Options are \code{"mean"} for element-wise mean log ratios, \code{"sum"} for summed RNA and DNA counts before computing log ratios, or \code{"none"} to calculate log ratios directly without aggregation. Default is \code{"mean"}.
#'
#' @details
#' The function calculates the log ratio as \eqn{log2(RNA + 1) - log2(DNA + 1)}. Aggregation settings allow flexibility:
#' \itemize{
#'   \item \code{"sum"}: Sum DNA and RNA counts for each element, then compute the log ratio.
#'   \item \code{"none"}: Compute log ratios directly without aggregation.
#'   \item \code{"mean"}: Calculate log ratios without aggregation, then average values for each element.
#' }
#'
#' If \code{aggregate = "mean"}, the function first computes log ratios for each replicate, then averages them by element identifier.
#'
#' @return A matrix of log ratios. If \code{aggregate = "mean"}, the matrix will contain averaged log ratios for each element.
#'
#' @examples
#' \dontrun{
#' # Compute log ratios with mean aggregation
#' log_ratios <- compute_logratio(mpra_object, aggregate = "mean")
#' }
#'
#' @export
compute_logratio <- function(object, aggregate = c("mean", "sum", "none")) {
    .is_mpra_or_stop(object)

    aggregate <- match.arg(aggregate)

    if (aggregate=="sum")  {
        dna <- getDNA(object, aggregate = TRUE)
        rna <- getRNA(object, aggregate = TRUE)
        logr <- log2(rna + 1) - log2(dna + 1)
    } else if (aggregate == "none") {
        dna <- getDNA(object, aggregate = FALSE)
        rna <- getRNA(object, aggregate = FALSE)
        logr <- log2(rna + 1) - log2(dna + 1)
    } else if (aggregate=="mean") {
        dna <- getDNA(object, aggregate = FALSE)
        rna <- getRNA(object, aggregate = FALSE)
        eid <- getEid(object)
        logr <- log2(rna + 1) - log2(dna + 1)

        by_out <- by(logr, eid, colMeans, na.rm = TRUE)
        logr <- do.call("rbind", by_out)
        rownames(logr) <- names(by_out)
    }
    return(logr)
}
