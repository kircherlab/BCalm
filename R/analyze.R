#' Plot Groups with Logratio Density
#'
#' Creates a normalized histogram with overlaid density plots of logratio values for given groups in the data.
#' If specified, it also plots percentile lines for a negative control label.
#'
#' @param mpra_fit A data frame containing at least two columns: \code{logratio} (numeric values to be plotted) and \code{label} (categorical group labels).
#' @param percentile A numeric value specifying the percentile threshold for the negative control. If \code{NULL}, defaults to 0.95 when \code{neg_label} is provided.
#' @param neg_label A character string specifying the label to use as the "negative control" group. If provided, percentile lines are added to the plot.
#' @param test_label A character string specifying the label for the test group to be compared. If provided, only entries with labels matching \code{neg_label} or \code{test_label} are plotted.
#'
#' @details
#' This function checks that the required package \code{ggplot2} is installed and loads it if available.
#' The function expects \code{mpra_fit} to have columns named \code{logratio} and \code{label}.
#' If both \code{neg_label} and \code{test_label} are provided, the plot includes only rows matching these labels.
#' Dashed vertical lines indicate the specified percentiles for the negative control group if \code{neg_label} is not \code{NULL}.
#'
#' @return A \code{ggplot} object displaying a histogram with overlaid density plot, and optional percentile lines.
#'
#' @note Requires the \code{ggplot2} package for plotting.
#'
#' @examples
#' \dontrun{
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'     mpra_data <- data.frame(logratio = rnorm(100), label = sample(c("control", "test"), 100, replace = TRUE))
#'     plot_groups(mpra_data, percentile = 0.95, neg_label = "control", test_label = "test")
#' }
#' }
#'
#' @import ggplot2
#' @export
plot_groups <- function(mpra_fit, percentile=NULL, neg_label=NULL, test_label=NULL) {
	if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("The 'ggplot2' package is required but not installed. Please install it.")
    }
	if (! "label" %in% names(mpra_fit)) {
		stop("Your mpra fit object should contain a label column.")
	}
	if (! is.null(neg_label) && ! neg_label %in% unique(mpra_fit$label)) {
		stop("The negative label you provided is not in the label column of the mpra fit object.")
	}
	if (! is.null(test_label) && ! test_label %in% unique(mpra_fit$label)) {
		stop("The test label you provided is not in the label column of the mpra fit object.")
	}
	if (is.null(percentile) && ! is.null(neg_label)) {
		percentile <- 0.95
		print("No percentile provided, using 0.95.")
	}

	df <- as.data.frame(mpra_fit)
	if (! is.null(test_label) && ! is.null(neg_label)) {
		df <- df[df$label %in% c(test_label, neg_label), ]
	}

	plot <- ggplot(df, aes(x = logratio, fill = label, y = after_stat(density))) +
		geom_histogram(alpha = 0.5, position = "identity", binwidth = 0.1) +
		geom_density(alpha = 0.2, adjust = 1) +
		theme_minimal() +
		labs(title = "Normalized Histogram of logratio Values", x = "logratio", y = "Density", color = NULL) +
		xlim(c(min(df$logratio), max(df$logratio)))

	if (!is.null(neg_label)) {

		percentile_up <- quantile(df$logratio[df$label == neg_label], percentile)
		up_label <- paste(percentile, "th percentile of negative controls")

		percentile_down <- quantile(df$logratio[df$label == neg_label], 1 - percentile)
		down_label <- paste(1 - percentile, "th percentile of negative controls")

		plot <- plot +
			geom_vline(aes(xintercept = percentile_up, color = up_label), linetype = "dashed", size = 1) +
			geom_vline(aes(xintercept = percentile_down, color = down_label), linetype = "dashed", size = 1) +
			#scale_fill_manual(values = c("reference" = "blue", neg_label = "red", "other" = "grey50")) +
			scale_color_manual(values = c(up_label = "green", down_label = "orange"), guide = guide_legend(override.aes = list(linetype = "dashed")))
	}

	return(plot)
}

#' Differential Logratio Testing for MPRA Data
#'
#' Applies a thresholding approach to identify significantly shifted logratios in test groups relative to a negative control.
#' This function allows for side-specific testing, with options to test logratios that are significantly higher, lower, or both, compared to the negative group.
#'
#' @param mpra_fit A data frame containing MPRA fit data with at least two columns: \code{logratio} (numeric values representing log ratios for each entry) and \code{label} (categorical labels for groups).
#' @param percentile Numeric, the threshold percentile to use for filtering significant shifts in logratio values. If \code{NULL}, defaults to 0.95.
#' @param neg_label A character string specifying the label to use as the "negative control" group. This parameter is required.
#' @param test_label A character string specifying the label for the test group to compare against the negative control. If \code{NULL}, all labels except the negative label are considered.
#' @param side Character string specifying the direction of the test. Accepts \code{"both"} (default), \code{"right"}, or \code{"left"} for two-sided, right-tailed, or left-tailed tests, respectively.
#'
#' @details
#' This function shifts the \code{logratio} values in \code{mpra_fit} by subtracting the mean of the negative control group (specified by \code{neg_label}) to standardize comparisons.
#' It then applies a threshold using the specified \code{percentile} on the logratios for the negative group, treating this as the cut-off for identifying shifts in the test group.
#' Side-specific testing is applied based on the \code{side} parameter, with options to filter for significant logratios in the positive (right), negative (left), or both directions.
#'
#' @return A data frame containing rows that meet the threshold criterion, with the following columns:
#' \itemize{
#'   \item \code{logFC}: Adjusted log fold change.
#'   \item \code{AveExpr}: Adjusted average expression.
#'   \item \code{name}: Sequence or gene names as extracted from row names.
#' }
#'
#' @note Requires the \code{edgeR} package for the \code{treat} function used in threshold testing.
#'
#' @examples
#' \dontrun{
#' mpra_data <- data.frame(logratio = rnorm(100), label = sample(c("control", "test"), 100, replace = TRUE))
#' significant_results <- mpra_treat(mpra_data, percentile = 0.95, neg_label = "control", test_label = "test", side = "both")
#' }
#'
#' @importFrom edgeR treat topTreat
#' @import dplyr
#' @export
mpra_treat <- function(mpra_fit, percentile=NULL, neg_label, test_label=NULL, side="both") {
	if (! "label" %in% names(mpra_fit)) {
		stop("Your mpra fit object should contain a label column.")
	}
	if (! "logratio" %in% names(mpra_fit)) {
		stop("Your mpra fit object should contain a logratio column. Did you create it using fit_elements?")
	}
	if (is.null(neg_label)) {
		stop("You need to provide the label of the negative class.")
	}
	if (! neg_label %in% unique(mpra_fit$label)) {
		stop("The negative label you provided is not in the label column of the mpra fit object.")
	}
	if (! is.null(test_label) && ! test_label %in% unique(mpra_fit$label)) {
		stop("The test label you provided is not in the label column of the mpra fit object.")
	}
	if (is.null(percentile)) {
		percentile <- 0.95
		print("No percentile provided, using 0.95.")
	}

	result <- NULL

	# We are shifting the logratios so that the mean of the negative class is 0.
	# This way we don't get into trouble with negative logFC thresholds.
	neg_mean <- mean(mpra_fit$logratio[mpra_fit$label == neg_label])
	shifted_fit <- mpra_fit
	shifted_fit$logratio <- mpra_fit$logratio - neg_mean
	shifted_fit$coefficients <- mpra_fit$coefficients - neg_mean

	if (! is.null(test_label)) {
		to_test <- shifted_fit[shifted_fit$label == test_label, ]
	} else {
		to_test <- shifted_fit[shifted_fit$label != neg_label, ]
	}

	if (side == "both" || side == "right") {
		percentile_up <- quantile(shifted_fit$logratio[shifted_fit$label == neg_label], percentile)
		tr_up <- treat(to_test, lfc=percentile_up)
		toptr_up <- topTreat(tr_up, coef = 1, number = Inf)
		# topTreat tests for logratios higher than the threshold, or lower than the negative of the threshold.
		# We instead want the upper and lower percentile, so we call topTreat twice, once for the upper and once for the lower percentile.
		# Here we filter for the upper percentile.
		toptr_up <- toptr_up %>% filter(logFC > 0)
		result <- toptr_up
	}

	if (side == "both" || side == "left") {
		percentile_down <- quantile(shifted_fit$logratio[shifted_fit$label == neg_label], 1 - percentile)
		tr_down <- treat(to_test, lfc=percentile_down)
		toptr_down <- topTreat(tr_down, coef = 1, number = Inf)
		# Here we filter for the lower percentile.
		toptr_down <- toptr_down %>% filter(logFC < 0)
		result <- rbind(result, toptr_down)
	}

	result$logFC <- result$logFC + neg_mean
	result$AveExpr <- result$AveExpr + neg_mean
	result$name <- row.names(result)

	return(result)
}