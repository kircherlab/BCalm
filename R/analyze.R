plot_groups <- function(object, percentile=NULL, neg_label=NULL, test_label=NULL) {
	if (is(object, "MPRASet")) {
		object <- rowData(object)
	} else {
		object <- as.data.frame(object)
	}
	if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("The 'ggplot2' package is required but not installed. Please install it.")
    }
	if (is.null(object$logFC)) {
		stop("Your MPRASet does not contain logFC values. Please run mpralm or fit_elements first.")
	}
	if (is.null(object$label) && (! is.null(neg_label) || ! is.null(test_label))) {
		stop("Your MPRASet should contain labels.")
	}
	if (! is.null(neg_label) && ! neg_label %in% unique(object$label)) {
		stop("The negative label you provided is not in the label column of the mpra fit object.")
	}
	if (! is.null(test_label) && ! test_label %in% unique(object$label)) {
		stop("The test label you provided is not in the label column of the mpra fit object.")
	}
	if (is.null(percentile) && ! is.null(neg_label)) {
		percentile <- 0.95
		print("No percentile provided, using 0.95.")
	}

	if (! is.null(test_label) && ! is.null(neg_label)) {
		object <- object[object$label %in% c(test_label, neg_label), ]
	} 

	plot <- ggplot(object, aes(x = logFC, fill = label, y = after_stat(density))) +
		geom_histogram(alpha = 0.5, position = "identity", binwidth = 0.1) +
		geom_density(alpha = 0.2, adjust = 1) + 
		theme_minimal() +
		labs(title = "Normalized Histogram of logratio Values", x = "Logratio", y = "Density", color = NULL) + 
		xlim(c(min(object$logFC), max(object$logFC))) 

	if (!is.null(neg_label)) {

		percentile_up <- quantile(object$logFC[object$label == neg_label], percentile)
		up_label <- paste(percentile, "th percentile of negative controls", sep="")

		percentile_down <- quantile(object$logFC[object$label == neg_label], 1 - percentile)
		down_label <- paste(1 - percentile, "th percentile of negative controls", sep="")

		plot <- plot + 
			geom_vline(aes(xintercept = percentile_up, color = up_label), linetype = "dashed", size = 1) +
			geom_vline(aes(xintercept = percentile_down, color = down_label), linetype = "dashed", size = 1) +
			scale_color_manual(values = setNames(c("green", "orange"), c(up_label, down_label)), 
				guide = guide_legend(override.aes = list(linetype = "dashed"))) 
	}

	return(plot)
}

mpra_treat <- function(mpra, percentile=NULL, neg_label, test_label=NULL, side="both") {
	if (is(mpra, "MPRASet")) {
		object <- attr(mpra, "MArrayLM")
		object$logFC <- rowData(mpra)$logFC
		object$label <- getLabel(mpra)
	} else {
		object <- mpra
	}
	if (is.null(object$label)) {
		stop("Your mpra fit object should contain a label column.")
	}
	if (is.null(object$logFC)) {
		stop("Your MPRASet does not contain logFC values. Please run mpralm or fit_elements first.")
	}
	if (is.null(neg_label)) {
		stop("You need to provide the label of the negative class.")
	}
	if (! neg_label %in% unique(object$label)) {
		stop("The negative label you provided is not in the label column of the mpra fit object.")
	}
	if (! is.null(test_label) && ! test_label %in% unique(object$label)) {
		stop("The test label you provided is not in the label column of the mpra fit object.")
	}
	if (is.null(percentile)) {
		percentile <- 0.95
		print("No percentile provided, using 0.95.")
	}
	if (side != "both" && side != "left" && side != "right") {
		stop("The side parameter should be either 'both', 'left' or 'right'.")
	}

	result <- NULL

	# We are shifting the logratios so that the mean of the negative class is 0.
	# This way we don't get into trouble with negative logFC thresholds.
	neg_mean <- mean(object$logFC[object$label == neg_label])
	shifted_fit <- object
	shifted_fit$logFC <- object$logFC - neg_mean
	shifted_fit$logFC <- object$logFC - neg_mean

	if (! is.null(test_label)) {
		to_test <- shifted_fit[shifted_fit$label == test_label, ]
	} else {
		to_test <- shifted_fit[shifted_fit$label != neg_label, ]
	}

	if (side == "both" || side == "right") {
		percentile_up <- quantile(shifted_fit$logFC[shifted_fit$label == neg_label], percentile)
		tr_up <- treat(to_test, lfc=percentile_up)
		toptr_up <- topTreat(tr_up, coef = 1, number = Inf)
		# topTreat tests for logratios higher than the threshold, or lower than the negative of the threshold.
		# We instead want the upper and lower percentile, so we call topTreat twice, once for the upper and once for the lower percentile.
		# Here we filter for the upper percentile.
		toptr_up <- toptr_up %>% filter(logFC > 0)
		result <- toptr_up
	}

	if (side == "both" || side == "left") {
		percentile_down <- quantile(shifted_fit$logFC[shifted_fit$label == neg_label], 1 - percentile)
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