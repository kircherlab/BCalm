plot_groups <- function(mpra_fit, percentile, neg_label=NULL, test_label=NULL) {
	if (! "label" %in% names(mpra_fit)) {
		stop("Your mpra fit object should contain a label column.")
	}
	if (! is.null(neg_label) && ! neg_label %in% unique(mpra_fit$label)) {
		stop("The negative label you provided is not in the label column of the mpra fit object.")
	}
	if (! is.null(test_label) && ! test_label %in% unique(mpra_fit$label)) {
		stop("The test label you provided is not in the label column of the mpra fit object.")
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

mpra_treat <- function(mpra_fit, percentile, neg_label, test_label, side="both") {
	if (! "label" %in% names(mpra_fit)) {
		stop("Your mpra fit object should contain a label column.")
	}
	if (! "logratio" %in% names(mpra_fit)) {
		stop("Your mpra fit object should contain a logratio column. Did you create it using fit_elements?")
	}
	if (is.null(neg_label)) {
		stop("You need to provide the label of the negative class.")
	}
	if (is.null(test_label)) {
		stop("You need to provide the label of the test class.")
	}
	if (! neg_label %in% unique(mpra_fit$label)) {
		stop("The negative label you provided is not in the label column of the mpra fit object.")
	}
	if (! test_label %in% unique(mpra_fit$label)) {
		stop("The test label you provided is not in the label column of the mpra fit object.")
	}

	result <- NULL
	if (side == "both" || side == "right") {
		percentile_up <- quantile(mpra_fit$logratio[mpra_fit$label == neg_label], percentile)
		tr_up <- treat(mpra_fit, lfc=percentile_up)
		toptr_up <- topTreat(tr_up, coef = 1, number = Inf)
		# topTreat tests for logratios higher than the threshold, or lower than the negative of the threshold.
		# We instead want the upper and lower percentile, so we call topTreat twice, once for the upper and once for the lower percentile.
		# Here we filter for the upper percentile.
		toptr_up <- toptr_up %>% filter(logFC > 0)
		result <- toptr_up
	}

	if (side == "both" || side == "left") {
		percentile_down <- quantile(mpra_fit$logratio[mpra_fit$label == neg_label], 1 - percentile)
		tr_down <- treat(mpra_fit, lfc=percentile_down)
		toptr_down <- topTreat(tr_down, coef = 1, number = Inf)
		# Here we filter for the lower percentile.
		toptr_down <- toptr_down %>% filter(logFC < 0)
		result <- rbind(result, toptr_down)
	}

	result$variant_id <- row.names(result)
	return(result)
}