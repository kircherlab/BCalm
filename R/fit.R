fit_elements <- function(object, normalize=TRUE, block = NULL, endomorphic=FALSE, normalizeSize=1e9, ...) {
	design <- data.frame(rep(1, ncol(object)))
	if (normalize) {
		if ("normalizeSize" %in% names(formals(normalize_counts))) {
			object <- normalize_counts(object, normalizeSize=normalizeSize, block=block)
		} else {
			object <- normalize_counts(object, block=block)
		}
	}
	if ("endomorphic" %in% names(formals(mpralm))) {
		mpralm_fit <- mpralm(object = object, design = design, aggregate = "none", 
						normalize = F, model_type = "indep_groups", 
						block = block, endomorphic = endomorphic, normalizeSize = normalizeSize, ...)
	} else {
		mpralm_fit <- mpralm(object = object, design = design, aggregate = "none", 
						normalize = F, model_type = "indep_groups", 
						block = block, ...)
	}
	if (! endomorphic) {
		mpralm_fit$label <- getLabel(object)
		mpralm_fit$logFC <- mpralm_fit$coefficients
	}
	return(mpralm_fit)
}

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
