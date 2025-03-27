replicate_fit <- function(mpra, block_vector){
    nr_reps <- length(unique(block_vector))
    dna <- getDNA(mpra)
    rna <- getRNA(mpra)
    labels <- getLabel(mpra)
    bcs <- ncol(dna)/nr_reps
    result <- data.frame(oligo_names = rownames(dna))
    for (i in 1:nr_reps){
        # create new MPRASets for each replicate
        rep_dna <- dna[,(bcs*(i-1)+1):(bcs*i)]
        rep_rna <- rna[,(bcs*(i-1)+1):(bcs*i)]
        rep_block_vec <- rep(1:nr_reps, each=bcs)
        # different cases if mpra data was variant or element testing data
        if (grepl("ref|alt", colnames(dna)[1])) {
            rep_mpra <- MPRASet(DNA = rep_dna, RNA = rep_rna, eid = row.names(rep_dna), barcode = NULL)
            rep_design <- data.frame(intcpt = 1, alt = grepl("alt", colnames(rep_mpra)))
            # fit replicate with "indep_groups"
            rep_mpralm_fit <- mpralm(object = rep_mpra, design = rep_design, aggregate = "none", normalize = TRUE, model_type = "indep_groups", plot = FALSE, block = rep_block_vec)
            top_var <- topTable(rep_mpralm_fit, coef = 2, number = Inf)
        } else {
            labels_vec <- labels[rownames(rep_dna)]
            rep_mpra <- MPRASet(DNA = rep_dna, RNA = rep_rna, eid = row.names(rep_dna), barcode = NULL)
            rep_mpralm_fit <- fit_elements(object = rep_mpra, normalize=TRUE, block = rep_block_vec, endomorphic = FALSE, plot = FALSE)
            top_var <- topTable(rep_mpralm_fit, coef = 1, number = Inf)
        }
        # put into requested dataframe form
        top_var <- top_var[rownames(dna), , drop = FALSE]
        col_name <- paste0("logFC_", i)
        result[[col_name]] <- top_var$logFC
    }
    return(result)
}

calculate_corr <- function(object, type = "pearson"){
    result <- object[,-1]
    rownames(result) <- object[,1]
    corr_matrix <- cor(result, method = type)
    return(corr_matrix)
}

plot_corr <- function(object, corr_matrix){
 
    # Extract logFC values for the chosen replicates
    xvals <- object[["logFC_1"]]
    yvals <- object[["logFC_2"]]

    data <- data.frame(x = xvals,y = yvals)
    # Get correlation value from the matrix
    corr_val <- corr_matrix[1, 2]

    # Generate scatter plot
    plot <- ggplot(data, aes(x = x, y = y)) +
        geom_point(alpha = 0.5, color = "blue") +  # Scatter points
        geom_smooth(method = "lm", color = "red", se = FALSE) +  # Correlation line
        ggtitle(paste("Scatterplot of Replicate", 1, "vs Replicate", 2)) +
        xlab("Replicate 1") +
        ylab("Replicate 2") +
        theme_minimal() +
        annotate("text", x = min(xvals, na.rm = TRUE), 
                 y = max(yvals, na.rm = TRUE), 
                 label = paste("r =", round(corr_val, 3)), 
                 hjust = 0, vjust = 1, size = 5, color = "red")

    dev.off()
    return(plot)
}
