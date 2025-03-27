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

plot_corr<- function(object, corr_matrix, output_file, axis_limits = NULL) {
    
    #Scan for all logfc columns
    logFC_cols <- grep("logFC_", colnames(object), value = TRUE)
    
    # create nececsary plot combinations
    pairs <- combn(logFC_cols, 2, simplify = FALSE)
    
    # min/max for axis limit
    if (is.null(axis_limits)) {
        all_vals <- unlist(object[logFC_cols], use.names = FALSE)
        axis_limits <- c(min(all_vals, na.rm = TRUE), max(all_vals, na.rm = TRUE))
    }
    
    # Dataset for plot
    plot_data <- data.frame()
    
    for (pair in pairs) {
        rep_x <- pair[1]
        rep_y <- pair[2]
        
        xvals <- object[[rep_x]]
        yvals <- object[[rep_y]]
        
        rep_num_x <- gsub("logFC_", "", rep_x)
        rep_num_y <- gsub("logFC_", "", rep_y)
        
        # Get correlations
        corr_val <- round(corr_matrix[as.numeric(rep_num_x), as.numeric(rep_num_y)], 3)
        
        # New data
        plot_data <- rbind(plot_data, data.frame(
            x = xvals, 
            y = yvals, 
            replicate_pair = paste("Rep", rep_num_x, "vs Rep", rep_num_y, "\nr =", corr_val)
        ))
    }
    
    #Scatterplot
    plot <- ggplot(plot_data, aes(x = x, y = y)) +
        geom_point(alpha = 0.5, color = "black") +
        geom_smooth(method = "lm", color = "blue", se = FALSE) +  
        geom_abline(slope = 1, intercept = 0, color = "black", linetype = "dashed") + 
        theme_minimal() +
        theme(
            panel.border = element_rect(color = "black", fill = NA, size = 1),  
            axis.line = element_line(color = "black", size = 1),
            axis.text = element_text(size = 10),
            axis.title = element_text(size = 12)
        ) +
        coord_fixed(ratio = 1) +  
        scale_x_continuous(limits = axis_limits) +  
        scale_y_continuous(limits = axis_limits) +
        facet_wrap(~replicate_pair, scales = "fixed")  # All in one figure
        
    # Save
    png(output_file)
    print(plot)  
    dev.off()
    
    return(plot)
}
