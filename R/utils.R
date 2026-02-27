# Declare global variables to avoid R CMD check NOTEs
utils::globalVariables(c(
    # Functions from mpra package
    "mpralm", "getDNA", "getRNA", "getEid", "getBarcode", "getEseq",
    "normalize_counts", "get_precision_weights", "compute_logratio",
    "MPRASet", "getLabel",
    # Functions from other packages
    "squeezeVar", "rowData",
    # ggplot2 functions (in Suggests)
    "ggplot", "aes", "geom_histogram", "geom_density", "theme_minimal",
    "labs", "xlim", "geom_vline", "scale_color_manual", "guide_legend",
    "after_stat",
    # NSE column names used in dplyr and ggplot2 operations
    "allele", "row_num", "n", "ID", "REF", "ALT", "max_bc", "label", "logFC",
    "Barcode", "new_idx"
))

.is_mpra_or_stop <- function(object) {
    if (!is(object, "MPRASet")) {
        stop("object is of class '", class(object), "', but needs to be of class 'MPRASet'")
    }
}

.onLoad <- function(libname, pkgname) {
    # Override the compute_logratio function in the mpra namespace
    ns <- base::getNamespace("mpra")
    base::unlockBinding("compute_logratio", ns)
    utils::assignInNamespace("compute_logratio", compute_logratio, ns = "mpra")
    base::lockBinding("compute_logratio", ns)
}
