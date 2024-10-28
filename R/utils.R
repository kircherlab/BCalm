.is_mpra_or_stop <- function(object) {
    if (!is(object, "MPRASet"))
        stop("object is of class '", class(object), "', but needs to be of class 'MPRASet'")
}

.onLoad <- function(libname, pkgname) {
  # Override the compute_logratio function in the mpra namespace
  assignInNamespace("compute_logratio", compute_logratio, ns = "mpra")
}