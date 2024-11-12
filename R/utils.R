#' Validate MPRA Object Class
#'
#' This function checks whether the provided object is of class `MPRASet`. If the object does not belong to this class, an error is raised.
#'
#' @param object An R object to be checked. It is expected to be an instance of the `MPRASet` class.
#'
#' @details
#' The function performs a class check on the input `object`. If the object is not of class `MPRASet`, it stops the execution and raises an informative error message, specifying the actual class of the object.
#' This ensures that subsequent operations requiring an `MPRASet` object are only executed on valid inputs.
#'
#' @return This function does not return a value. If the class of `object` is not `MPRASet`, an error is thrown.
#'
#' @examples
#' \dontrun{
#' # Assuming `mpra_object` is an MPRA object
#' .is_mpra_or_stop(mpra_object) # This will pass if mpra_object is of class MPRASet
#' }
#'
#' @keywords internal
.is_mpra_or_stop <- function(object) {
    if (!is(object, "MPRASet"))
        stop("object is of class '", class(object), "', but needs to be of class 'MPRASet'")
}


#' Initialize Package Environment on Load
#'
#' This function is called automatically when the package is loaded. It overrides the `compute_logratio` function in the `mpra` namespace with a custom implementation.
#'
#' @param libname A character string representing the library name. This is passed by the R environment during package loading.
#' @param pkgname A character string representing the package name. This is passed by the R environment during package loading.
#'
#' @details
#' The function is invoked when the package is loaded using `library()`. It assigns a custom implementation of the `compute_logratio` function to the `mpra` namespace, ensuring that the modified version of `compute_logratio` is used in the context of the loaded package.
#'
#' This is useful when you want to override functions in a package namespace at the time of loading, without modifying the original package source code.
#'
#' @return This function does not return a value. It is executed for side effects during the loading of the package.
#'
#' @examples
#' \dontrun{
#' # This function is executed automatically when the package is loaded
#' # No direct call is needed.
#' }
#'
#' @keywords internal
.onLoad <- function(libname, pkgname) {
  # Override the compute_logratio function in the mpra namespace
  assignInNamespace("compute_logratio", compute_logratio, ns = "mpra")
}