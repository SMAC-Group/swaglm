#' print.swaglm
#'
#' Print a \code{swaglm} object
#' @name print.swaglm
#' @method print swaglm
#' @param x An object of class \code{swaglm}.
#' @param ... Additional arguments
#' @example /inst/examples/eg_print.swaglm.R
#' @return None.
#' @export
#'
print.swaglm <- function(x, ...) {
  cat("SWAGLM results :\n")
  cat("-----------------------------------------\n")
  cat("Input matrix dimension: ", dim(x$X), "\n")
  cat("Number of explored models: ", dim(plyr::rbind.fill.matrix(x$lst_var_mat))[1], "\n")
  # cat("Number of selected models: ", dim(plyr::rbind.fill.matrix(x$lst_selected_models))[1], "\n")
  cat("Number of dimensions explored: ", length(x$lst_AIC), "\n")
}
