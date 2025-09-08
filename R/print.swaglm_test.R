#' print.swaglm_test
#'
#' Print a \code{swaglm_test} object
#' @param x An object of class \code{swaglm_test}.
#' @param ... Additional arguments
#' @example  /inst/examples/eg_print.swaglm_test.R
#' @export
#' 
print.swaglm_test <- function(x, ...) {
  cat("SWAGLM Test Results:\n")
  cat("----------------------\n")
  cat("p-value (Eigen):", x$p_value_eigen, "\n")
  cat("p-value (Freq):", x$p_value_freq, "\n")
}
