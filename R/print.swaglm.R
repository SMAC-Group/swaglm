#' print.swaglm
#'
#' Print a \code{swaglm} object
#'
#' @param x An object of class \code{swaglm}.
#' @export
#' @examples
#' n <- 2000
#' p <- 50
#' # create design matrix and vector of coefficients
#' Sigma <- diag(rep(1/p, p))
#' X <- MASS::mvrnorm(n = n, mu = rep(0, p), Sigma = Sigma)
#' beta = c(-15,-10,5,10,15, rep(0,p-5))
#' z <- 1 + X%*%beta
#' pr <- 1/(1 + exp(-z))
#' y <- as.factor(rbinom(n, 1, pr))
#' y = as.numeric(y)-1
#' quantile_alpha = .15
#' p_max = 20
#' swag_obj = swaglm::swaglm(X=X, y = y, p_max = p_max, family = stats::binomial(),
#'                           alpha = quantile_alpha, verbose = TRUE, seed = 123)
#' print(swag_obj)
#' 
#' 
print.swaglm <- function(x, ...) {
  cat("SWAGLM results :\n")
  cat("-----------------------------------------\n")
  cat("Input matrix dimension: ", dim(x$X), "\n")
  cat("Number of explored models: ", dim(plyr::rbind.fill.matrix(x$lst_var_mat))[1], "\n")
  cat("Number of selected models: ", dim(plyr::rbind.fill.matrix(x$lst_selected_models))[1], "\n")
  cat("Number of dimensions explored: ", length(x$lst_AIC), "\n")
}
