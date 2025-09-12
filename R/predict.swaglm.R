#' print.swaglm
#'
#' Predict method for a \code{swaglm} object
#' @param x An object of class \code{swaglm}.
#' @param newdata A data.frame that contains the same variable as the training dataset
#' @example /inst/examples/eg_predict.swaglm.R
#' @return TODO.
#'
predict.swaglm <- function(object, newdata) {
  
  



  # # Parameters for data generation
  set.seed(12345)
  n <- 2000
  p <- 100
  # create design matrix and vector of coefficients
  Sigma <- diag(rep(1/p, p))
  X <- MASS::mvrnorm(n = n, mu = rep(0, p), Sigma = Sigma)
  beta = c(-15,-10,5,10,15, rep(0,p-5))

  # --------------------- generate from logistic regression with an intercept of one
  z <- 1 + X%*%beta
  pr <- 1/(1 + exp(-z))
  y <- as.factor(rbinom(n, 1, pr))
  y = as.numeric(y)-1

  # define swag parameters
  quantile_alpha = .15
  p_max = 20
  object = swaglm::swaglm(X=X, y = y, p_max = p_max, family = stats::binomial(),
                            alpha = quantile_alpha, verbose = TRUE, seed = 123)
  
  
  newdata = X
  
  # verify that new data has the same structure as X
  if(colnames(newdata) != colnames(X)){
    stop("'newdata' should have the same column names as the design matrix X used in the swaglm algorithm")
  }
  
  # compute summary on swaglm object to extract selected models 
  summary_swaglm = summary(object)
  
  # compute linear combination of predictors for all selected models on newdata
  
  
  
  
  
  object$family$linkinv(10)
  

  
}
