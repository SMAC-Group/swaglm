#' predict.swaglm
#'
#' Predict method for a \code{swaglm} object
#'
#' Computes predictions from a fitted \code{swaglm} object on new data. 
#' The function returns the linear predictors (eta) for each selected model 
#' and the corresponding predicted responses using the model's inverse link function.
#'
#' @param object An object of class \code{swaglm} returned by \code{\link[swaglm]{swaglm}}.
#' @param newdata A data.frame or matrix containing the same predictors (columns) as used in the original training dataset.
#' @param ... Further arguments passed to or from other methods.
#'
#' @return A list with two elements:
#' \describe{
#'   \item{mat_eta_prediction}{A numeric matrix of linear predictors (\eqn{\eta = X \beta}) with rows corresponding to observations in \code{newdata} and columns to the selected models.}
#'   \item{mat_reponse_prediction}{A numeric matrix of predicted responses (after applying the inverse link function) with the same dimensions as \code{mat_eta_prediction}.}
#' }
#'
#' @details
#' The function performs the following steps:
#' \itemize{
#'   \item Checks that \code{newdata} has the same number of columns as the design matrix used in the fitted \code{swaglm} object.
#'   \item Computes the linear predictors (\eqn{\eta = X \beta}) for each selected model in the swaglm object.
#'   \item Applies the inverse of the model's link function to compute predicted responses.
#' }
#'
#' @example /inst/examples/eg_predict.swaglm.R
#' @export
#'
predict.swaglm <- function(object, newdata, ...) {
  
  # # 
  # # 
  # # 
  # # 
  # # # Parameters for data generation
  # set.seed(12345)
  # n <- 2000
  # p <- 100
  # # create design matrix and vector of coefficients
  # Sigma <- diag(rep(1/p, p))
  # X <- MASS::mvrnorm(n = n, mu = rep(0, p), Sigma = Sigma)
  # beta = c(-15,-10,5,10,15, rep(0,p-5))
  # 
  # # --------------------- generate from logistic regression with an intercept of one
  # z <- 1 + X%*%beta
  # pr <- 1/(1 + exp(-z))
  # y <- as.factor(rbinom(n, 1, pr))
  # y = as.numeric(y)-1
  # 
  # # define swag parameters
  # quantile_alpha = .15
  # p_max = 20
  # object = swaglm::swaglm(X=X, y = y, p_max = p_max, family = stats::binomial(),
  #                         alpha = quantile_alpha, verbose = TRUE, seed = 123)
  # 
  # 
  # newdata = X

  # verify that new data has the same structure as X
  if(dim(newdata)[2] != dim(object$X)[2]){
    stop("'newdata' should have the same number of columns as the design matrix X used in the swaglm algorithm")
  }
  
  # compute summary on swaglm object to extract selected models
  summary_swaglm = summary(object)
  
  # define number of new poits to predict
  n_newdata = nrow(newdata)
  
  # define number of selected models
  n_selected_model = dim(summary_swaglm$mat_selected_model)[1]
  
  # create matrix of eta prediction
  mat_eta_prediction = matrix(NA, ncol=n_selected_model, nrow = n_newdata)
  
  for(model_i in seq(n_selected_model)){
    # find variable that compose this model
    var_model_i = as.vector(na.omit(summary_swaglm$mat_selected_model[model_i, ]))
    # subset design matrix
    X_model_i = cbind(rep(1, n_newdata), newdata[, var_model_i])
    estimated_beta_model_i = as.vector(na.omit(summary_swaglm$mat_beta_selected_model[model_i, ]))
    # eta equals X %*% beta, the linear predictor
    eta_model_i = X_model_i %*% estimated_beta_model_i
    mat_eta_prediction[, model_i] = eta_model_i
  }
  
  # compute the response
  mat_reponse_prediction = object$family$linkinv(mat_eta_prediction)

  # structure return object
  res = list("mat_eta_prediction" = mat_eta_prediction,
             "mat_reponse_prediction" = mat_reponse_prediction)
  # return object
  return(res)
}
