stack_models <- function(lst) {
  # drop NULLs
  lst <- Filter(Negate(is.null), lst)

  # find the maximum number of columns
  max_cols <- max(sapply(lst, ncol))

  # pad each matrix with NA to have the same number of columns
  padded <- lapply(lst, function(m) {
    cbind(m, matrix(NA, nrow(m), max_cols - ncol(m)))
  })

  # stack them into one big data frame
  do.call(rbind, padded)
}

extract_by_index <- function(lst_values, lst_index) {
  Map(function(v, idx) v[idx], lst_values, lst_index)
}


#' postprocessing_swaglm
#' @param x An object of class \code{swaglm}.
#' @example /inst/examples/eg_postprocessing_swaglm.R
#' @importFrom stats median
#' @return A list with three elements:
#' \describe{
#'   \item{mat_selected_model}{A matrix where each row represents a selected model and the columns indicate the indices of the variables included in that model.}
#'   \item{vec_aic_selected_model}{A numeric vector containing the AIC values corresponding to each selected model in \code{mat_selected_model}.}
#'   \item{estimated_beta_per_variable}{A named list where each element corresponds to a variable (named \code{V<index>}) and contains a numeric vector of estimated beta coefficients for that variable across all selected models in which it appears.}
#' }
#' @export
postprocessing_swaglm <- function(x) {
  # # Parameters for data generation
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
  # x = swaglm::swaglm(X=X, y = y, p_max = p_max, family = stats::binomial(),
  #                           alpha = quantile_alpha, verbose = TRUE, seed = 123)
  #

  # check that it is indeed a swaglm object
  if (!inherits(x, "swaglm")) {
    stop("Provided object 'x' needs to be of class 'swaglm'")
  }

  # find the lower median AIC
  vec_median_aic <- unlist(lapply(x$lst_AIC, FUN = median))
  id_minimum <- which.min(vec_median_aic)
  min_median_aic <- vec_median_aic[id_minimum]

  # find all model that are below this AIC value
  index_model_below_median_aic <- lapply(x$lst_AIC, function(x) which(x <= min_median_aic))

  # extract models and associated estimated beta
  max_dim <- length(x$lst_var_mat)
  selected_models <- list()
  beta_selected_model <- list()
  for (dim_i in seq(max_dim)) {
    if (length(index_model_below_median_aic[[dim_i]]) == 0) {
      selected_models[[dim_i]] <- NULL
    } else {
      # create matrix to store variables of the selected model at that dimension
      mat_selected_model_dim_i <- matrix(NA, nrow = length(index_model_below_median_aic[[dim_i]]), ncol = dim_i)
      mat_selected_model_dim_i <- x$lst_var_mat[[dim_i]][index_model_below_median_aic[[dim_i]], ] + 1
      # transform to always matrix
      if (is.vector(mat_selected_model_dim_i)) {
        mat_selected_model_dim_i <- t(as.matrix(mat_selected_model_dim_i, nrow = 1))
      }
      # create matrix to store beta
      # beta_selected_model_dim_i = matrix(NA, nrow = length(index_model_below_median_aic[[dim_i]]), ncol = dim_i+1)
      beta_selected_model_dim_i <- x$lst_estimated_beta[[dim_i]][index_model_below_median_aic[[dim_i]], ]
      if (is.vector(beta_selected_model_dim_i)) {
        beta_selected_model_dim_i <- t(as.matrix(beta_selected_model_dim_i, nrow = 1))
      }
      # append to list
      beta_selected_model[[dim_i]] <- beta_selected_model_dim_i
      selected_models[[dim_i]] <- mat_selected_model_dim_i
    }
  }

  # create stacked matrix of selected models
  mat_stacked_model <- stack_models(selected_models)
  # mat_stacked_model =  stack_models(lst_selected_models)

  # get AIC of selected models
  vec_aic_selected_model <- unlist(extract_by_index(x$lst_AIC, index_model_below_median_aic))

  # # verify
  # n_models = dim(mat_stacked_model)[1]
  # vec_aic_refit = vector(mode = "numeric", length(n_models))
  # for(i_model in seq(n_models)){
  #   # i_model = 1
  #   selected_var_i_model = na.omit(mat_stacked_model[i_model, ])
  #   fit = fastglm::fastglm(x = cbind(rep(1, length(y)), x$X[,selected_var_i_model]), y=x$y, family = binomial())
  #   vec_aic_refit[i_model] =  fit$aic
  # }
  # all.equal(vec_aic_selected_model, vec_aic_refit)


  # extract beta
  # beta_selected_model
  # selected_models
  #
  # create a list that contain the estimated beta for each variable that is part of the selected models

  # Get all unique variable indices used across selected models
  all_vars <- unique(na.omit(as.vector(mat_stacked_model)))
  all_vars <- sort(all_vars)

  # Initialize result list
  estimated_beta_per_variable <- vector("list", length(all_vars))
  names(estimated_beta_per_variable) <- paste0("V", all_vars)

  # iterate over each variable
  for (v in all_vars) {
    vals <- c()
    # iterate over elements of list
    for (i in seq_along(beta_selected_model)) {
      if (!is.null(beta_selected_model[[i]])) {
        # case where there are only one model estimated so beta_selected model is a matrix of 1 by something
        if (nrow(beta_selected_model[[i]]) == 1) {
          # identify the column assocaited with the variable
          idx <- which(selected_models[[i]] == v)
          # extract beta (plus 1 because there is always the intercept)
          vals <- c(vals, beta_selected_model[[i]][, (idx + 1)])
        } else if (nrow(beta_selected_model[[i]]) > 1) {
          for (row_i in seq(nrow(beta_selected_model[[i]]))) {
            # row_i=1
            idx <- which(selected_models[[i]][row_i, ] == v)
            # extract beta (plus 1 because there is always the intercept)
            vals <- c(vals, beta_selected_model[[i]][row_i, (idx + 1)])
          }
        }
      }
    }
    estimated_beta_per_variable[[paste0("V", v)]] <- vals
  }


  res <- list(
    "mat_selected_model" = mat_stacked_model,
    "vec_aic_selected_model" = vec_aic_selected_model,
    "estimated_beta_per_variable" = estimated_beta_per_variable
  )
  return(res)
}
