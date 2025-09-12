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


#' summary.swaglm
#' @param object An object of class \code{swaglm}.
#' @param ... Additional parameters
#' @name summary.swaglm
#' @method summary swaglm
#' @example /inst/examples/eg_summary.swaglm.R
#' @importFrom stats median
#'
#' @return A list of class \code{summary_swaglm} with five elements:
#' \describe{
#'   \item{mat_selected_model}{A matrix where each row represents a selected model. Columns give the indices of the variables included in that model. Models are padded with \code{NA} to match the largest model size.}
#'   \item{mat_beta_selected_model}{A matrix containing the estimated regression coefficients (including intercept) of the selected models, stacked across all dimensions. Only models with AIC less than or equal to the **lowest median AIC across all dimensions** are included. Each row corresponds to a model, columns correspond to the coefficients. Rows are padded with \code{NA} to match the largest model size.}
#'   \item{mat_p_value_selected_model}{A matrix containing the p-values associated with the estimated regression coefficients in \code{beta_selected_model}, stacked in the same order. Each row corresponds to a model, columns correspond to the coefficients. Rows are padded with \code{NA} to match the largest model size.}
#'   \item{vec_aic_selected_model}{A numeric vector containing the AIC values of all models in \code{mat_selected_model}, stacked across dimensions. These are the AIC values for the selected models that passed the threshold described above.}
#'   \item{lst_estimated_beta_per_variable}{A named list where each element corresponds to a variable (named \code{V<index>}). Each element is a numeric vector containing all estimated beta coefficients for that variable across all selected models in which it appears. This summarizes the distribution of effects for each variable across the selected models.}
#'   \item{lst_p_value_per_variable}{A named list where each element corresponds to a variable (named \code{V<index>}). Each element is a numeric vector containing all estimated p-values for that variable across all selected models in which it appears.}
#' }
#'
#' \strong{Model selection criterion:}
#' For each model dimension (number of variables in the model), the median AIC across all models of that dimension is computed. The **smallest median AIC across all dimensions** is identified. Then, **all models with AIC less than or equal to this value** are selected. This ensures that only relatively well-performing models across all dimensions are retained for summarization.
#' @export
summary.swaglm <- function(object, ...) {
  # Parameters for data generation
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
  #                           alpha = quantile_alpha, verbose = TRUE, seed = 123)
  # 
  

  # check that it is indeed a swaglm object
  if (!inherits(object, "swaglm")) {
    stop("Provided object 'object' needs to be of class 'swaglm'")
  }
  # find the lower median AIC
  vec_median_aic <- unlist(lapply(object$lst_AIC, FUN = median))
  id_minimum <- which.min(vec_median_aic)
  min_median_aic <- vec_median_aic[id_minimum]

  # find all model that are below this AIC value
  index_model_below_median_aic <- lapply(object$lst_AIC, function(x) which(x <= min_median_aic))

  # extract models and associated estimated beta and p value
  max_dim <- length(object$lst_var_mat)
  selected_models <- list()
  beta_selected_model <- list()
  p_value_selected_model <- list()

  for (dim_i in seq(max_dim)) {
    if (length(index_model_below_median_aic[[dim_i]]) == 0) {
      selected_models[[dim_i]] <- NULL
      beta_selected_model[[dim_i]] <- NULL
      p_value_selected_model[[dim_i]] <- NULL
    } else {
      # create matrix to store variables of the selected model at that dimension
      # mat_selected_model_dim_i <- matrix(NA, nrow = length(index_model_below_median_aic[[dim_i]]), ncol = dim_i)
      mat_selected_model_dim_i <- object$lst_var_mat[[dim_i]][index_model_below_median_aic[[dim_i]], ] 
      # transform to always matrix
      if (is.vector(mat_selected_model_dim_i)) {
        mat_selected_model_dim_i <- t(as.matrix(mat_selected_model_dim_i, nrow = 1))
      }
      # create matrix to store beta
      # beta_selected_model_dim_i = matrix(NA, nrow = length(index_model_below_median_aic[[dim_i]]), ncol = dim_i+1)
      beta_selected_model_dim_i <- object$lst_estimated_beta[[dim_i]][index_model_below_median_aic[[dim_i]], ]
      p_value_selected_model_dim_i <- object$lst_p_value[[dim_i]][index_model_below_median_aic[[dim_i]], ]
      if (is.vector(beta_selected_model_dim_i)) {
        beta_selected_model_dim_i <- t(as.matrix(beta_selected_model_dim_i, nrow = 1))
        p_value_selected_model_dim_i <- t(as.matrix(p_value_selected_model_dim_i, nrow = 1))
      }
      # append to list
      beta_selected_model[[dim_i]] <- beta_selected_model_dim_i
      p_value_selected_model[[dim_i]] <- p_value_selected_model_dim_i
      selected_models[[dim_i]] <- mat_selected_model_dim_i
    }
  }


  # create stacked matrix of selected models
  mat_stacked_model <- stack_models(selected_models)

  # get AIC of selected models
  vec_aic_selected_model <- unlist(extract_by_index(object$lst_AIC, index_model_below_median_aic))


  # # # --------------------------------------------------test
  # n_models <- nrow(mat_stacked_model)
  # 
  # vec_aic_refit   <- numeric(n_models)
  # beta_refit_list <- vector("list", n_models)
  # pval_refit_list <- vector("list", n_models)
  # 
  # for (i_model in seq_len(n_models)) {
  #   # i_model = 1
  #   selected_var_i_model <- na.omit(mat_stacked_model[i_model, ])
  #   X_mat <- cbind(1, object$X[, selected_var_i_model, drop = FALSE])
  # 
  #   fit <- fastglm::fastglm(
  #     x = X_mat,
  #     y = object$y,
  #     family = binomial()
  #   )
  # 
  #   # AIC
  #   vec_aic_refit[i_model] <- fit$aic
  # 
  #   # coefficients & p-values
  #   summ <- summary(fit)
  #   beta_refit_list[[i_model]] <- summ$coefficients[, 1]   # estimates
  #   pval_refit_list[[i_model]] <- summ$coefficients[, 4]   # p-values
  # }
  # 
  # # Compare AIC
  # all.equal(vec_aic_selected_model, vec_aic_refit)
  # 
  # 
  # stack_matrices <- function(lst, fill = NA_real_) {
  #   if (!is.list(lst)) stop("'lst' must be a list")
  # 
  #   # remove NULL entries
  #   lst2 <- lst[!vapply(lst, is.null, logical(1))]
  #   if (length(lst2) == 0) return(matrix(nrow = 0, ncol = 0))
  # 
  #   # reject nested lists (ask user to convert them)
  #   is_bad <- vapply(lst2, function(el) is.list(el) && !is.data.frame(el), logical(1))
  #   if (any(is_bad)) {
  #     stop("List contains nested lists. Convert those elements to a matrix, data.frame or numeric vector first.")
  #   }
  # 
  #   # coerce each element to a matrix:
  #   mats <- lapply(lst2, function(el) {
  #     if (is.data.frame(el)) return(as.matrix(el))
  #     if (is.matrix(el)) return(el)
  #     if (is.numeric(el) && is.null(dim(el))) return(matrix(el, nrow = 1))
  #     stop("Unsupported element type: ", paste(class(el), collapse = ", "))
  #   })
  # 
  #   # pad to the maximum number of columns
  #   max_cols <- max(vapply(mats, ncol, integer(1)))
  #   padded <- lapply(mats, function(m) {
  #     nc <- ncol(m)
  #     if (nc < max_cols) {
  #       cbind(m, matrix(fill, nrow = nrow(m), ncol = max_cols - nc))
  #     } else m
  #   })
  # 
  #   # bind and return
  #   do.call(rbind, padded)
  # }
  # 
  # test1 = stack_matrices(beta_selected_model)
  # test2 = stack_matrices(beta_refit_list)
  # all.equal(test1,test2)
  # 
  # test3 = stack_matrices(p_value_selected_model)
  # test4 = stack_matrices(pval_refit_list)
  # all.equal(test3,test4)


  # --------------------------------------------------test


  # Get all unique variable indices used across selected models
  all_vars <- unique(na.omit(as.vector(mat_stacked_model)))
  all_vars <- sort(all_vars)

  # Initialize result list
  lst_estimated_beta_per_variable <- vector("list", length(all_vars))
  lst_p_value_per_variable <- vector("list", length(all_vars))
  
  names(lst_estimated_beta_per_variable) <- paste0("V", all_vars)
  names(lst_p_value_per_variable) <- paste0("V", all_vars)
  
  # iterate over each variable
  for (v in all_vars) {

    # v=1
    beta_vals <- c()
    p_value_vals <- c()
    # iterate over elements of list
    for (i in seq_along(beta_selected_model)) {

      if (!is.null(beta_selected_model[[i]])) {
        # case where there are only one model estimated so beta_selected model is a matrix of 1 by something
        if (nrow(beta_selected_model[[i]]) == 1) {
          # identify the column assocaited with the variable
          idx <- which(selected_models[[i]] == v)
          # extract beta (plus 1 because there is always the intercept)
          beta_vals <- c(beta_vals, beta_selected_model[[i]][, (idx + 1)])
          p_value_vals <- c(p_value_vals, p_value_selected_model[[i]][, (idx + 1)])
          
          
        } else if (nrow(beta_selected_model[[i]]) > 1) {
          
          
          for (row_i in seq(nrow(beta_selected_model[[i]]))) {
            # row_i=1
            idx <- which(selected_models[[i]][row_i, ] == v)
            # extract beta (plus 1 because there is always the intercept)
            beta_vals <- c(beta_vals, beta_selected_model[[i]][row_i, (idx + 1)])
            p_value_vals <- c(p_value_vals, p_value_selected_model[[i]][row_i, (idx + 1)])
         
          }
        }
      }
    }
    lst_estimated_beta_per_variable[[paste0("V", v)]] <- beta_vals
    lst_p_value_per_variable[[paste0("V", v)]] <- p_value_vals
  }

  # stack estimated beta and estimated p value
  mat_beta_selected_model <- stack_models(beta_selected_model)
  mat_p_value_selected_model <- stack_models(p_value_selected_model)

  # return object
  res <- list(
    "mat_selected_model" = mat_stacked_model,
    "mat_beta_selected_model" = mat_beta_selected_model,
    "mat_p_value_selected_model" = mat_p_value_selected_model,
    "vec_aic_selected_model" = vec_aic_selected_model,
    "lst_estimated_beta_per_variable" = lst_estimated_beta_per_variable,
    "lst_p_value_per_variable" = lst_p_value_per_variable
  )
  # define class
  class(res) <- "summary_swaglm"

  # return
  return(res)
}




#
#
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
# x = swaglm::swaglm(X=X, y = y, p_max = p_max, family = stats::binomial(),
#                           alpha = quantile_alpha, verbose = TRUE, seed = 123)
# test=summary(x)
