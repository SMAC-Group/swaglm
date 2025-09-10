rm(list=ls())
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
swag_obj = swaglm::swaglm(X=X, y = y, p_max = p_max, family = stats::binomial(), alpha = quantile_alpha, verbose = TRUE, seed = 123)

# verify the p value for dim 1
swag_obj$lst_p_value[[1]]
test_p_value_dim_1 = matrix(NA, ncol=2, nrow=p)
for(i in 1:p){

      X_mat_i = cbind(rep(1, n), X[,i])
      fit = fastglm(x =X_mat_i, y=y, family=binomial())
      test_p_value_dim_1[i,] = summary(fit)$coefficient[,4]

    }
all.equal(test_p_value_dim_1, swag_obj$lst_p_value[[1]])

# test p value for dim 2 and greater
matrix_of_variables = swag_obj$lst_var_mat[[2]]+1
test_p_value_dim_2 = matrix(NA, ncol=3, nrow=dim(matrix_of_variables)[1])
for(i in 1:(dim(matrix_of_variables)[1])){

  column_to_subset = matrix_of_variables[i, ]
  X_mat_i = cbind(rep(1, n), X[,column_to_subset])
  fit = fastglm(x =X_mat_i, y=y, family=binomial())
  test_p_value_dim_2[i,] = summary(fit)$coefficient[,4]
  
}
all.equal(test_p_value_dim_2, swag_obj$lst_p_value[[2]])
# check for all dimension

# Check p-values across all dimensions

n <- nrow(X)
n_dim <- length(swag_obj$lst_var_mat)

for (d in seq_len(n_dim)) {
  # d=1
  matrix_of_variables <- swag_obj$lst_var_mat[[d]] + 1L  # +1 for R indexing
  n_models <- nrow(matrix_of_variables)
  n_coef <- ncol(matrix_of_variables) + 1  # intercept + predictors
  
  test_p_values <- matrix(NA_real_, nrow = n_models, ncol = n_coef)
  
  for (i in seq_len(n_models)) {
    # i=1
    cols <- matrix_of_variables[i, ]
    X_mat_i <- cbind(1, X[, cols])
    
    fit <- fastglm::fastglm(x = X_mat_i, y = y, family = binomial())
    test_p_values[i, ] <- summary(fit)$coefficients[, 4]
  }
  
  msg <- paste0("Dimension ", d, ": ",
                isTRUE(all.equal(test_p_values, swag_obj$lst_p_value[[d]])))
  message(msg)
}


# verify that var mat and selected model do match
for (j in seq_len(max_dim)) {
  index_selected_model <- swag_obj$lst_index_selected_models[[j]]
  index_selected_model_R <- index_selected_model + 1
  
  mat_selected_var_1 <- (swag_obj$lst_var_mat[[j]] + 1)[index_selected_model_R, ]
  mat_selected_var_2 <- (swag_obj$lst_selected_models[[j]] + 1)
  
  # Always flatten to vector before comparing
  res <- all.equal(as.vector(mat_selected_var_1), as.vector(mat_selected_var_2))
  
  print(paste0("j=", j, " -> ", res))
}





