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
swag_obj$lst_estimated_beta
swag_obj$lst_AIC
swag_obj$lst_var_mat
swag_obj$lst_selected_models
swag_obj$lst_index_selected_models
?swaglm

for (j in seq_len(max_dim)) {
  index_selected_model <- swag_obj$lst_index_selected_models[[j]]
  index_selected_model_R <- index_selected_model + 1
  
  mat_selected_var_1 <- (swag_obj$lst_var_mat[[j]] + 1)[index_selected_model_R, ]
  mat_selected_var_2 <- (swag_obj$lst_selected_models[[j]] + 1)
  
  # Always flatten to vector before comparing
  res <- all.equal(as.vector(mat_selected_var_1), as.vector(mat_selected_var_2))
  
  print(paste0("j=", j, " -> ", res))
}





