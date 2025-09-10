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
swaglm_obj = swaglm::swaglm(X=X, y = y, p_max = p_max, family = stats::binomial(),
                          alpha = quantile_alpha, verbose = TRUE, seed = 123)
swaglm_obj
swaglm_summ = summary(swaglm_obj)
swaglm_summ$mat_selected_model
swaglm_summ$mat_beta_selected_model
swaglm_summ$mat_p_value_selected_model
swaglm_summ$vec_aic_selected_model
swaglm_summ$lst_estimated_beta_per_variable


# plot distribution of estimated beta with respect to true value
for(i in seq_along(swaglm_summ$lst_estimated_beta_per_variable)){
  # get variable name

  var_name_i = names(swaglm_summ$lst_estimated_beta_per_variable)[i]
  var_index_i =  as.numeric(substr(var_name_i, 2, nchar(var_name_i)))
  boxplot(swaglm_summ$lst_estimated_beta_per_variable[[i]], 
          main=paste0(var_name_i, " ", "true value=", beta[i]))
}
