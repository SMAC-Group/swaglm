# --------------------------------------------
# --- dev swaglm


rm(list=ls())

# library(fastglm)
library(swaglm)

# Parameters for data generation
set.seed(1)
n <- 2000
p <- 20
Sigma <- diag(rep(1/p, p))

X_mat_predictors = MASS::mvrnorm(n = n, mu = rep(0, p), Sigma = Sigma)
beta_predictors = c(2,5,6,19,70, rep(0,p-5))
X_generate <- cbind(rep(1, n), X_mat_predictors)
beta = c(1, beta_predictors)

# ----------------------------- gaussian iid regression

sigma2 = 5
y = X_generate %*% beta + rnorm(n, sd = sqrt(sigma2))
run_estimation_model_one_dimension_cpp(X_mat_predictors, y = y, family = gaussian())


#---------------------------------------- logistic regression
z <- 1 + X_generate%*%beta
pr <- 1/(1 + exp(-z))
y <- as.factor(rbinom(n, 1, pr))
y = as.numeric(y)-1




run_estimation_model_one_dimension = function(X, y, family=binomial()){
  
 # X=X_mat_predictors
  p = dim(X)[2]
  n = dim(X)[1]
  
  # Initialize beta matrices to store beta coef and AIC values
  mat_AIC_dim_1 = matrix(nrow=p, ncol=1)
  mat_beta_dim_1 = matrix(nrow=p, ncol=2)
  
  for(i in 1:p){
    X_mat_i = cbind(rep(1, n), X[,i])
    fit = fastglm::fastglm(x =X_mat_i, y=y, family=family)
    mat_AIC_dim_1[i,1] = fit$aic
    mat_beta_dim_1[i, ] = fit$coefficients
  }
  ret=list(mat_AIC_dim_1,
           mat_beta_dim_1)
  return(
    ret
  )
}

res1 = run_estimation_model_one_dimension(X=X_mat_predictors, y = y)
res2 = run_estimation_model_one_dimension_cpp(X_mat_predictors, y = y)
res1
res2
