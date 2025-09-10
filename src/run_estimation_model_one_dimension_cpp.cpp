#include "run_estimation_model_one_dimension_cpp.h"

#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace Rcpp;
// Function definition

// This function runs an estimation model in one dimension using the fastglm function
// [[Rcpp::export]]
Rcpp::List run_estimation_model_one_dimension_cpp(const arma::mat& X, const arma::vec& y,  Nullable<List> family, int method) {
  int p = X.n_cols;
  int n = X.n_rows;
  
  // Initialize matrices to store AIC and beta coefficients
  arma::mat mat_AIC_dim_1(p, 1);
  arma::mat mat_beta_dim_1(p, 2);
  arma::mat mat_p_value_dim_1(p, 2);
  
  // Load the fastglm function from the fastglm package
  Environment fastglm_env = Environment::namespace_env("fastglm");
  Function fastglm = fastglm_env["fastglm"];
  Function summary_fastglm = fastglm_env["summary.fastglm"];
  
  // Load the binomial function from the stats package
  Environment stats_env = Environment::namespace_env("stats");
  Function binomial = stats_env["binomial"];
  
  // Set default family to binomial if not provided
  List family_obj;
  if (family.isNull()) {
    family_obj = binomial();
  } else {
    family_obj = as<List>(family);
  }
  
  // Create the design matrix with an intercept column, to avoid recomputing it every time
  arma::mat X_with_intercept(n, p + 1, arma::fill::ones);  // Add a column of ones for the intercept
  X_with_intercept.cols(1, p) = X;  // Fill the rest with the columns from X
  
  for (unsigned int i = 0; i < p; ++i) {
    // Create the design matrix with an intercept column
    // arma::mat X_mat_i = join_rows(arma::ones<arma::vec>(n), X.col(i));
    // Select only the intercept and the ith column
    arma::uvec idx = {0, i + 1};  // Create an index for the first and i-th column
    arma::mat X_mat_i = X_with_intercept.cols(idx); 
    
    // Fit the model using fastglm
    // This is not optimal, fastglm is a R function implemented in Rcpp Eigen, but still we call its R implementation
    List fit = fastglm(Named("x") = X_mat_i, Named("y") = y, Named("family") = family_obj, Named("method") = method);
    List summary_fit = summary_fastglm(fit);
    
    // Extract coefficients matrix as arma::mat
    arma::mat coef_mat = as<arma::mat>(summary_fit["coefficients"]);
    
    
    // Extract AIC, estimated coefficients and p value
    mat_AIC_dim_1(i, 0) = as<double>(fit["aic"]);
    mat_beta_dim_1(i, 0) = as<arma::vec>(fit["coefficients"])(0);
    mat_beta_dim_1(i, 1) = as<arma::vec>(fit["coefficients"])(1);
    
    // P-values (last column of coefficients matrix)
    mat_p_value_dim_1(i, 0) = coef_mat(0, 3); // intercept p-value
    mat_p_value_dim_1(i, 1) = coef_mat(1, 3); // slope p-value
    
    
    
  }
  
  // create matrix of variables for that dimension
  arma::mat matrix_variables(p,1, arma::fill::zeros);
  
  // Fill the first column with values from 0 to n-1
  matrix_variables.col(0) = arma::linspace<arma::vec>(0, p-1, p);
  
  return List::create(Named("mat_AIC_dim_1") = mat_AIC_dim_1,
                      Named("mat_beta_dim_1") = mat_beta_dim_1,
                      Named("matrix_of_variables") = matrix_variables,
                      Named("mat_p_value_dim_1") = mat_p_value_dim_1);
}






/*** R


# # Parameters for data generation
# set.seed(1)
# n <- 2000
# p <- 20
# Sigma <- diag(rep(1/p, p))
# # 
# X = MASS::mvrnorm(n = n, mu = rep(0, p), Sigma = Sigma)
# beta = c(3,5,6,19,70,rep(0,15))
# z <- 1 + X%*%beta
# pr <- 1/(1 + exp(-z))
# y <- as.factor(rbinom(n, 1, pr))
# y = as.numeric(y)-1
# test = run_estimation_model_one_dimension_cpp(X = X, y = y, family = binomial(), method =0 )
# test$mat_p_value_dim_1
# # 
# # 
# # 
# # 
# run_estimation_model_one_dimension = function(X, y, family=binomial()){
# 
#   # X=X_mat_predictors
#   p = dim(X)[2]
#   n = dim(X)[1]
# 
#   # Initialize beta matrices to store beta coef and AIC values
#   mat_AIC_dim_1 = matrix(nrow=p, ncol=1)
#   mat_beta_dim_1 = matrix(nrow=p, ncol=2)
#   mat_p_value_dim_1 = matrix(nrow=p, ncol=2)
# 
# 
# 
#   for(i in 1:p){
#     X_mat_i = cbind(rep(1, n), X[,i])
#     fit = fastglm(x =X_mat_i, y=y, family=family)
#     mat_AIC_dim_1[i,1] = fit$aic
#     mat_beta_dim_1[i, ] = fit$coefficients
#     mat_p_value_dim_1[i,] = summary(fit)$coefficients[,4]
#     
#   }
#   ret=list("mat_AIC_dim_1" = mat_AIC_dim_1,
#            "mat_beta_dim_1" = mat_beta_dim_1,
#            "mat_p_value_dim_1" =mat_p_value_dim_1)
#   return(
#     ret
#   )
# }
# # 
# # 
# res1 = run_estimation_model_one_dimension(X=X, y = y, family=binomial())
# res2 = run_estimation_model_one_dimension_cpp(X=X, y = y, family = binomial(), method = 0)
# res1$mat_AIC_dim_1
# res2$mat_AIC_dim_1
# all.equal(res1$mat_AIC_dim_1, res2$mat_AIC_dim_1)
# all.equal(res1$mat_beta_dim_1, res2$mat_beta_dim_1)
# all.equal(res1$mat_p_value_dim_1,res2$mat_p_value_dim_1)
# # res3 = run_estimation_model_one_dimension_cpp(X=X_mat_predictors, y = y)
# #
# # all.equal(res1, res2)
# # all.equal(res2, res3)
# # all.equal(res1, res2)
# # microbenchmark::microbenchmark(res1 = run_estimation_model_one_dimension(X=X_mat_predictors, y = y),
# #                                res2 = run_estimation_model_one_dimension_cpp(X=X_mat_predictors, y = y), times = 10)
# 
# 
# 
# 

*/

