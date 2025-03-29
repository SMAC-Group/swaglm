#include "f3.h"

#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace Rcpp;


// Function to estimate all model combinations
// [[Rcpp::export]]
List estimate_all_model_combinations_cpp(const arma::mat& X, const arma::vec& y, const arma::mat& matrix_of_variables,
                                         Nullable<List> family, int method) {
  // Load the fastglm function from the fastglm package
  Environment fastglm_env = Environment::namespace_env("fastglm");
  Function fastglm = fastglm_env["fastglm"];
  
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
  
  // Get the number of models and dimensions
  int n_models = matrix_of_variables.n_rows;
  int dimension_of_models = matrix_of_variables.n_cols;
  int n = X.n_rows;
  
  // Create the design matrix with an intercept column (of ones)
  int p = X.n_cols;
  arma::mat X_with_intercept(n, p + 1, arma::fill::ones);  // Add column of ones for intercept
  X_with_intercept.cols(1, p) = X;  // Fill the rest with columns from X
  
  
  // Initialize matrices to store AIC values and beta coefficients
  arma::mat mat_AIC(n_models, 1);
  arma::mat mat_beta(n_models, dimension_of_models + 1);
  
  // Iterate over each model combination
  for (int i = 0; i < n_models; ++i) {
    // Extract the variables for the current model
    arma::vec variables_mod_i = matrix_of_variables.row(i).t();
    
    // Create an index vector for the intercept (0) and the selected variables
    arma::uvec idx = {0};  // Start with the intercept (0-th column)
    
    // Convert variables_mod_i to uvec and add 1 (to account for intercept) 
    arma::uvec variables_mod_i_uvec = arma::conv_to<arma::uvec>::from(variables_mod_i);
    arma::uvec idx2 = join_cols(idx, variables_mod_i_uvec + 1);  // Add variables to the index (shifted by 1)
    
    // Create a subview of X_with_intercept for the intercept and selected columns
    arma::mat X_mat_mod_i = X_with_intercept.cols(idx2);
    
    // Create the design matrix with an intercept column, here we get the columns minus 1 ! so index should start at 1 if so
    // arma::mat X_mat_mod_i = join_rows(arma::ones<arma::vec>(n), X.cols(arma::conv_to<arma::uvec>::from(variables_mod_i - 1)));
    
    // Fit the model using fastglm
    List fit = fastglm(Named("x") = X_mat_mod_i, Named("y") = y, Named("family") = family_obj, Named("method") = method);
    
    // Store the AIC and coefficients
    mat_AIC(i, 0) = as<double>(fit["aic"]);
    mat_beta.row(i) = as<arma::vec>(fit["coefficients"]).t();
  }
  
  // Return the results as a list
  return List::create(Named("mat_AIC") = mat_AIC, Named("mat_beta") = mat_beta);
}














// 
// /*** R
// 
// 
// # Parameters for data generation
// set.seed(1)
// n <- 2000
// p <- 20
// Sigma <- diag(rep(1/p, p))
// 
// X <- MASS::mvrnorm(n = n, mu = rep(0, p), Sigma = Sigma)
// beta = c(-10,5,6,19,70,rep(0,15))
// 
// # --------------------- logistic reg
// z <- 1 + X%*%beta
// pr <- 1/(1 + exp(-z))
// y <- as.factor(rbinom(n, 1, pr))
// y = as.numeric(y)-1
// 
// matrix_of_variables = matrix(c(sample(x = 1:p, size=14)), ncol=2, byrow = T)
// 
// estimate_all_model_combinations = function(X, y, matrix_of_variables, family=stats::binomial(), method=0){
// 
// 
// 
// 
// 
//   # each row of matrix_of_variable is a combination of variables, there are nrow(matrix_of_variable) models to evaluate
//   n_models = nrow(matrix_of_variables)
//   dimension_of_models = ncol(matrix_of_variables)
//   n=dim(X)[1]
//   # Initialize beta matrices to store beta coef and AIC values
//   mat_AIC = matrix(nrow=n_models, ncol=1)
//   mat_beta = matrix(nrow=n_models, ncol=dimension_of_models+1)
// 
//   for(i in 1:n_models){
//     # i=1
//     variables_mod_i = matrix_of_variables[i,]
// 
//     # create X matrix
//     X_mat_mod_i = cbind(rep(1, n),
//                         X[,variables_mod_i])
//     fit = fastglm::fastglm(x =X_mat_mod_i, y=y, family=family, method = method)
//     mat_AIC[i,1] = fit$aic
//     mat_beta[i, ] = fit$coefficients
//   }
// 
//   ret=list("mat_AIC" = mat_AIC,
//            "mat_beta" = mat_beta)
// 
// 
//   return(ret)
// 
// 
// }
// 
// 
// 
// 
// res1 = estimate_all_model_combinations(X = X, y, matrix_of_variables = (matrix_of_variables+1), family =binomial(), method = 0)
// res2 = estimate_all_model_combinations_cpp(X = X, y, matrix_of_variables = matrix_of_variables, family =binomial(), method = 0)
// 
// all.equal(res1, res2)
// res1$mat_AIC
// res2$mat_AIC
// all.equal(res1, res2)
// # 
// # microbenchmark::microbenchmark(res1 = estimate_all_model_combinations(X = X, y, matrix_of_variables = matrix_of_variables, family =binomial(), method = 0),
// #                                res2 = estimate_all_model_combinations_cpp(X = X, y, matrix_of_variables = matrix_of_variables, family =binomial(), method = 0))
// 
// 
// 
// 
// 
// */

