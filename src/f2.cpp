#include "f2.h"

#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace Rcpp;



// Function to compute the quantile of a vector, similar to R's quantile function
double compute_quantile(const arma::vec& x, double alpha) {
  arma::vec sorted_x = x;
  std::sort(sorted_x.begin(), sorted_x.end());
  
  int n = sorted_x.n_elem;
  double pos = alpha * (n - 1);
  
  int k = static_cast<int>(std::floor(pos));
  double d = pos - k;
  
  if (k + 1 < n) {
    return (1 - d) * sorted_x[k] + d * sorted_x[k + 1];
  } else {
    return sorted_x[k];
  }
}


// Function to identify selected combinations
// given a matrix of variables and a matrix of one column/vector of scores associated with each combination of variables (rows of the first matrix), return the matrix of combinations for which the criterion is smaller than the alpha-th quantile of the score of these models 
// [[Rcpp::export]]
List identify_selected_combinations_cpp(const arma::mat& mat_of_variables, const arma::mat& mat_criterion, double alpha) {
  // Calculate the quantile of the first column of mat_criterion
  double q_star = compute_quantile(mat_criterion.col(0), alpha);
  
  // Identify rows where the first column's value is less than or equal to q_star
  arma::uvec id_selected_models = find(mat_criterion.col(0) <= q_star);
  
  // Extract the selected rows from mat_of_variables
  arma::mat mat_of_variables_selected_models = mat_of_variables.rows(id_selected_models);
  
  // Return the results as a list
  return List::create(
    Named("id_selected_models") = id_selected_models , 
    Named("mat_of_variables_selected_models") = mat_of_variables_selected_models
  );
}

