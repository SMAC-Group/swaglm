#ifndef identify_selected_combinations_cpp_H
#define identify_selected_combinations_cpp_H

#include <RcppArmadillo.h>

// Function declaration for compute_quantile
double compute_quantile(const arma::vec& x, double alpha);

// Function declaration for identify_selected_combinations_cpp
Rcpp::List identify_selected_combinations_cpp(const arma::mat& mat_of_variables, const arma::mat& mat_criterion, double alpha = 0.1);

#endif // identify_selected_combinations_cpp_H