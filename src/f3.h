#ifndef F3_H
#define F3_H

#include <RcppArmadillo.h>

// Function declaration
Rcpp::List estimate_all_model_combinations_cpp(const arma::mat& X, const arma::vec& y, const arma::mat& matrix_of_variables, Rcpp::Nullable<Rcpp::List> family = R_NilValue, int method = 0);

#endif // F3_H