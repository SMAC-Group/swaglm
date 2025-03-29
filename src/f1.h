#ifndef F1_H
#define F1_H

#include <RcppArmadillo.h>

// Function declaration
Rcpp::List run_estimation_model_one_dimension_cpp(const arma::mat& X, const arma::vec& y, Rcpp::Nullable<Rcpp::List> family = R_NilValue, int method = 0);

#endif // F1_H