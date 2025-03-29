#ifndef F4_H
#define F4_H

#include <RcppArmadillo.h>

// Function declarations
arma::mat removeRowsWithDuplicates(const arma::mat& inputMatrix);
arma::mat removeDuplicateRowsRegardlessOfOrder(const arma::mat& inputMatrix);
arma::mat sort_rows(const arma::mat& X);
arma::mat compute_all_possible_variable_combinations_cpp(const arma::mat& originalMatrix, const arma::vec& idScreening);
int binomial_coefficient(int n, int k);

#endif // F4_H
