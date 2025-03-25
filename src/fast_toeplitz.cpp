#include <numeric>
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

//' Generate a Toeplitz Matrix from a Vector
//'
//' This function generates a Toeplitz matrix from a given vector using the Armadillo library.
//' @name fast_toeplitz_matrix_from_vector_cpp
//' @param v A numeric vector.
//' @return A numeric matrix that is the Toeplitz matrix generated from \code{v}.
//' @export
// [[Rcpp::export]]
arma::mat fast_toeplitz_matrix_from_vector_cpp(const arma::vec& v) {
  return arma::toeplitz(v);
}
 