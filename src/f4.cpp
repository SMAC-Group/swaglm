#include "f4.h"
#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace Rcpp;


// Given a matrix, remove rows where all entries are not unique
// [[Rcpp::export]]
arma::mat removeRowsWithDuplicates(const arma::mat& inputMatrix) {
  // Create a new matrix to store rows without duplicates
  arma::mat uniqueRowsMatrix;
  
  // Check each row for duplicates and add to uniqueRowsMatrix if none are found
  for (size_t i = 0; i < inputMatrix.n_rows; ++i) {
    arma::rowvec row = inputMatrix.row(i);
    std::unordered_set<int> seenValues(row.begin(), row.end());
    
    // If the size of the set matches the number of elements in the row, there are no duplicates
    if (seenValues.size() == row.n_elem) {
      uniqueRowsMatrix.insert_rows(uniqueRowsMatrix.n_rows, row);
    }
  }
  
  return uniqueRowsMatrix;
}




// [[Rcpp::export]]
arma::mat removeDuplicateRowsRegardlessOfOrder(const arma::mat& inputMatrix) {
  // Create a set to track unique rows
  std::unordered_set<std::string> uniqueRows;
  arma::mat uniqueMatrix;
  
  // Iterate over each row
  for (size_t i = 0; i < inputMatrix.n_rows; ++i) {
    arma::rowvec row = inputMatrix.row(i);
    std::vector<int> sortedRow(row.begin(), row.end());
    std::sort(sortedRow.begin(), sortedRow.end());
    
    // Convert the sorted row to a string for hashing
    std::string rowKey;
    for (size_t j = 0; j < sortedRow.size(); ++j) {
      rowKey += std::to_string(sortedRow[j]) + ",";
    }
    
    // Check if the row is unique
    if (uniqueRows.find(rowKey) == uniqueRows.end()) {
      uniqueRows.insert(rowKey);
      uniqueMatrix.insert_rows(uniqueMatrix.n_rows, row);
    }
  }
  
  return uniqueMatrix;
}


// [[Rcpp::export]]
arma::mat sort_rows(const arma::mat& X) {
  arma::mat Y = X; // Copy input matrix to avoid modifying the original
  
  for (size_t i = 0; i < Y.n_rows; ++i) {
    Y.row(i) = arma::sort(Y.row(i)); // Sort each row
  }
  
  return Y;
}

// Function that given a matrix of selected variabls combination and id of screened variables in dimension 1, return the matrix of possible model combination for next dimension 
// [[Rcpp::export]]
arma::mat compute_all_possible_variable_combinations_cpp(const arma::mat& originalMatrix, const arma::vec& idScreening) {
  int nrv = originalMatrix.n_rows;
  int d = originalMatrix.n_cols + 1; // New number of columns
  int numScreening = idScreening.n_elem;
  
  // Create a new matrix with nrv * numScreening rows and d columns
  arma::mat stackedMatrix(nrv * numScreening, d);
  
  // Fill the first d-1 columns with replicated originalMatrix
  for (int i = 0; i < numScreening; ++i) {
    stackedMatrix.submat(i * nrv, 0, (i + 1) * nrv - 1, d - 2) = originalMatrix;
  }
  
  // Fill the last column with repeated idScreening values
  for (int i = 0; i < numScreening; ++i) {
    stackedMatrix.col(d - 1).subvec(i * nrv, (i + 1) * nrv - 1).fill(idScreening[i]);
  }
  // remove identical rows
  arma::mat ret_mat = removeRowsWithDuplicates(stackedMatrix);
  // remove rows that are not combinations of unique variables (if there is more than once the variable per row)
  arma::mat ret_mat_2 = removeDuplicateRowsRegardlessOfOrder(ret_mat);
  // sort the variables in each row, not strictly necessary but to compare exactly with algorithm in R
  arma::mat ret_mat_3 = sort_rows(ret_mat_2);
  return(ret_mat_3);
  
}



// binomial coefficient

// Function to compute the binomial coefficient (n choose k)
// [[Rcpp::export]]
int binomial_coefficient(int n, int k) {
  if (k > n) return 0;
  if (k == 0 || k == n) return 1;
  
  k = std::min(k, n - k); // Take advantage of symmetry
  long long c = 1;
  for (int i = 0; i < k; ++i) {
    c = c * (n - i) / (i + 1);
  }
  return static_cast<int>(c);
}









// 
//  // You can include R code blocks in C++ files processed with sourceCpp
//  // (useful for testing and development). The R code will be automatically
//  // run after the compilation.
//  //
// 
//  /*** R
// 
// 
// compute_all_possible_variable_combinations = function(id_var, id_screening){
// 
//   # id_var =matrix(c(1,2,3,4,5,6,7,8,9), ncol=3, byrow = T)
//   # id_screening=3:12
//   # d=4
// 
//   # get number of potential combinations
//   nrv = dim(id_var)[1]
//   d = dim(id_var)[2]+1
//   # create matrix with number of row corresponding to number of selected model from previous dimension plus always all the screened variables from first step
//   A <- matrix(nr = nrv*length(id_screening), nc = d)
//   # create all first 1 to 2 columns the extracted model so basically replicating the matrix id_var
//   A[, 1:(d - 1)] <- kronecker(cbind(rep(1, length(id_screening))), id_var)
//   # fill in last column with repeating nrv times each values in id screening
//   A[, d] <- rep(id_screening, each = nrv)
//   B <- unique(t(apply(A, 1, sort))) # deletes the repeated rows
//   id_ndup <- which(apply(B, 1, anyDuplicated) == 0) # removes the models with same Xi
//   B[id_ndup,]
// 
// 
// }
// 
// 
// 
// id_var =matrix(c(1,2,3,4,5,23), ncol=2, byrow = T)
// id_var
// id_screening = c(1,2,3,7,8)
// id_screening
// res1 = compute_all_possible_variable_combinations(id_var = id_var, id_screening = id_screening)
// res2 = compute_all_possible_variable_combinations_cpp(originalMatrix = id_var, idScreening = id_screening)
// res1
// res2
// all.equal(res1, res2)
// 
// 
// microbenchmark::microbenchmark(res1 = compute_all_possible_variable_combinations(id_var = id_var, id_screening = id_screening),
//                                res2 = compute_all_possible_variable_combinations_cpp(originalMatrix = id_var, idScreening = id_screening))
// 
// # test if first dimension
// id_var2 = matrix(seq(10), nrow=10)
// id_screening=c(15:25)
// compute_all_possible_variable_combinations(id_var = id_var2, id_screening = id_screening)
// 
// 
// 
// # areAllRowsInMatrixRegardlessOfOrder <- function(mat1, mat2) {
// #   # Sort each row of both matrices
// #   sortedMat1 <- t(apply(mat1, 1, sort))
// #   sortedMat2 <- t(apply(mat2, 1, sort))
// #
// #   # Check if all rows in sortedMat1 are present in sortedMat2
// #   allRowsPresent <- all(apply(sortedMat1, 1, function(row) {
// #     any(apply(sortedMat2, 1, function(row2) {
// #       all(row == row2)
// #     }))
// #   }))
// #
// #   return(allRowsPresent)
// # }
// # areAllRowsInMatrixRegardlessOfOrder(res1, res2)
// */
// 
// 
// 

