// // [[Rcpp::depends(RcppArmadillo)]]
// #include <RcppArmadillo.h>
// using namespace Rcpp;
// 
// // set seed
// // [[Rcpp::export]]
// void my_set_seed(unsigned int seed) {
//   Rcpp::Environment base_env("package:base");
//   Rcpp::Function set_seed_r = base_env["set.seed"];
//   set_seed_r(seed);  
// }
// 
// // [[Rcpp::export]]
// arma::uvec generate_permutation(int n, int m, int seed=123){
//   my_set_seed(seed);
//   arma::uvec random_indices = arma::randperm(n, m);
//   return(random_indices);
// }
// 
// // Generate random indices for sampling `m` rows
// 
// 
// 
// /*** R
// 
// generate_permutation(n = 10, m = 3, seed =  1234)
// 
// */
