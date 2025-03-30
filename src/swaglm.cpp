// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;


// Load all required functions

// ------------------------------------------------ estimation model one dimension / screening


// run_estimation_model_one_dimension_cpp

#include "f1.h"
 
// ------------------------------------------------------- identify selected models to keep for next dimension
 
// identify_selected_combinations_cpp

#include "f2.h"
 
// ---------------------------------------------------------------- estimate all models combinations provided

// estimate_all_model_combinations_cpp

#include "f3.h"


// ------------------- create matrix of combinations of models to explore

#include "f4.h"

// ---------------------------------- SWAG algorithm

//' sawglm 
//'
//' Run the SWAG algorithm supporting Generalized Linear Models specified by a family object and using fastglm library.
//' @name swaglm
//' @param X A numeric \code{matrix} of predictors.
//' @param y A numeric \code{vector} of responses.
//' @param family A \code{family} object. Default is binomial.
//' @param p_max A \code{integer} specifying the maximum dimension to explore
//' @param method An \code{integer} scalar with value 0 for the column-pivoted QR decomposition, 1 for the unpivoted QR decomposition, 2 for the LLT Cholesky, or 3 for the LDLT Cholesky
//' @param alpha a \code{double} specifying the quantile of the criterion used to select models which are employed to construct models to explore at the next dimension
//' @param verbose A \code{boolean} used to control verbose
//' @param seed A \code{int} that is the random seed used when creating the set of model to explore for the next dimension
//' @return A list of list of matrices. The first list corresponds to the estimated coefficients per models per dimension. At each dimension, the matrix of estimated coefficient including the intercept for each models.
//' models are listed per lines and estimated coefficients on each columns. Second list is the list of AIC criterion values for evaluated models at eahc dimension. Each entries of that List is a matrix of a single column where for each row is indicated the criterion value evaluated at this models.
//' Third list is the list of matrices that specify the combinations of variables explored at each dimension. On this list, each entry correspond to a dimension, the matrix for each dimension indicates on the row the combinations of variables explored and on the columns the variables composing each model represented per each row.
//' Fouth entry on the List of Lists is the List corresponding to xxx.
//' @export
// [[Rcpp::export]]
List swaglm(const arma::mat& X, const arma::vec& y, int p_max=2, Nullable<List> family = R_NilValue, int method = 0, double alpha=0.3, bool verbose = false, int seed  = 123) {
  
  // create List objects of the estimated beta, the AICs, the varmat and the selected models
  List lst_estimated_beta;
  List lst_AIC;
  List lst_var_mat;
  List lst_selected_models;
  List lst_index_selected_models;
  
  //-------------- first we run screening procedure
  List res_screening = run_estimation_model_one_dimension_cpp(X, y, family, method);
  
  // identify the selected variables from screening
  List res2 = identify_selected_combinations_cpp(as<arma::mat>(res_screening["matrix_of_variables"]),  as<arma::mat>(res_screening["mat_AIC_dim_1"]), alpha);
  
  // //  Store estimated beta and AIC for dimension 1
  lst_estimated_beta.push_back(res_screening["mat_beta_dim_1"]);
  lst_AIC.push_back(res_screening["mat_AIC_dim_1"]);
  lst_var_mat.push_back(res_screening["matrix_of_variables"]);
  lst_selected_models.push_back(res2["mat_of_variables_selected_models"]);
  lst_index_selected_models.push_back(res2["id_selected_models"]);
  
  // save variables identified in screening
  arma::vec variables_screening = as<arma::vec>(res2["mat_of_variables_selected_models"]).col(0);
  
  // Compute the number of variables (p) in X and define the number of maximum explored model per dimension as number of selected variables in first screening step choose 2, so the number of models explored at dimension 2
  int p = X.n_cols;
  int m = binomial_coefficient(variables_screening.n_elem, 2); // Now using C++ function
  
  if (verbose) {
    Rcpp::Rcout << "Completed models of dimension 1" << "\n";
  }
  
  // run procedure for 2 to pmax
  for (int dimension_model = 2; dimension_model <= p_max; ++dimension_model) {
    // create matrix of new combinations of models to explore
    arma::mat res3 = compute_all_possible_variable_combinations_cpp(as<arma::mat>(res2["mat_of_variables_selected_models"]), variables_screening);
    // Select at most m models to evaluate randomly
    if (res3.n_rows > m) {
      // Set the seed for reproducibility
      arma::arma_rng::set_seed(seed + dimension_model);
      // Generate random indices for sampling `m` rows
      arma::uvec random_indices = arma::randperm(res3.n_rows, m);
      res3 = res3.rows(random_indices); // Select random rows based on indices
    }
    
    // Check if res3 is empty and break if true
    if (res3.n_rows == 0) {
      break;
    }
    
    // run estimation on all possible models
    List res4 = estimate_all_model_combinations_cpp(X, y, res3, family, method);
    
    // identify selected combination of variables based on quantile
    List new_res2 = identify_selected_combinations_cpp(res3, as<arma::mat>(res4["mat_AIC"]), alpha);
    
    // update res 2
    res2 = new_res2;
    
    // save for each dimension
    lst_estimated_beta.push_back( res4["mat_beta"]);
    lst_AIC.push_back(res4["mat_AIC"]);
    lst_var_mat.push_back(res3);
    lst_selected_models.push_back(res2["mat_of_variables_selected_models"]);
    lst_index_selected_models.push_back(res2["id_selected_models"]);
    
    // print verbose 
    if (verbose) {
      Rcpp::Rcout << "Completed models of dimension " << dimension_model << "\n";
    }
  }
  // create return object
  List ret =  List::create(Named("lst_estimated_beta") = lst_estimated_beta,
                           Named("lst_AIC") = lst_AIC,
                           Named("lst_var_mat") = lst_var_mat,
                           Named("lst_selected_models") = lst_selected_models,
                           Named("lst_index_selected_models")= lst_index_selected_models
  );
  // output return object
  return(ret);
}



/*** R
# Parameters for data generation
# set.seed(1)
# n <- 2000
# p <- 50
# Sigma <- diag(rep(1/p, p))
# 
# X <- MASS::mvrnorm(n = n, mu = rep(0, p), Sigma = Sigma)
# beta = c(-10,5,6,19,70,rep(0,p-5))
# 
# # --------------------- logistic reg
# z <- 1 + X%*%beta
# pr <- 1/(1 + exp(-z))
# y <- as.factor(rbinom(n, 1, pr))
# y = as.numeric(y)-1

# 
# 
# 
# # 
# # 
# swag2 = function(X, y, p_max, family=binomial(), method=0, alpha=.3){
# 
#   # p_max = 5
#   # family=binomial()
#   # method=0
#   # alpha=.3
# 
# 
#   estimated_beta = list()
#   lst_AIC = list()
#   var_mat = list()
#   selected_models = list()
#   index_selected_models = list()
# 
#   res_screening = run_estimation_model_one_dimension_cpp(X, y, family, method);
# 
#   #  identify selected combination from first dimension
#   res2 = identify_selected_combinations_cpp(mat_of_variables = res_screening$matrix_of_variables, mat_criterion = res_screening$mat_AIC_dim_1, alpha = alpha);
#   
# 
#   
#   # store for dimension 1
#   estimated_beta[[1]] = res_screening$mat_beta_dim_1
#   lst_AIC[[1]] = res_screening$mat_AIC_dim_1
#   p=ncol(X)
#   var_mat[[1]] = matrix(0:(p-1), nrow=p)
#   selected_models[[1]] = res2$mat_of_variables_selected_models
#   index_selected_models[[1]] = res2$id_selected_models
# 
#   #  save variable screened at first dimension
#   variables_screening = res2$mat_of_variables_selected_models
#   # define m the maximum number of models which are screened at each dimension
#   p = dim(X)[2]
#   m = choose(n =  length(variables_screening), k =2)
#   for(dimension_model in 2:p_max){
# 
#     # dimension_model = 2
# 
#     # create matrix of row is a combination of variable and each value in the corresponding row is the index of a variables
#     res3 = compute_all_possible_variable_combinations_cpp(res2$mat_of_variables_selected_models, variables_screening)
#     
#     # select randomly m only row of this matrix if the number of row is greater than m
#     if(nrow(res3)>m){
#       res3 = res3[1:m, ]
#     }
#     
#     
#     
#     if(nrow(res3)==0){
#       break
#     }
# 
#     # run estimation on all possible models provided by res3
#     res4 = estimate_all_model_combinations_cpp(X, y, res3, family, method);
# 
#     # identify the sbest models in this dimension
#     res2 = identify_selected_combinations_cpp(res3,  res4$mat_AIC, alpha = alpha);
# 
#     # save at each dimension the selected model explored, their aic and etc
#     estimated_beta[[dimension_model]] = res4$mat_beta
#     lst_AIC[[dimension_model]] = res4$mat_AIC
#     var_mat[[dimension_model]] = res3
#     selected_models[[dimension_model]] = res2$mat_of_variables_selected_models
#     index_selected_models[[dimension_model]] = res2$id_selected_models
# 
#   }
# 
#   # create return
#   ret = list( "lst_estimated_beta" = estimated_beta,
#               "lst_AIC" = lst_AIC,
#               "lst_var_mat" = var_mat,
#               "lst_selected_models" = selected_models,
#               "lst_index_selected_models" = index_selected_models)
# 
#   return(ret)
# 
# 
# }
# 
# 
# 
# # 
# res1 = swaglm::swaglm(X = X, y = y, family = binomial(), method = 0,p_max = 15, alpha = .3, verbose = T)

# res
# res$lst_var_mat[2]

# res2 = swag2(X = X, y = y, p_max = 10, family = binomial(), method = 0, alpha = .4)
# 
# res$lst_var_mat
# res2$lst_var_mat
# # res2
# all.equal(res, res2)
# binomial_coefficient(1050,2)
# choose(1050,2)
# # for dimension 2
# alpha=.4
# res3 = compute_all_possible_variable_combinations_cpp(res2$mat_of_variables_selected_models, res$selected_models[[1]])
# res4 = estimate_all_model_combinations_cpp(X, y, res3, family=binomial(), method=0)
# res2 = identify_selected_combinations_cpp(res3, res4$mat_AIC, alpha = alpha)
# 
# # for dimension 3
# res3 = compute_all_possible_variable_combinations_cpp(res$selected_models[[1]], res$selected_models[[1]])





# 
# 
# res2  =var_mat
# res2-1 == res[[1]]
# res
# 
# 
# swag <- function(y, X, p_max, q, m = choose(floor(q*ncol(X)), 2), family = "gaussian", eval_func = AIC, seed = 123) {
#   
#   
#   p_max=10
#   family = binomial()
#   eval_func = AIC
#   seed = 123
#   q=0.4
#   m=choose(floor(q*ncol(X)), 2)
#   
#   if(p_max > ncol(X)) stop("p_max is larger than the number of predictors")
#   
#   p <- ncol(X)
#   index_screen <- 1:p
#   
#   # Initializing result storage
#   criteria <- list()
#   group <- list()
#   selected_group <- list()
#   
#   
#   ####################################
#   # Screening Step
#   ####################################
#   
#   crit <- rep(NA, p)
#   
#   for(i in seq_along(index_screen)) {
#     
#     # index of group of variables
#     fit <- glm(y ~ X[, i], family = family)
#     
#     # RMSE for each models
#     crit[i] = eval_func(fit)
#     
#   }
#   
#   criteria[[1]] <- crit
#   group[[1]] <- seq_along(crit)
#   id_screening <- selected_group[[1]] <- which(crit <= quantile(crit, q))
#   
#   
#   ####################################
#   # General Step
#   ####################################
#   
#   # if(progress == T) {
#   #   
#   #   pb <- txtProgressBar(min = 2, max = p_max, initial = 2, style = 3) 
#   #   
#   # }
#   
#   for(d in 2:p_max) {
#     
#     d=2
#     
#     id_row <- selected_group[[d - 1]] # indices of models in the prev. step with smaller error
#     
#     if(d == 2) {
#       
#       id_var <- group[[d - 1]][id_row] # group[[1]] is always a vector
#       nrv <- length(id_var)
#       
#     } else {
#       
#       if(length(id_row) == 1) { # only one model selected from previous dimension
#         
#         id_var <- as.matrix(t(group[[d - 1]][id_row, ]))
#         nrv <- nrow(id_var)
#         
#       } else {
#         
#         id_var <- group[[d - 1]][id_row,]
#         nrv <- nrow(id_var)
#         
#       }
#       
#     }
#     
#     # build all possible models
#     A <- matrix(nr = nrv*length(id_screening), nc = d)
#     A[, 1:(d - 1)] <- kronecker(cbind(rep(1, length(id_screening))), id_var)
#     A[, d] <- rep(id_screening, each = nrv)
#     B <- unique(t(apply(A, 1, sort))) # deletes the repeated rows
#     id_ndup <- which(apply(B, 1, anyDuplicated) == 0) # removes the models with same Xi
#     
#     if(length(id_ndup) == 1) {
#       
#       var_mat <- as.matrix(t(B[id_ndup, ]))
#       
#     } else {
#       
#       var_mat <- B[id_ndup, ] # all possible combinations of size d
#       
#     }
#     
#     rm(list=c("A", "B"))
#     
#     ##randomly selecting the models of size d
#     if(nrow(var_mat) > m) {
#       
#       set.seed(seed + d)
#       
#       group[[d]] <- var_mat[sample.int(nrow(var_mat), m), ]
#       
#     } else {
#       
#       group[[d]] <- var_mat
#       
#     }
#     
#     var_mat <- group[[d]]
#     
#     crit <- rep(NA, nrow(var_mat))
#     
#     ##training
#     for(i in seq_along(crit)){
#       
#       # index of group of variables
#       fit <- glm(y ~ X[, var_mat[i,]], family = family)
#       
#       # RMSE for each models
#       crit[i] = eval_func(fit)
#       
#     }
#     
#     criteria[[d]] <- na.omit(crit)
#     selected_group[[d]] <- which(crit <= quantile(crit, probs = q, na.rm = T))
#     
#     #if(progress == T) setTxtProgressBar(pb, d)
#     
#   }
#   
#   #close(pb)
#   
#   out <- list("group" = group, "selected_group" = selected_group, "criteria" = criteria, "id_screening" = id_screening, "p_max" = p_max)
#   class(out) <- "swag"
#   
#   return(out)
#   
# }
# 










# 
# 
# identify_selected_combinations = function(mat_of_variables, mat_criterion, alpha = 0.1){
#   # identify row that have smaller values at first column than quantile alpha
#   q_star = quantile(mat_criterion[,1], alpha)
#   id_selected_models = which(mat_criterion[,1] <= q_star)
#   ret = list(
#     "id_selected_models" = id_selected_models,
#     "mat_of_variables_selected_models" = mat_of_variables[id_selected_models, ]
#   )
#   return(ret)
# }

# identify_selected_combinations(mat_of_variables = res$matrix_of_variables, mat_criterion = res$mat_AIC_dim_1, alpha = .1)
# identify_selected_combinations_cpp(mat_of_variables = res$matrix_of_variables, mat_criterion = res$mat_AIC_dim_1, alpha = .2)
# 
# 
# identify_selected_combinations(res2$)


*/

