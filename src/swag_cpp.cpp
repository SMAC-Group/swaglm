// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;





// ------------------------------------------------ estimation model one dimension / screening


// This function runs an estimation model in one dimension using the fastglm function
// [[Rcpp::export]]
Rcpp::List run_estimation_model_one_dimension_cpp(const arma::mat& X, const arma::vec& y,  Nullable<List> family = R_NilValue, int method = 0) {
  int p = X.n_cols;
  int n = X.n_rows;
  
  // Initialize matrices to store AIC and beta coefficients
  arma::mat mat_AIC_dim_1(p, 1);
  arma::mat mat_beta_dim_1(p, 2);
  
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
  
  // Create the design matrix with an intercept column, to avoid recomputing it every time
  arma::mat X_with_intercept(n, p + 1, arma::fill::ones);  // Add a column of ones for the intercept
  X_with_intercept.cols(1, p) = X;  // Fill the rest with the columns from X
  
  for (unsigned int i = 0; i < p; ++i) {
    // Create the design matrix with an intercept column
    // arma::mat X_mat_i = join_rows(arma::ones<arma::vec>(n), X.col(i));
    // Select only the intercept and the ith column
    arma::uvec idx = {0, i + 1};  // Create an index for the first and i-th column
    arma::mat X_mat_i = X_with_intercept.cols(idx); 
    
    // Fit the model using fastglm
    // This is not optimal, fastglm is a R function implemented in Rcpp Eigen, but still we call its R implementation
    List fit = fastglm(Named("x") = X_mat_i, Named("y") = y, Named("family") = family_obj, Named("method") = method);
    
    // Extract AIC and coefficients
    mat_AIC_dim_1(i, 0) = as<double>(fit["aic"]);
    mat_beta_dim_1(i, 0) = as<arma::vec>(fit["coefficients"])(0);
    mat_beta_dim_1(i, 1) = as<arma::vec>(fit["coefficients"])(1);
  }
  
  // create matrix of variables for that dimension
  arma::mat matrix_variables(p,1, arma::fill::zeros);
  
  // Fill the first column with values from 0 to n-1
  matrix_variables.col(0) = arma::linspace<arma::vec>(0, p-1, p);
  
  return List::create(Named("mat_AIC_dim_1") = mat_AIC_dim_1,
                      Named("mat_beta_dim_1") = mat_beta_dim_1,
                      Named("matrix_of_variables") = matrix_variables);
}

 
 
 
// ------------------------------------------------------- identify selected models
 
 
 
 
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
List identify_selected_combinations_cpp(const arma::mat& mat_of_variables, const arma::mat& mat_criterion, double alpha = 0.1) {
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



// ---------------------------------------------------------------- estimate all models combinations



// Function to estimate all model combinations
// [[Rcpp::export]]
List estimate_all_model_combinations_cpp(const arma::mat& X, const arma::vec& y, const arma::mat& matrix_of_variables,
                                         Nullable<List> family = R_NilValue, int method = 0) {
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



// ------------------- create matrix of model

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



// ---------------------------------- SWAG algorithm




//' sawglm 
//'
//' This function runs the SWAG algorithm.
//' @name sawglm
//' @param X A numeric matrix of predictors.
//' @param y A numeric vector of responses.
//' @param family A glm family object (default is binomial).
//' @param p_max maximum dimension to explore
//' @param method an integer scalar with value 0 for the column-pivoted QR decomposition, 1 for the unpivoted QR decomposition, 2 for the LLT Cholesky, or 3 for the LDLT Cholesky
//' @param alpha lower quantile to of criterion
//' @param verbose boolean for verbose
//' @param seed A \code{int} that is the random seed used when creating the set of model to explore for the next dimension
//' @return A list of list of matrices. Each list correspond to the estimated coefficients from the estimated models, xxx.
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
  List res2 = identify_selected_combinations_cpp(as<arma::mat>(res_screening["matrix_of_variables"]), 
                                                 as<arma::mat>(res_screening["mat_AIC_dim_1"]), alpha = alpha);
  
  
  
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
    List new_res2 = identify_selected_combinations_cpp(res3,  as<arma::mat>(res4["mat_AIC"]), alpha = alpha);
    
    // update res 2
    res2 = new_res2;
    
    // save
    lst_estimated_beta.push_back( res4["mat_beta"]);
    lst_AIC.push_back(res4["mat_AIC"]);
    lst_var_mat.push_back(res3);
    lst_selected_models.push_back(res2["mat_of_variables_selected_models"]);
    lst_index_selected_models.push_back(res2["id_selected_models"]);
    
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

