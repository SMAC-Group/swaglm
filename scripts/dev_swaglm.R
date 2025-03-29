# --------------------------------------------
# --- dev swaglm


rm(list=ls())

# library(fastglm)
library(swaglm)


# Parameters for data generation
set.seed(1)
n <- 2000
p <- 20
Sigma <- diag(rep(1/p, p))

X <- MASS::mvrnorm(n = n, mu = rep(0, p), Sigma = Sigma)
beta = c(-10,5,6,19,70,rep(0,15))

# --------------------- logistic reg
z <- 1 + X%*%beta
pr <- 1/(1 + exp(-z))
y <- as.factor(rbinom(n, 1, pr))
y = as.numeric(y)-1


# demo fit
fit=glm(y~X+1, family = stats::binomial())
fit


# ----------------------------- gaussian iid regression
# 
# sigma2 = 5
# y = 1 + X %*% beta + rnorm(n, sd = sqrt(sigma2))
# # run_estimation_model_one_dimension_cpp(X, y = y, family = gaussian())




# -------------- define R function
# an integer scalar with value 0 for the column-pivoted QR decomposition, 1 for the unpivoted QR decomposition, 2 for the LLT Cholesky, or 3 for the LDLT Cholesky
# run_estimation_model_one_dimension = function(X, y, family=binomial(), method = 0){
#   
#  # X=X_mat_predictors
#   p = dim(X)[2]
#   n = dim(X)[1]
#   
#   # Initialize beta matrices to store beta coef and AIC values
#   mat_AIC_dim_1 = matrix(nrow=p, ncol=1)
#   mat_beta_dim_1 = matrix(nrow=p, ncol=2)
#   
#   for(i in 1:p){
#     X_mat_i = cbind(rep(1, n), X[,i])
#     fit = fastglm::fastglm(x =X_mat_i, y=y, family=family, method = method)
#     mat_AIC_dim_1[i,1] = fit$aic
#     mat_beta_dim_1[i, ] = fit$coefficients
#   }
#   ret=list("mat_AIC_dim_1" = mat_AIC_dim_1,
#            "mat_beta_dim_1" = mat_beta_dim_1)
#   return(
#     ret
#   )
# }

# res1 = run_estimation_model_one_dimension(X=X, y = y)
# res2 = run_estimation_model_one_dimension_cpp(X, y = y)
# res1
# res2
# all.equal(res1, res2)

# -----------benchmark
# microbenchmark::microbenchmark(res1 = run_estimation_model_one_dimension(X=X, y = y),
#                                res2 = run_estimation_model_one_dimension_cpp(X=X, y = y), times = 10)



#-----------------------------------------------------------------



swag <- function(y, X, p_max, q, m = choose(floor(q*ncol(X)), 2), family = "gaussian", eval_func = AIC, seed = 123) {
  
  # p_max=10
  # family = binomial()
  # eval_func = AIC
  # seed = 123
  # q=0.4
  # m=choose(floor(q*ncol(X)), 2)
  
  if(p_max > ncol(X)) stop("p_max is larger than the number of predictors")
  
  p <- ncol(X)
  index_screen <- 1:p
  
  # Initializing result storage
  criteria <- list()
  group <- list()
  selected_group <- list()
  
  
  ####################################
  # Screening Step
  ####################################
  
  crit <- rep(NA, p)
  
  
  # compute criterion for all single variable models
  for(i in seq_along(index_screen)) {
    
    # index of group of variables
    fit <- glm(y ~ X[, i], family = family)
    
    # RMSE for each models
    crit[i] = eval_func(fit)
    
  }
  
  criteria[[1]] <- crit
  group[[1]] <- seq_along(crit)
  id_screening <- selected_group[[1]] <- which(crit <= quantile(crit, q))
  
  
  ####################################
  # General Step
  ####################################
  
  # if(progress == T) {
  #   
  #   pb <- txtProgressBar(min = 2, max = p_max, initial = 2, style = 3) 
  #   
  # }
  
  for(d in 2:p_max) {
    
    # d=2
    
    id_row <- selected_group[[d - 1]] # indices of models in the prev. step with smaller error
    
    if(d == 2) {
      
      id_var <- group[[d - 1]][id_row] # group[[1]] is always a vector
      nrv <- length(id_var)
      
    } else {
      
      if(length(id_row) == 1) { # only one model selected from previous dimension
        
        id_var <- as.matrix(t(group[[d - 1]][id_row, ]))
        nrv <- nrow(id_var)
        
      } else {
        
        id_var <- group[[d - 1]][id_row,]
        nrv <- nrow(id_var)
        
      }
      
    }
    
    # build all possible models
    
    # create matrix with number of row corresponding to number of selected model from previous dimension plus always all the screened variables from first step
    A <- matrix(nr = nrv*length(id_screening), nc = d)
    # create all first 1 to 2 columns the extracted model so basically replicating the matrix id_var
    A[, 1:(d - 1)] <- kronecker(cbind(rep(1, length(id_screening))), id_var)
    # fill in last column with repeating nrv times each values in id screening
    A[, d] <- rep(id_screening, each = nrv)
    B <- unique(t(apply(A, 1, sort))) # deletes the repeated rows
    id_ndup <- which(apply(B, 1, anyDuplicated) == 0) # removes the models with same Xi
    
    if(length(id_ndup) == 1) {
      
      var_mat <- as.matrix(t(B[id_ndup, ]))
      
    } else {
      
      var_mat <- B[id_ndup, ] # all possible combinations of size d
      
    }
    
    rm(list=c("A", "B"))
    
    ##randomly selecting the models of size d
    if(nrow(var_mat) > m) {
      
      set.seed(seed + d)
      
      group[[d]] <- var_mat[sample.int(nrow(var_mat), m), ]
      
    } else {
      
      group[[d]] <- var_mat
      
    }
    
    var_mat <- group[[d]]
    
    crit <- rep(NA, nrow(var_mat))
    
    ##training
    for(i in seq_along(crit)){
      
      # index of group of variables
      fit <- glm(y ~ X[, var_mat[i,]], family = family)
      
      # RMSE for each models
      crit[i] = eval_func(fit)
      
    }
    
    criteria[[d]] <- na.omit(crit)
    selected_group[[d]] <- which(crit <= quantile(crit, probs = q, na.rm = T))
    
    #if(progress == T) setTxtProgressBar(pb, d)
    
  }
  
  #close(pb)
  
  out <- list("group" = group, "selected_group" = selected_group, "criteria" = criteria, "id_screening" = id_screening, "p_max" = p_max)
  class(out) <- "swag"
  
  return(out)
  
}

res = swaglm::swaglm(X=X, y = y, p_max = 10, family = binomial(), method = 0, alpha = .4, verbose = T)
res
# 
# microbenchmark::microbenchmark(res1 = swag(X=X, y = y, p_max = 10, q = .4, family = binomial(), eval_func = AIC),
#                                res2 = swaglm::swaglm(X=X, y = y, p_max = 10, family = binomial(), method = 0, alpha = .4, verbose = F), times=10)
# 
# 
# 
# res1$group
# res1$selected_group
