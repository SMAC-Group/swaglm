# --------------------------------------------
# --- dev swaglm


rm(list=ls())

# library(fastglm)
library(swaglm)


# Parameters for data generation
set.seed(12345)
n <- 2000
p <- 100

# create design matrix and vector of coefficients
Sigma <- diag(rep(1/p, p))
X <- MASS::mvrnorm(n = n, mu = rep(0, p), Sigma = Sigma)
beta = c(-15,-10,5,10,15, rep(0,p-5))

# --------------------- generate from logistic regression with an intercept of one
z <- 1 + X%*%beta
pr <- 1/(1 + exp(-z))
y <- as.factor(rbinom(n, 1, pr))
y = as.numeric(y)-1

# define swag parameters
quantile_alpha = .1
p_max = 10
swag_obj = swaglm::swaglm(X=X, y = y, p_max = p_max, family = stats::binomial(), alpha = quantile_alpha, verbose = T, seed = 123)
names(swag_obj)
swag_network = compute_network(swag_obj)
plot(swag_network, scale_vertex = .4)
# # demo fit
# fit=glm(y~X+1, family = stats::binomial())
# fit


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

#

quantile_alpha = .4
p_max  =20
res1 = swag(y = y, X = X, p_max = p_max, q = quantile_alpha, family = binomial(), eval_func = AIC)
res2 = swaglm::swaglm(X=X, y = y, p_max = p_max, family = binomial(),alpha = quantile_alpha, verbose = T)
# res1$group
# res2$lst_var_mat
# microbenchmark::microbenchmark(res1 = swag(y = y, X = X, p_max = p_max, q = quantile_alpha, family = binomial(), eval_func = AIC),
#                                res2 = swaglm::swaglm(X=X, y = y, p_max = p_max, family = binomial(),
#                                                      alpha = quantile_alpha, verbose = F), times=10)







# test to make sure the new implementation of swag works
res2$lst_estimated_beta
# lets take model of 5 dimension first model
res2$lst_var_mat[[5]][1,]
Xtemp = cbind(rep(1,n), X[, (res2$lst_var_mat[[5]][1,] +1)])
fit = glm(y~Xtemp-1, family = binomial())
all.equal(fit$coefficients, res2$lst_estimated_beta[[5]][1,])
res2$lst_AIC[[5]][1,]
AIC(fit)

# check also for the selected variables
res2$lst_var_mat[[4]][res2$lst_index_selected_models[[4]] +1, ]  == 
res2$lst_selected_models[[4]]

# good, now lets add the function compute network



swag_network <- function(obj, mode = "undirected", weighted = F, show_net = T) {
  
  
  obj = res2
  mode = "undirected"
  weighted = F
  show_net = T
  
  
  
  
  p_max=length(obj$lst_AIC)
  # 
  require(plyr)
  require(igraph)
  require(gdata)
  
  selected_ind <- list()
  
  selected_ind[[1]] <- obj$lst_var_mat[[1]][obj$lst_index_selected_models[[1]]+1, ]
  ## add plus one because variables in R
  selected_ind[[1]] = selected_ind[[1]]+1
  
  for(i in 2:p_max) {

    
    if(length(obj$lst_index_selected_models[[i]]) != 1) {
      
      selected_ind[[i]] <- obj$lst_var_mat[[i]][as.matrix(obj$lst_index_selected_models[[i]])+1, ]
      selected_ind[[i]] = selected_ind[[i]] +1
    } else {
      selected_ind[[i]] = obj$lst_var_mat[[i]][obj$lst_index_selected_models[[i]]+1, ]
      selected_ind[[i]] = t(selected_ind[[i]] +1)
    }
    
  }
  
  
  models <- rbind.fill.matrix(selected_ind)
  
  models2 = rbind.fill.matrix(obj$lst_selected_models) +1
  all.equal(models, models2)
  #### intensity matrix
  
  selected_var <- c()
  
  for(i in 1:ncol(models)) {
    
    selected_var <- c(selected_var, models[, i])
    
  }
  
  selected_var <- sort(na.omit(unique(selected_var)))
  
  A <- matrix(0, nrow = ncol(models), ncol = ncol(models))
  intensity <- matrix(0, nrow = length(selected_var), ncol = length(selected_var))
  a <- list()
  
  for(i in 1:(length(selected_var) - 1)) {
    
    for(j in (i + 1):length(selected_var)) {
      
      for(k in 1:(ncol(models) - 1)) {
        
        a[[i]] <- which(models[, k] == selected_var[i])
        
        for(n in (k + 1):(ncol(models))) {
          
          A[k, n] <- length(which(models[a[[i]], n] == selected_var[j]))
          
        }
        
      }
      
      intensity[j, i] <- intensity[i, j] <- sum(A)
      
    }
    
  }
  
  colnames(intensity) <- obj$lst_selected_models[[1]]+1
  rownames(intensity) <- obj$lst_selected_models[[1]]+1
  
  #relation matrix for pairwise connection
  
  relation_mat = matrix(NA, nrow = choose(length(selected_var),2), ncol = 3)
  c = rep(selected_var[1],length(selected_var)-1)
  for(i in (length(selected_var)-2):1){
    a = rep(selected_var[length(selected_var)-i],i)
    c = c(c,a)
  }
  relation_mat[,1]=c
  
  bb = selected_var[-1]
  for(i in 3:length(selected_var)){
    b1 = selected_var[i:length(selected_var)]
    bb = c(bb,b1)
  }
  relation_mat[,2]=bb
  
  relation_mat[,3] = upperTriangle(intensity,byrow = T)
  
  g <- graph_from_adjacency_matrix(intensity, mode = mode, weighted = NULL)
  
  vertex_degrees_obs <- degree(g)
  
  E(g)$weight <- 1  
  
  # Simplify the graph, summing the weights of multiple edges
  g_simplified_obs <- simplify(g, edge.attr.comb = list(weight = "sum"))
  
  # Set the edge width based on the combined weights
  E(g_simplified_obs)$width <- E(g_simplified_obs)$weight
  
  
  if(show_net == T) plot(g_simplified_obs, layout = layout.circle, vertex.color = "skyblue", edge.color = "black", vertex.size = 0.1*degree(g), edge.width = 0.1*E(g_simplified_obs)$width )
  
  return(list(g = g, models = models))
  
}











# 
# 
# res1$group
# res1$selected_group
