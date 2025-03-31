library(DescTools)

swag <- function(y, X, p_max, q, m = choose(floor(q*ncol(X)), 2), family = "gaussian", eval_func = AIC, seed = 123) {
  
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
    A <- matrix(nr = nrv*length(id_screening), nc = d)
    A[, 1:(d - 1)] <- kronecker(cbind(rep(1, length(id_screening))), id_var)
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

swag_network <- function(obj, mode = "undirected", weighted = F, show_net = T) {
  
  
  # obj = swag_obs
  # mode = "undirected"
  # weighted = F
  # show_net = T
  # 
  require(plyr)
  require(igraph)
  require(gdata)
  
  selected_ind <- list()
  
  selected_ind[[1]] <- obj$group[[1]][obj$selected_group[[1]]]
  
  for(i in 2:obj$p_max) {
    
    if(length(obj$selected_group[[i]]) != 1) {
      
      selected_ind[[i]] <- obj$group[[i]][as.matrix(obj$selected_group[[i]]), ]
      
    } else {
      
      selected_ind[[i]] <- t(as.matrix(obj$group[[i]][obj$selected_group[[i]], ]))
      
    }
    
  }
  
  
  models <- rbind.fill.matrix(selected_ind)
  
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
  
  colnames(intensity) <- obj$id_screening
  rownames(intensity) <- obj$id_screening
  
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





# Function to calculate optimal bandwidth using Silverman's rule of thumb
optimal_bandwidth <- function(data) {
  
  n <- length(data)
  s <-   sd(data)                # Standard deviation of data
  IQR <- IQR(data)               # Interquartile range of data
  h <- 0.9 * min(s, IQR / 1.34) * n^(-1/5)
  
  return(h)
  
}

# Smoothed bootstrap function with bandwidth
smoothed_bootstrap <- function(data, H = 1000, bandwidth = NULL) {
  
  n <- length(data)
  
  # If bandwidth is not provided, use Silverman's rule to select it
  if (is.null(bandwidth)) {
    bandwidth <- optimal_bandwidth(data)
  }
  
  # Matrix to store bootstrap samples
  boot_samples <- matrix(NA, nrow = H, ncol = n)
  
  # Perform B bootstrap replicates
  for (i in 1:H) {
    
    # Resample the data (with replacement)
    resample <- sample(data, size = n, replace = TRUE)
    
    # Add smoothed noise with standard deviation based on bandwidth
    noise <- rnorm(n, mean = 0, sd = bandwidth)
    
    # Smoothed bootstrap sample
    boot_samples[i, ] <- resample + noise
    
  }
  
  return(boot_samples)
  
}



swag_significance_test <- function(y, X, swag_obs, net_obs, p_max, q, sign_level = 0.05, B = 100) {
  
  library(DescTools)
  
  # Frequency table for observed data
  frequency <- table(net_obs$models)
  variable <- swag_obs$id_screening
  freq_obs <- cbind(variable, frequency)
  freq_obs <- freq_obs[order(-freq_obs[, "frequency"]), ]
  
  # Compute observed statistics
  entropy_freq_obs <- Entropy(freq_obs)
  entropy_eigen_obs <- Entropy(eigen_centrality(net_obs$g)$vector)
  
  # Initialize vectors for null statistics
  entropy_freq_null <- numeric(B)
  entropy_eigen_null <- numeric(B)
  
  for (b in 1:B) {
    set.seed(b)
    # Generate response under null
    y_null <- sample(y, length(y), replace = TRUE) 
    
    # Run SWAG under null
    swag_null <- swag(y_null, X, p_max = p_max, q = q, family = "binomial")
    net_null <- swag_network(swag_null)
    net_null$g <- delete_vertices(net_null$g, V(net_null$g)[degree(net_null$g) == 0])
    
    # Frequency table under null
    frequency <- table(net_null$models)
    variable <- swag_null$id_screening
    freq_null <- cbind(variable, frequency)
    freq_null <- freq_null[order(-freq_null[, "frequency"]), ]
    
    # Compute null statistics
    entropy_freq_null[b] <- Entropy(freq_null)
    entropy_eigen_null[b] <- Entropy(eigen_centrality(net_null$g)$vector)
  }
  
  # smoothed bootstrap
  
  entropy_freq_null = as.vector(smoothed_bootstrap(entropy_freq_null))
  entropy_eigen_null = as.vector(smoothed_bootstrap(entropy_eigen_null))
  
  # Compute p-values
  p_value_eigen <- mean(entropy_eigen_null < entropy_eigen_obs)
  p_value_freq <- mean(entropy_freq_null < entropy_freq_obs)
  
  # Determine significance
  eigen_significance <- p_value_eigen < sign_level
  freq_significance <- p_value_freq < sign_level
  
  # Return results
  return(list(
    eigenvector_centrality = list(
      result = ifelse(eigen_significance, "SWAG is significant by eigenvector centrality", "SWAG is not significant by eigenvector centrality"),
      p_value = p_value_eigen
    ),
    frequency = list(
      result = ifelse(freq_significance, "SWAG is significant by frequency", "SWAG is not significant by frequency"),
      p_value = p_value_freq
    )
  ))
}

### Example

# Parameters for data generation
set.seed(1)
n <- 2000
p <- 20
Sigma <- diag(rep(1/p, p))

X <- MASS::mvrnorm(n = n, mu = rep(0, p), Sigma = Sigma)
beta = c(-10,5,6,19,70,rep(0,15))
z <- 1 + X%*%beta
pr <- 1/(1 + exp(-z))
y <- as.factor(rbinom(n, 1, pr))



# run procedure
swag_obs = swag(y, X, p_max=5, q = 0.7, family = "binomial")
swag_obs$group[[1]]
swag_obs$group[[2]]



swag_net = swag_network(swag_obs)
swag_test = swag_significance_test(y,X, swag_obs = swag_obs, net_obs = swag_net,  p_max=5, q = 0.7, sign_level = 0.1, B=10)
