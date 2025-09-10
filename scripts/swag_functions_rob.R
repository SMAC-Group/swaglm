# Ideally we would like to get the distribution of betas from this function
swag <- function(y, X, p_max, q, m = choose(floor(q*ncol(X)), 2), family = "gaussian", eval_func = AIC, seed = 123, progress = T, include_var = NULL) {
  
  if(p_max > ncol(X)) stop("p_max is larger than the number of predictors")
  
  p <- ncol(X)
  index_screen <- 1:p
  
  # Initializing result storage
  criteria <- list()
  group <- list()
  selected_group <- list()
  betas <- list()
  pvals <- list()
  
  
  ####################################
  # Screening Step
  ####################################
  
  crit <- rep(NA, p)
  beta_mat <- t(as.matrix(rep(NA, ncol(X))))
  pval_mat <- t(as.matrix(rep(NA, ncol(X))))
  colnames(beta_mat) <- colnames(pval_mat) <- colnames(X)
  
  for(i in seq_along(index_screen)) {
    
    X_subset <- as.data.frame(X[, i, drop = FALSE])
    colnames(X_subset) <- paste0("V", i)
    
    # Combine with y
    dat <- data.frame(y = y, X_subset)
    
    # Create formula dynamically
    formula <- as.formula(paste("y ~", paste(colnames(X_subset), collapse = " + ")))
    
    # index of group of variables
    fit <- glm(formula, data = dat, family = binomial())
    
    # RMSE for each models
    crit[i] = eval_func(fit)
    
    # Update beta matrix
    beta <- rep(NA, ncol(X))
    beta[i] <- fit$coefficients[-1]
    beta_mat <- rbind(beta_mat, beta)
    
    # Update pval matrix
    pval <- rep(NA, ncol(X))
    pval[i] <- coef(summary(fit))[-1 , "Pr(>|z|)"]
    pval_mat <- rbind(pval_mat, pval)
    
  }
  
  criteria[[1]] <- crit
  group[[1]] <- seq_along(crit)
  betas[[1]] <- beta_mat[-1, ]
  pvals[[1]] <- pval_mat[-1, ]
  id_screening <- selected_group[[1]] <- which(crit <= quantile(crit, q))
  
  if(!is.null(include_var)) id_screening <- selected_group[[1]] <- unique(c(include_var, id_screening))
  
  
  ####################################
  # General Step
  ####################################
  
  if(progress == T) {
    
    pb <- txtProgressBar(min = 2, max = p_max, initial = 2, style = 3) 
    
  }
  
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
    beta_mat <- t(as.matrix(rep(NA, ncol(X))))
    pval_mat <- t(as.matrix(rep(NA, ncol(X))))
    colnames(beta_mat) <- colnames(pval_mat) <- colnames(X)
    
    ##training
    for(i in seq_along(crit)){
      
      X_subset <- as.data.frame(X[, var_mat[i,], drop = FALSE])
      colnames(X_subset) <- paste0("V", var_mat[i,])
      
      # Combine with y
      dat <- data.frame(y = y, X_subset)
      
      # Create formula dynamically
      formula <- as.formula(paste("y ~", paste(colnames(X_subset), collapse = " + ")))
      
      # index of group of variables
      fit <- glm(formula, data = dat, family = binomial())
      
      # RMSE for each models
      crit[i] = eval_func(fit)
      
      # Update beta matrix
      beta <- rep(NA, ncol(X))
      beta[var_mat[i,]] <- fit$coefficients[-1]
      beta_mat <- rbind(beta_mat, beta)
      
      # Update pval matrix
      pval <- rep(NA, ncol(X))
      pval[var_mat[i,]] <- coef(summary(fit))[-1, "Pr(>|z|)"]
      pval_mat <- rbind(pval_mat, pval)
      
    }
    
    criteria[[d]] <- na.omit(crit)
    selected_group[[d]] <- which(crit <= quantile(crit, probs = q, na.rm = T))
    betas[[d]] <- beta_mat[-1, ]
    pvals[[d]] <- pval_mat[-1, ]
    
    if(progress == T) setTxtProgressBar(pb, d)
    
  }
  
  if(progress == T) close(pb)
  
  out <- list("group" = group, "selected_group" = selected_group, "criteria" = criteria, "id_screening" = id_screening, "p_max" = p_max, "family" = family, "betas" = betas, "p_values" = pvals)
  class(out) <- "swag"
  
  return(out)
  
}

# ------------------ test swag fct 
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

test1 = swag(y = y, X = X, p_max = 15, q = .1, m = 100, family = binomial())
test1$group
test1$selected_group
test1$criteria
test1$id_screening
test1$selected_group
test2 = swaglm::swaglm(y = y, X = X, p_max = 15, alpha = .1, family = binomial())


  
# this function just add the argument best_model to the swag object... 
swag_rashomon <- function(obj, bound = NA, delta = 0.01) {
  
  # obj = test1 
  # delta = 0.01
  
  if(is.na(bound)) {
    
    q_d <- quantile(obj$criteria[[which.min(unlist(lapply(obj$criteria, median)))]], probs = delta, na.rm = T)
    
  } else {
    
    q_d <- bound
    
  }
  
  best_model <- c(1, which.min(obj$criteria[[1]]), min(obj$criteria[[1]]))
  
  for(i in seq_along(obj$p_max)) {
    
    obj$selected_group[[i]] <- which(obj$criteria[[i]] <= q_d)
    
    if(best_model[3] > min(obj$criteria[[i]])) best_model <- c(i, which.min(obj$criteria[[i]]), obj$criteria[[i]])
    
    
  }
  
  obj$best_model <- best_model
  
  return(obj)
  
}



# swag_rashomon(test1)


get_betas <- function(obj) {
  
  beta_mat <- rep(NA, ncol(obj$betas[[1]]))
  pval_mat <- rep(NA, ncol(obj$p_values[[1]]))
  
  for(i in 1:length(obj$betas)) {
    
    beta <- obj$betas[[i]][obj$selected_group[[i]], ]
    beta_mat <- rbind(beta_mat, beta)
    
    pval <- obj$p_values[[i]][obj$selected_group[[i]], ]
    pval_mat <- rbind(pval_mat, pval)
    
  }
  
  return(list("betas" = beta_mat[-1, ], "p_vals" = pval_mat[-1, ]))
  
}


get_betas(test1)

summary.swag <- function(obj) {
  
  # number of models (total + per-dimension)
  
  # range of model sizes
  
  # frequency of variables
  
  # distribution of betas for each variable
  
  # summary of errors/criteria
  
}

plot.swag <- function(obj) {
  
  require(plotrix)
  
  m_vector <- sapply(obj$criteria, function(x) summary(x)[4])
  l_vector <- sapply(obj$criteria, function(x) summary(x)[1])
  u_vector <- sapply(obj$criteria, function(x) summary(x)[6])
  
  plotCI(1:length(obj$criteria), m_vector, ui = u_vector, li = l_vector, scol = "grey", col="red", pch = 16, main = "Ranges of criterion", ylab = "Range", xlab = "Model Size")
  
  
}

swag_network <- function(obj, mode = "undirected", weighted = NULL, show_net = T) {
  
  require(plyr)
  require(igraph)
  
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
  
  g <- graph_from_adjacency_matrix(intensity, mode = mode, weighted = weighted)
  
  if(show_net == T) plot(g)
  
  return(g)
  
}

swag_edge <- function(obj, sort_down = T) {
  
  require(plyr)
  
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
  
  E <- matrix(NA, ncol = 3, nrow = choose(length(obj$id_screening), 2))
  
  E[, 1:2] <- t(combn(obj$id_screening, 2))
  E[,3] <- 0
  
  for(j in 1:nrow(E)) {
    
    for(i in 1:nrow(models)) {
      
      if(sum(E[j, 1:2] %in% models[i, ]) == 2) {
        
        E[j, 3] <- E[j, 3] + 1
        
      }
      
    }
    
  }
  
  E <- E[order(E[, 3], decreasing = sort_down), ]
  
  return(E)
  
}

swag_freq <- function(obj, sort_down = T) {
  
  require(plyr)
  
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
  
  freq_table <- sort(table(models), decreasing = sort_down)
  
  return(freq_table)
  
}

pred_swag <- function(obj, y_train, x_train, y_test, x_test, threshold = 0.5) {
  
  require(plyr)
  
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
  
  pred <- matrix(NA, length(y_test), nrow(models))
  
  for(i in seq_len(nrow(models))) {
    
    index <- na.omit(models[i, ])
    
    # index of group of variables
    train_data <- data.frame(y_train = y_train, x_train[, index])
    fit <- glm(y_train ~ ., family = obj$family, data = train_data)
    
    new_data <- as.data.frame(x_test[, index])
    names(new_data) <- names(train_data)[-1]
    
    # Predictions
    if(obj$family == "binomial") {
      
      pred[, i] = predict(fit, newdata = new_data, type = "response") > threshold
      
    } else {
      
      pred[, i] = predict(fit, newdata = new_data, type = "response")
      
    }
    
  }
  
  if(obj$family == "binomial") {
    
    class_error <- function(y_hat, y) 1 - sum(diag(table(y, y_hat)))/length(y)
    
    pred <- apply(pred, c(1, 2), function(x) {
      factor(x, levels = c(FALSE, TRUE), labels = levels(y_test))
    })
    error <- apply(pred, 2, class_error, y = y_test)
    
  } else {
    
    l2_error <- function(y_hat, y) sum((y - y_hat)^2)
    
    error <- apply(pred, 2, l2_error, y = y_test)
    
  }
  
  return(list("predictions" = pred, "error" = error))
  
}