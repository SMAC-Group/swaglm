# library(swaglm)
# n <- 2000
# p <- 50
#
# # create design matrix and vector of coefficients
# Sigma <- diag(rep(1/p, p))
# X <- MASS::mvrnorm(n = n, mu = rep(0, p), Sigma = Sigma)
# beta = c(-15,-10,5,10,15, rep(0,p-5))
#
# # --------------------- generate from logistic regression with an intercept of one
# z <- 1 + X%*%beta
# pr <- 1/(1 + exp(-z))
# y <- as.factor(rbinom(n, 1, pr))
# y = as.numeric(y)-1
#
# # define swag parameters
# quantile_alpha = .15
# p_max = 20
# swag_obj = swaglm::swaglm(X=X, y = y, p_max = p_max, family = stats::binomial(),
#                           alpha = quantile_alpha, verbose = TRUE, seed = 123)
# swag_network = compute_network(swag_obj)
# plot(swag_network, scale_vertex = .05)



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


#' swag_test
#'
#' Compute significance of identified set of variables
#'
#' @param swag_obj An object of class \code{swaglm}.
#' @importFrom igraph  eigen_centrality V delete_vertices degree
#' @importFrom DescTools Entropy
#' @importFrom progress progress_bar
#' @importFrom stats rnorm sd
#' @param significance_level A \code{double} between 0 and 1 indicating the specified significance level. Default is 0.05.
#' @param B a \code{integer} specifying the number of swag procedures to generate a distribution of the network statistics under the null.
#' @param verbose A \code{boolean} used to control verbose
#' @export
#' @examples
#' n <- 2000
#' p <- 50
#' 
#' # create design matrix and vector of coefficients
#' Sigma <- diag(rep(1/p, p))
#' X <- MASS::mvrnorm(n = n, mu = rep(0, p), Sigma = Sigma)
#' beta = c(-15,-10,5,10,15, rep(0,p-5))
#' z <- 1 + X%*%beta
#' pr <- 1/(1 + exp(-z))
#' y <- as.factor(rbinom(n, 1, pr))
#' y = as.numeric(y)-1
#' quantile_alpha = .15
#' p_max = 20
#' swag_obj = swaglm::swaglm(X=X, y = y, p_max = p_max, family = stats::binomial(),
#'                           alpha = quantile_alpha, verbose = TRUE, seed = 123)
#' swag_test(swag_obj, significance_level = .05, B = 50, verbose = TRUE)
#' 
swag_test <- function(swag_obj, significance_level = 0.05, B = 50, verbose = FALSE) {

  # ------------------------------ extract parameters from swag object to provide them later 
  y = swag_obj$y
  X = swag_obj$X
  p_max = swag_obj$p_max
  family= swag_obj$family
  alpha = swag_obj$alpha
  method = swag_obj$method
  
  # ----------------------- compute network on swag obj
  net_obj = compute_network(swag_obj)
  
  # Frequency table for observed data
  models =  net_obj$models
  frequency <- table(models)
  variable <- swag_obj$lst_selected_models[[1]]+1
  freq_obs <- cbind(variable, frequency)
  freq_obs <- freq_obs[order(-freq_obs[, "frequency"]), ]

  # Compute observed statistics on network
  entropy_freq_obs <- DescTools::Entropy(freq_obs)
  entropy_eigen_obs <- DescTools::Entropy(igraph::eigen_centrality(net_obj$g)$vector)

  # Initialize vectors for null statistics
  entropy_freq_null <- numeric(B)
  entropy_eigen_null <- numeric(B)
  
  # start progress bar
  if(verbose){
    pb <- progress_bar$new(total = B)
  }
  
  # start bootstrap procedure
  for (b in 1:B) {
    seed_b = 123 + b

    # Generate response under null by resampling y 
    y_null <- sample(y, length(y), replace = TRUE)

    # Run SWAG under null
    swag_null <- swaglm(y = y_null, X =  X, p_max = p_max, alpha =  alpha, 
                        family = family, method = method, seed = seed_b, verbose=FALSE)
    net_null <- compute_network(swag_null)
    net_null$g <- igraph::delete_vertices(net_null$g, igraph::V(net_null$g)[igraph::degree(net_null$g) == 0])

    # Frequency table under null
    frequency <- table(net_null$models)
    variable <- swag_null$lst_selected_models[[1]]+1
    freq_null <- cbind(variable, frequency)
    freq_null <- freq_null[order(-freq_null[, "frequency"]), ]

    # Compute null statistics
    entropy_freq_null[b] <- DescTools::Entropy(freq_null)
    entropy_eigen_null[b] <- DescTools::Entropy(igraph::eigen_centrality(net_null$g)$vector)
    
    #- print verbose if specified
    if(verbose){
      pb$tick()
    }
   
  }

  # smoothed bootstrap
  entropy_freq_null = as.vector(smoothed_bootstrap(entropy_freq_null))
  entropy_eigen_null = as.vector(smoothed_bootstrap(entropy_eigen_null))

  # Compute p-values
  p_value_eigen <- mean(entropy_eigen_null < entropy_eigen_obs)
  p_value_freq <- mean(entropy_freq_null < entropy_freq_obs)

  # # Determine significance
  # eigen_significance <- p_value_eigen < significance_level
  # freq_significance <- p_value_freq < significance_level

  # Return results
  ret = list(
    "p_value_eigen" = p_value_eigen,
    "p_value_freq" = p_value_freq,
    "significance_level" = significance_level
  )
  # return
  return(ret)
}
