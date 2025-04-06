
#' plot.swaglm_network
#'
#' This function plots a \code{swaglm_network} object.
#'
#' @param x An object of class \code{swaglm_network}.
#' @param scale_vertex a \code{double} that specify the size of the nodes of the graph representing variables
#' @param ... Additional graphical parameters
#' @importFrom igraph layout.circle E degree
#' @importFrom fastglm fastglm
#' @examples
#' # Parameters for data generation
#' set.seed(12345)
#' n <- 2000
#' p <- 100
#' 
#' # create design matrix and vector of coefficients
#' Sigma <- diag(rep(1/p, p))
#' X <- MASS::mvrnorm(n = n, mu = rep(0, p), Sigma = Sigma)
#' beta = c(-15,-10,5,10,15, rep(0,p-5))
#' 
#' # --------------------- generate from logistic regression with an intercept of one
#' z <- 1 + X%*%beta
#' pr <- 1/(1 + exp(-z))
#' y <- as.factor(rbinom(n, 1, pr))
#' y = as.numeric(y)-1
#' 
#' # define swag parameters
#' quantile_alpha = .15
#' p_max = 20
#' swag_obj = swaglm::swaglm(X=X, y = y, p_max = p_max, family = stats::binomial(), 
#' alpha = quantile_alpha, verbose = TRUE, seed = 123)
#' names(swag_obj)
#' swag_network = compute_network(swag_obj)
#' plot(swag_network, scale_vertex = .05)
#' 
#' @export
plot.swaglm_network = function(x, scale_vertex = .1, ...){
  
  if(!inherits(x, "swaglm_network")){
    stop("Provided object 'x' needs to be of class 'swaglm_network'")
  }
  
  
  
  
  # x = swag_network
  plot(x$g_simplified_obs, layout = igraph::layout_in_circle(x$g_simplified_obs) , 
       vertex.color = "#e3ddff", edge.color = "black",  vertex.size = scale_vertex*igraph::degree(x$g), 
       edge.width = 0.1*igraph::E(x$g_simplified_obs)$width)
}




