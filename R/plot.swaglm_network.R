#' Plot method for swaglm_network objects
#'
#' Visualizes a SWAG network with discretized vertex size, optional edge width scaling,
#' and edge coloring based on a correlation matrix.
#'
#'
#' @param x An object of class \code{swaglm_network}
#' @param scale_vertex Optional vertex scaling parameter
#' @param bins Number of bins for vertex size discretization (default = 5)
#' @param ... Additional arguments passed to \code{igraph::plot}
#' @importFrom igraph layout.circle E degree
#' @importFrom fastglm fastglm
#' @examples
#' # Parameters for data generation
#' set.seed(12345)
#' n <- 2000
#' p <- 50
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
#' plot(swag_network)
#' 
#' @export
#' @export
#' 
#
#' 
#' 
library(igraph)
library(scales)  # for rescale
library(RColorBrewer)
library(fields)  # for image.plot

plot.swaglm_network = function(x, bins = 5, scale_edge = NULL, size_range = c(8, 30), ...) {
  
  if (!inherits(x, "swaglm_network")) {
    stop("Provided object 'x' needs to be of class 'swaglm_network'")
  }
  
  degs <- igraph::degree(x$g)
  
  # --- Discretize vertex sizes based on degree
  degree_bins <- cut(degs,
                     breaks = quantile(degs, probs = seq(0, 1, length.out = bins + 1), na.rm = TRUE),
                     include.lowest = TRUE, labels = FALSE)
  
  size_seq <- seq(size_range[1], size_range[2], length.out = bins)
  vertex_sizes <- size_seq[degree_bins]
  
  # --- Edge width scaling: automatic if not specified
  edge_weights <- igraph::E(x$g_simplified_obs)$width
  
  if (is.null(scale_edge)) {
    edge_min <- min(edge_weights)
    edge_max <- max(edge_weights)
    if (edge_max == edge_min) {
      edge_scaled <- rep(1, length(edge_weights))  # All equal
    } else {
      edge_scaled <- 0.5 + 2.5 * (edge_weights - edge_min) / (edge_max - edge_min)
    }
  } else {
    edge_scaled <- scale_edge * edge_weights
  }
  
  plot(
    x$g_simplified_obs,
    layout = igraph::layout_in_circle(x$g_simplified_obs),
    vertex.color = "#e3ddff",
    edge.color = "black",
    vertex.size = vertex_sizes,
    edge.width = edge_scaled,
    ...
  )
}




