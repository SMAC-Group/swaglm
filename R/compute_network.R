#' compute_network
#'
#' Compute a network representation of the selected models from a \code{swaglm_obj} object
#'
#' @param x An object of class \code{swaglm_obj}.
#' @param mode Character string specifying the mode of the network. Default is "undirected".
#' @importFrom igraph  graph_from_adjacency_matrix E
#' @importFrom plyr rbind.fill.matrix
#' @importFrom stats na.omit
#' @importFrom gdata upperTriangle
#' @return A list of class  \code{swaglm_network_obj} containing:
#' 
#' - \code{g}: The computed graph object.
#' - \code{models}: The selected models matrix.
#' - \code{g_simplified_obs}: The simplified network graph.
#' 
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
#' @export
compute_network <- function(x, mode = "undirected") {
  
  # x = res
  # mode = "undirected"
  # weighted = F
  # show_net = T
  
  if(!inherits(x, "swaglm_obj")){
    stop("Provided object 'x' needs to be of class 'swaglm_obj'")
  }
  


  # define number of maximum explored dimension
  p_max=length(x$lst_AIC)
  
  # selected_ind <- list()
  # 
  # selected_ind[[1]] <- obj$lst_var_mat[[1]][obj$lst_index_selected_models[[1]]+1, ]
  # ## add plus one because variables in R
  # selected_ind[[1]] = selected_ind[[1]]+1
  # 
  # for(i in 2:p_max) {
  # 
  # 
  #   if(length(obj$lst_index_selected_models[[i]]) != 1) {
  # 
  #     selected_ind[[i]] <- obj$lst_var_mat[[i]][as.matrix(obj$lst_index_selected_models[[i]])+1, ]
  #     selected_ind[[i]] = selected_ind[[i]] +1
  #   } else {
  #     selected_ind[[i]] = obj$lst_var_mat[[i]][obj$lst_index_selected_models[[i]]+1, ]
  #     selected_ind[[i]] = t(selected_ind[[i]] +1)
  #   }
  # 
  # }
  # 
  # 
  # models <- rbind.fill.matrix(selected_ind)

  models = plyr::rbind.fill.matrix(x$lst_selected_models) +1
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

  colnames(intensity) <- x$lst_selected_models[[1]]+1
  rownames(intensity) <- x$lst_selected_models[[1]]+1

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

  relation_mat[,3] = gdata::upperTriangle(intensity,byrow = TRUE)

  g <- igraph::graph_from_adjacency_matrix(intensity, mode = mode, weighted = NULL)

  vertex_degrees_obs <- igraph::degree(g)

  igraph::E(g)$weight <- 1

  # Simplify the graph, summing the weights of multiple edges
  g_simplified_obs <- igraph::simplify(g, edge.attr.comb = list(weight = "sum"))

  # Set the edge width based on the combined weights
  igraph::E(g_simplified_obs)$width <- igraph::E(g_simplified_obs)$weight
  
  # define return object
  ret = list("g" = g,
             "models" = models,
             "g_simplified_obs" = g_simplified_obs)

  # define class of return object
  class(ret) = "swaglm_network_obj"
  # if(show_net == T) plot(g_simplified_obs, layout = layout.circle, vertex.color = "skyblue", edge.color = "black", vertex.size = 0.1*degree(g), edge.width = 0.1*E(g_simplified_obs)$width )

  # return return object
  return(ret)

}





#' Plot a \code{swaglm_network_obj} object
#'
#' This function plots a \code{swaglm_network_obj} object.
#'
#' @param x An object of class \code{swaglm_network_obj}.
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
plot.swaglm_network_obj = function(x, scale_vertex = .1, ...){
  
  if(!inherits(x, "swaglm_network_obj")){
    stop("Provided object 'x' needs to be of class 'swaglm_network_obj'")
  }
  
  
  # x = swag_network
  plot(x$g_simplified_obs, layout = igraph::layout_in_circle(x$g_simplified_obs) , 
       vertex.color = "#e3ddff", edge.color = "black",  vertex.size = scale_vertex*igraph::degree(x$g), 
       edge.width = 0.1*igraph::E(x$g_simplified_obs)$width)
}




