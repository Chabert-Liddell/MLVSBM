

#' Create a MLVSBM object from observed data
#'
#' @param X A list of 2 squares binary matrices,
#' the first one being the individual or lower level
#' the second one being the organizational or upper level
#' @param A A matrix the affiliation matrix with individuals in rows and
#' organizations in columns
#' @param directed A list of 2 boolean are the upper and lower level
#' directed or not. Default will check if the matrix are symmetric or not.
#' @param distribution A list for the distribution of X,
#' only "bernoulli" is implemented
#'
#' @return An unfitted MLVSBM object corresponding to the multilevel network
#' @importFrom stats rbinom
#' @export
#' @examples
#' ind_adj <- matrix(stats::rbinom(n = 10**2, size = 1, prob = .2),
#'                   nrow = 10, ncol = 10)
#' org_adj <- matrix(stats::rbinom(n = 10**2, size = 1, prob = .3),
#'                   nrow = 10, ncol = 10)
#' affiliation <- diag(1, 10)
#' my_mlvsbm <- mlvsbm_create_network(X = list(I = ind_adj, O = org_adj),
#'                                    directed = list(I = FALSE, O = FALSE),
#'                                    A = affiliation)
mlvsbm_create_network <-
  function(X, A, directed = NULL,
           distribution = list("bernoulli", "bernoulli")) {
    if (! is.matrix(A)) {
      warning(paste0("A must be a binary matrix!!!"))
    }
    if ( any(rowSums(A) != 1) |
         any(A != 0 & A != 1)) {
      warning(paste0("All rows of A must have exactly one 1!!!"))
    }
    if (! is.matrix(X[[1]]) |
        any(X[[1]] != 0 & X[[1]] != 1, na.rm = TRUE) |
        ncol(X[[1]]) != nrow(X[[1]])) {
      warning(paste0("X[[1]] must be a square binary matrix!!!"))
    }
    if (! is.matrix(X[[2]]) |
        any(X[[2]] != 0 & X[[2]] != 1, na.rm = TRUE)|
        ncol(X[[2]]) != nrow(X[[2]])) {
      warning(paste0("X[[2]] must be a square binary matrix!!!"))
    }
    if (ncol(X[[1]]) != nrow(A) |
        ncol(X[[2]]) != ncol(A)) {
      stop(paste0("A, X[[1]] and X[[2]]'s dimensions are not compatible!!!"))
    }
    if (is.null(directed)) {
      directed  <-  list(I = ! isSymmetric(X[[1]]),
                         O = ! isSymmetric(X[[2]]))
    }
    new_mlvsbm <-
      MLVSBM$new(
        n = list(I = nrow(A),
                 O = ncol(A)),
        X = list(I = X[[1]],
                 O = X[[2]]),
        A = A,
        directed = directed,
        distribution = list(I = distribution[[1]],
                            O = distribution[[2]])
        )
    new_mlvsbm$min_clusters <- list(I = 1, O = 1)
    new_mlvsbm$max_clusters <- list(I = floor(sqrt(nrow(A))),
                                    O = floor(sqrt(ncol(A))))
    return (new_mlvsbm)
  }

#' Create a simulated multilevel network (MLVSBM object)
#'
#' @param n A list of 2 positive integers,
#' the number of individuals and organisations
#' @param Q A list of 2 positive integers,
#' the number of clusters of individuals and organisations
#' @param pi A vector of probabilities of length Q_O,
#' the mixture parameter for the organisations
#' @param alpha A list of 2 matrices, a Q_IxQ_I matrix giving the connectivity probabilities
#' of the individuals and a Q_OxQ_O matrix giving the connectivity probabilities
#' of the organisations
#' @param directed A list of 2 booleans. Is the individual level a directed network ?
#' Is the organisational level a directed network ?
#' @param gamma A Q_IxQ_O matrix with each column suming to one,
#' the mixture parameters for the individuals
#' @param affiliation The distribution under which the affiliation matrix is
#' simulated in c("uniform", "preferential")
#' @param no_empty_org A boolean with FALSE as default, should
#' every organisation have at least one affiliated individuals?
#' Needs to have n_I >= n_O
#' @param distribution  A list for the distribution of X,
#' only "bernoulli" is implemented
#' @param no_isolated_node A boolean, if TRUE then the network is resimulated
#' until all nodes are connected
#'
#' @return An MLVSBM object, a simulated multilevel network with levels,
#' affiliations and memberships
#' @export
#'
#' @examples
#' my_mlvsbm <- MLVSBM::mlvsbm_simulate_network(
#'   n = list(I = 10, O = 20), # Number of nodes for the lower level and the upper level
#'   Q = list(I = 2, O = 2), # Number of blocks for the lower level and the upper level
#'   pi = c(.3, .7), # Block proportion for the upper level, must sum to one
#'   gamma = matrix(c(.9, .2,   # Block proportion for the lower level,
#'                    .1, .8), # each column must sum to one
#'                  nrow = 2, ncol = 2, byrow = TRUE),
#'   alpha = list(I = matrix(c(.8, .2,
#'                             .2, .1),
#'                           nrow = 2, ncol = 2, byrow = TRUE), # Connection matrix
#'                O = matrix(c(.99, .3,
#'                             .3, .1),
#'                           nrow = 2, ncol = 2, byrow = TRUE)),# between blocks
#'   directed = list(I = FALSE, O = FALSE)) # Are the upper and lower level directed
mlvsbm_simulate_network <-
  function (n, Q, pi, gamma, alpha,
            directed, affiliation = "uniform",
            distribution = list("bernoulli", "bernoulli"),
            no_empty_org = FALSE,
            no_isolated_node = FALSE) {
    if (n[[1]] < 1 | n[[2]] < 1 | n[[1]]%%1 != 0 | n[[2]] %% 1 != 0) {
      stop(paste0("n[[1]] and n[[2]] must be positive integers!!!"))
    }
    if (! is.matrix(gamma) | ! is.matrix(alpha[[1]]) | ! is.matrix(alpha[[2]])) {
      stop(paste0("gamma, alpha[[1]] and alpha[[2]] must be matrices!!!"))
    }
    if (any(pi > 1) | any(pi < 0) | sum(pi) != 1) {
      warning(paste0("pi is a probability vector,
                 its coefficients must sum to one!!!"))
    }
    if ( Q[[2]] != length(pi) |
         Q[[2]] != ncol(gamma) |
         Q[[2]] != nrow(alpha[[2]]) |
         Q[[2]] != ncol(alpha[[2]]) |
         Q[[1]] != nrow(alpha[[1]]) |
         Q[[1]] != ncol(alpha[[1]]) |
         Q[[1]] != nrow(gamma)) {
      stop(paste0("Number of clusters and parameters dimension are
                 not compatible!!!"))
    }
    if (any(gamma > 1) | any(gamma < 0) | any(colSums(gamma) != 1)) {
      stop(paste0("Any column of gamma must be a probability vector!!!"))
    }
    if (any(alpha[[1]] > 1) | any(alpha[[1]] < 0)) {
      stop(paste0("Any coefficient of alpha[[1]] must be between 0 and 1!!!"))
    }
    if (any(alpha[[2]] > 1) | any(alpha[[2]] < 0)) {
      stop(paste0("Any coefficient of alpha[[1]] must be between 0 and 1!!!"))
    }
    if (! directed[[1]] & !isSymmetric(alpha[[1]])) {
      stop(paste0("alpha[[1]] is not symmetric but level is undirected!!!"))
    }
    if (! directed[[2]] & ! isSymmetric(alpha[[2]])) {
      stop(paste0("alpha[[2]] is not symmetric but level is undirected!!!"))
    }
    new_mlvsbm <-
      MLVSBM$new(n = list(I = n[[1]],
                          O = n[[2]]),
                 directed = list(I = directed[[1]],
                                 O = directed[[2]]),
                 sim_param = list(alpha = list(I = alpha[[1]],
                                               O = alpha[[2]]),
                                  gamma = gamma,
                                  pi = list(O = pi),
                                  Q = list(I = Q[[1]],
                                           O = Q[[2]]),
                                  affiliation = affiliation,
                                  no_empty_org = no_empty_org,
                                  no_isolated_node = no_isolated_node),
                 distribution = list(I = distribution[[1]],
                                     O = distribution[[2]]))
    new_mlvsbm$simulate()
    return(new_mlvsbm)
  }

#' Infer a multilevel network (MLVSBM object), the original object is modified
#'
#' @description The inference use a greedy algorithm to navigate between model
#' size. For a given model size, the inference is done via a variational EM
#' algorithm. The returned model is the one with the highest ICL criterion among
#' all visited models.
#'
#' By default the algorithm fits a single level SBM for each level, before
#' inferring the multilevel network. This step can be skipped by specifying an
#' initial clustering with the \code{init_clustering}. Also, a given model size
#' can be force by setting the parameters \code{nb_clusters} to a given value.
#'
#' @param mlv A MLVSBM object, the network to be inferred
#' @param nb_clusters A list of 2 integers, the model size.
#' If left to \code{NULL}, the algorithm
#' will navigate freely. Otherwise it will navigate between the specified model
#' size and its neighbours.
#' @param init_clustering A list of 2 vectors of integers of the same length as
#' the number of node of each level. If specified, the algorithm will start from
#' this clustering, then navigate freely.
#' @param nb_cores An integer, the number of cores to use. Default to \code{1}
#' for Windows and \code{detectCores()/2} for Linux and MacOS
#'
#' @return A FitMLVSBM object, the best inference of the network
#' @export

#' @examples
#' my_mlvsbm <- MLVSBM::mlvsbm_simulate_network(
#'   n = list(I = 10, O = 20), # Number of nodes for the lower level and the upper level
#'   Q = list(I = 2, O = 2), # Number of blocks for the lower level and the upper level
#'   pi = c(.3, .7), # Block proportion for the upper level, must sum to one
#'   gamma = matrix(c(.9, .2,   # Block proportion for the lower level,
#'                    .1, .8), # each column must sum to one
#'                  nrow = 2, ncol = 2, byrow = TRUE),
#'   alpha = list(I = matrix(c(.8, .2,
#'                             .2, .1),
#'                           nrow = 2, ncol = 2, byrow = TRUE), # Connection matrix
#'                O = matrix(c(.99, .3,
#'                             .3, .1),
#'                           nrow = 2, ncol = 2, byrow = TRUE)),# between blocks
#'   directed = list(I = FALSE, O = FALSE), # Are the upper and lower level directed or not ?
#'   affiliation = "preferential") # How the affiliation matrix is generated
#' fit <- MLVSBM::mlvsbm_estimate_network(mlv = my_mlvsbm, nb_cores = 1)
mlvsbm_estimate_network <-
  function(mlv, nb_clusters = NULL, init_clustering = NULL, nb_cores = NULL) {
    if (! "MLVSBM" %in% class(mlv)) {
      stop("Object mlv must be of class MLVSBM,
            please use the function mlvsbm_create_network to create one")
    }
    if (is.null(nb_clusters)) {
      if (is.null(init_clustering)) {
        lower_fit <- mlv$estimate_level(level = "lower", nb_cores = nb_cores)
        upper_fit <- mlv$estimate_level(level = "upper", nb_cores = nb_cores)
        nb_clusters <- list("I" = which.max(lower_fit$ICL),
                            "O" = which.max(upper_fit$ICL))
        init_clustering <-
          list("I" = lower_fit$models[[which.max(lower_fit$ICL)]]$Z,
               "O" = upper_fit$models[[which.max(upper_fit$ICL)]]$Z)
        fit <- mlv$estimate_all_bm(Q = nb_clusters,
                                   Z = init_clustering,
                                   nb_cores = nb_cores)
        print(paste0("ICL for independent levels : ",
                     max(lower_fit$ICL) + max(upper_fit$ICL)))
        print(paste0("ICL for interdependent levels : ",
                     fit$ICL))
        if (max(lower_fit$ICL) + max(upper_fit$ICL) <= fit$ICL) {
          print("=====Interdependence is detected between the two levels!=====")
        } else {
          print("=====The levels of this network are independent!=====")
        }
        return(fit)
      } else {
        nb_clusters <- list("I" = max(init_clustering[[1]]),
                            "O" = max(init_clustering[[2]]))
        fit <- mlv$estimate_all_bm(Q = nb_clusters,
                                   Z = init_clustering,
                                   nb_cores = nb_cores)
        print(paste0("ICL for interdependent levels : ",
                     fit$ICL))
        return(fit)
      }
    } else {
      fit <- mlv$estimate_one(Q = nb_clusters,
                       Z = init_clustering,
                       nb_cores = nb_cores)
      print(paste0("ICL for interdependent levels : ",
                   fit$ICL))
      return(fit)
    }
  }
