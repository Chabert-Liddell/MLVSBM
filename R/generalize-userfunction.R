#' Create a GenMLVSBM object from observed data
#'
#' @param X A list of \eqn{L} square matrices with binary or count data,
#' ordered from upper to lower level
#' @param A A list of \eqn{L-1} matrices the affiliation matrix, linking X[[l]]
#' and X[[l+1]] with \code{nrow(A[[l]]) == nrow(X[[l+1]])} and
#' \code{ncol(A[[l]]) == nrow(X[[l]])}. Rows may have any number of affiliations
#' but will be rescaled so that each rows sum either to 1 or 0.
#' @param directed A vector of L booleans are the levels
#' directed or not. Default will check if the matrix are symmetric or not.
#' @param distribution A vector of length L, the distribution of X,
#' only "bernoulli" and "poisson" are implemented.
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
mlvsbm_create_generalized_network <-
  function(X, A, directed = NULL,
           distribution = NULL) {
    for (l in seq_along(A)) {
      if ( any(! rowSums(A[[l]]) %in% c(0,1))) {
#        warning(paste0("All rows of A must have exactly one 1!!!"))
        message(paste0("Affiliation has been normalized so that any rows sum to zero or one!"))
        A[[l]][rowSums(A[[l]]) != 0,] <- A[[l]][rowSums(A[[l]]) != 0,]/
          rowSums(A[[l]][rowSums(A[[l]]) != 0,])
      }
      if ( any(rowSums(A[[l]]) == 0)) {
        message(paste0("Some rows of A have no non 0 entry!!!"))
      }
    }
    for (l in seq_along(X)) {
      if (! is.matrix(X[[l]]) |
          ncol(X[[l]]) != nrow(X[[l]])) {
        stop(paste0("X[[l]] must be a square matrix!!!"))
      }
    }
    for (l in seq_along(A)) {
      if (ncol(X[[l+1]]) != nrow(A[[l]]) |
          ncol(X[[l]]) != ncol(A[[l]])) {
        stop(paste0("A[[l]], X[[l]] and X[[l+1]]'s dimensions are not compatible!!!"))
      }
    }

    if (is.null(directed)) {
      directed  <-  vapply(
        X, function(x) ! isSymmetric.matrix(x), FUN.VALUE = TRUE
      )
    }
    new_genmlvsbm <-
      GenMLVSBM$new(
        n = vapply(X, ncol, FUN.VALUE = 1),
        X = X,
        A = A,
        L = length(X),
        directed = directed,
        distribution = distribution
      )
    new_genmlvsbm$min_clusters <- rep(1, length(X))
    new_genmlvsbm$max_clusters <- vapply(X, function(x) floor(sqrt(nrow(x))),
                                      FUN.VALUE = 1)
    return (new_genmlvsbm)
  }

#' Create a simulated multilevel network (MLVSBM object)
#'
#' @param n A list of 2 positive integers,
#' the number of individuals and organizations.
#' @param Q A list of 2 positive integers,
#' the number of clusters of individuals and organizations.
#' @param pi A vector of probabilities of length Q_O,
#' the mixture parameter for the organizations.
#' @param alpha A list of 2 matrices, a \eqn{Q_I \times Q_I} matrix giving the
#' connectivity probabilities of the individuals and a \eqn{Q_O \times Q_O}
#' matrix giving the connectivity probabilities of the organizations.
#' @param directed A list of 2 logical. Is the individual level a directed
#' network ? Is the inter-organizational level a directed network?
#' @param gamma A \eqn{Q_I \times Q_O} matrix with each column summing to one,
#' the mixture parameters for the individuals
#' @param affiliation The distribution under which the affiliation matrix is
#' simulated in c("uniform", "preferential").
#' @param no_empty_org A logical with FALSE as default, should
#' every organizations have at least one affiliated individual?
#' Needs to have \eqn{n_I \geq n_O}.
#' @param distribution  A list for the distribution of X,
#' only "bernoulli" is implemented.
#' @param no_isolated_node A logical, if TRUE then the network is simulated
#' again until all nodes are connected.
#'
#' @return An MLVSBM object, a simulated multilevel network with levels,
#' affiliations and memberships.
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
                 L = L,
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

#' Infer a generalized multilevel network (GenMLVSBM object),
#' the original object is modified
#'
#' @description The inference use a greedy algorithm to navigate between model
#' size. For a given model size, the inference is done via a variational EM
#' algorithm. The returned model is the one with the highest ICL criterion among
#' all visited models.
#'
#' By default the algorithm fits a single level SBM for each level, before
#' inferring the generalized multilevel network. This step can be skipped by
#' specifying an initial clustering with the \code{init_clustering}. Also,
#' a given model size can be force by setting the parameters \code{nb_clusters}
#' to a given value.
#'
#' @param gmlv A GenMLVSBM object, the network to be inferred.
#' @param nb_clusters A vector of L integer, the model sizes.
#' If left to \code{NULL}, the algorithm
#' will navigate freely. Otherwise it will navigate between the specified model
#' size and its neighbors.
#' @param init_clustering A list of L vectors of integers of the same length as
#' the number of node of each level. If specified, the algorithm will start from
#' this clustering, then navigate freely.
#' @param nb_cores An integer, the number of cores to use. Default to \code{1}
#' for Windows and \code{detectCores()/2} for Linux and MacOS
#' @param init_method One of "hierarchical" (the default) or "spectral",
#' "spectral" might be more efficient but can lead to some numeric errors.
#' Not used when int_clustering is given.
#'
#' @return A FitGenMLVSBM object, the best inference of the network
#'
#' @importFrom blockmodels BM_bernoulli
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
#' \donttest{fit <- MLVSBM::mlvsbm_estimate_network(mlv = my_mlvsbm, nb_cores = 1)}
mlvsbm_estimate_generalized_network <-
  function(gmlv, nb_clusters = NULL, init_clustering = NULL, nb_cores = NULL,
           init_method = "hierarchical") {
    if (! "GenMLVSBM" %in% class(gmlv)) {
      stop("Object mlv must be of class GenMLVSBM,
            please use the function mlvsbm_create_network to create one")
    }
    os <- Sys.info()["sysname"]
    if (is.null(nb_cores)) {
      if (os != 'Windows') {
        nb_cores <-
          max(parallel::detectCores(all.tests = FALSE, logical = TRUE) %/% 2, 1)
        if (is.na(nb_cores)) nb_cores <- 1
      } else {
        nb_cores <-  1
      }
    }
   # browser()
    if (is.null(nb_clusters)) {
      if (is.null(init_clustering)) {
        fit_sbm <-
          bettermc::mclapply(
            X = seq_along(gmlv$adjacency_matrix),
            FUN = function(l) {
              sbm::estimateSimpleSBM(
                model = gmlv$distribution[l],
                netMat = gmlv$adjacency_matrix[[l]],
                estimOptions = list(verbosity = 0,
                                    plot = FALSE, nbCores = 1L,
                                    exploreMin = gmlv$max_clusters[l]))
            }, mc.cores = nb_cores
          )
        sbm_fit  <- lapply(X = seq(length(gmlv$nb_nodes)),
                           FUN = function(l) {
                             gmlv$estimate_level(level = l,
                                                 Z = fit_sbm[[l]]$memberships,
                                                 nb_cores = nb_cores,
                                                 init = "merged_split")
                           })
        # sbm_fit <- lapply(X = seq(length(gmlv$nb_nodes)),
        #                   FUN = function(l) {
        #                     gmlv$estimate_level(level = l,
        #                                         nb_cores = nb_cores,
        #                                        init = init_method)
        #                   })
        nb_clusters <- vapply(sbm_fit, function(bm) which.max(bm$ICL),
                              FUN.VALUE = 1)
        init_clustering <-
          lapply(seq_along(sbm_fit),
                 function(m) sbm_fit[[m]]$models[[nb_clusters[m]]]$Z)
        fit <- gmlv$estimate_all_bm(Q = nb_clusters,
                                    Z = init_clustering,
                                    nb_cores = nb_cores)
        ICL_sbm <- sum(vapply(sbm_fit,
                              function(bm) max(bm$ICL), FUN.VALUE = .1))
        print(paste0("ICL for independent levels : ", ICL_sbm))
        print(paste0("ICL for interdependent levels : ",
                     fit$ICL))
        if (ICL_sbm <= fit$ICL) {
          print("=====Interdependence is detected between levels!=====")
        } else {
          print("=====The levels of this network are independent!=====")
        }
        return(fit)
      } else {
        nb_clusters <- vapply(init_clustering, max, FUN.VALUE = 1)
        fit <- gmlv$estimate_all_bm(Q = nb_clusters,
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


#' Compute the complete log likelihood of a multilevel network for a given
#' clustering of the nodes.
#'
#' @description This function is useful to compute the likelihood for clusters
#' obtained by different methods.
#'
#' @param mlv A MLVSBM object, the network data
#' @param clustering A list of 2 vectors of integers of the same length as
#' the number of node of each level.
#' @return A numeric, the log likelihood of the multilevel network
#' for the given clustering.
#' @export
#'
#' @examples
#' my_mlvsbm <- MLVSBM::mlvsbm_simulate_network(
#'   n = list(I = 40, O = 20), # Number of nodes for the lower level and the upper level
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
#' mlvsbm_log_likelihood(mlv = my_mlvsbm, clustering = my_mlvsbm$memberships)
mlvsbm_log_likelihood <- function(mlv, clustering) {
  fit <- FitMLVSBM$new(Q = list(I = max(clustering[[1]]),
                                O = max(clustering[[2]])),
                       X = mlv$adjacency_matrix,
                       A = mlv$affiliation_matrix,
                       directed = mlv$directed,
                       distribution = mlv$distribution)
  fit$init_clustering(method = "merge_split",
                      Z = list(I = clustering[[1]], O = clustering[[2]]))
  fit$m_step()
  fit$complete_likelihood
}
