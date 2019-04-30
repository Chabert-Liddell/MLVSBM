

#' Create a MLVSBM object from observed data
#'
#' @param A A matrix the affiliation matrix with individuals in rows and
#' organisations in columns
#' @param X A list of 2 squares binary matrices,
#' the first one being the individual or lower level
#' the second one being the organisational or upper level
#' @param distribution A list for the distribution of X,
#' only "bernoulli" is implemented
#'
#' @return An unfitted MLVSBM object corresponding to the multilevel network
#' @export
#'
#' @examples
mlvsbm_create_network <-
  function(X, A, distribution = list("bernoulli", "bernoulli")) {
    if (! is.matrix(A)) {
      cat(paste0("A must be a binary matrix!!!"))
    }
    if ( any(rowSums(A) != 1) |
         any(A != 0 & A != 1)) {
      cat(paste0("All rows of A must have exactly one 1!!!"))
    }
    if (! is.matrix(X[[1]]) |
        any(X[[1]] != 0 & X[[1]] != 1) |
        ncol(X[[1]]) != nrow(X[[1]])) {
      cat(paste0("X[[1]] must be a square binary matrix!!!"))
    }
    if (! is.matrix(X[[2]]) |
        any(X[[2]] != 0 & X[[1]] != 1)|
        ncol(X[[2]]) != nrow(X[[2]])) {
      cat(paste0("X[[2]] must be a square binary matrix!!!"))
    }
    if (ncol(X[[1]]) != nrow(A) |
        ncol(X[[2]]) != ncol(A)) {
      cat(paste0("A, X[[1]] and X[[2]]'s dimensions are not compatible!!!"))
    }
    new_mlvsbm <-
      MLVSBM$new(
        n = list(I = nrow(A),
                 O = ncol(A)),
        X = list(I = X[[1]],
                 O = X[[2]]),
        A = A,
        directed = list(I = ! isSymmetric(X[[1]]),
                        O = ! isSymmetric(X[[2]])),
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
#' @param distribution  list for the distribution of X,
#' only "bernoulli" is implemented
#'
#' @return An MLVSBM object, a simulated multilevel network with levels,
#' affiliations and memberships
#' @export
#'
#' @examples
mlvsbm_simulate_network <-
  function (n, Q, pi, gamma, alpha,
            directed, affiliation = "uniform",
            distribution = list("bernoulli", "bernoulli"),
            no_empty_org = FALSE) {
    if (n[[1]] < 1 | n[[2]] < 1 | n[[1]]%%1 != 0 | n[[2]] %% 1 != 0) {
      cat(paste0("n[[1]] and n[[2]] must be positive integers!!!"))
    }
    if (! is.matrix(gamma) | ! is.matrix(alpha[[1]]) | ! is.matrix(alpha[[2]])) {
      cat(paste0("gamma, alpha[[1]] and alpha[[2]] must be matrices!!!"))
    }
    if (any(pi > 1) | any(pi < 0) | sum(pi) != 1) {
      cat(paste0("pi is a probability vector,
                 its coefficients must sum to one!!!"))
    }
    if ( Q[[2]] != length(pi) |
         Q[[2]] != ncol(gamma) |
         Q[[2]] != nrow(alpha[[2]]) |
         Q[[2]] != ncol(alpha[[2]]) |
         Q[[1]] != nrow(alpha[[1]]) |
         Q[[1]] != ncol(alpha[[1]]) |
         Q[[1]] != nrow(gamma)) {
      cat(paste0("Number of clusters and parameters dimension are
                 not compatible!!!"))
    }
    if (any(gamma > 1) | any(gamma < 0) | any(colSums(gamma) != 1)) {
      cat(paste0("Any column of gamma must be a probability vector!!!"))
    }
    if (any(alpha[[1]] > 1) | any(alpha[[1]] < 0)) {
      cat(paste0("Any coefficient of alpha[[1]] must be between 0 and 1!!!"))
    }
    if (any(alpha[[2]] > 1) | any(alpha[[2]] < 0)) {
      cat(paste0("Any coefficient of alpha[[1]] must be between 0 and 1!!!"))
    }
    if (directed[[1]] & !isSymmetric(alpha[[1]])) {
      cat(paste0("alpha[[1]] is not symmetric but level is directed!!!"))
    }
    if (directed[[2]] & !isSymmetric(alpha[[2]])) {
      cat(paste0("alpha[[2]] is not symmetric but level is directed!!!"))
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
                                  no_empty_org = no_empty_org),
                 distribution = list(I = distribution[[1]],
                                     O = distribution[[2]]))
    new_mlvsbm$simulate()
    return(new_mlvsbm)
  }
