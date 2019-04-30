

#' Create a MLVSBM object from observed data
#'
#' @param X_I A matrix square binary adjacency, the individual level
#' @param X_O A matrix square binary adjacency, the organisational level
#' @param A A matrix the affiliation matrix with individuals in rows and
#' organisations in columns
#'
#' @return An unfitted MLVSBM object corresponding to the multilevel network
#' @export
#'
#' @examples
mlvsbm_create_network <-
  function(X_I, X_O, A) {
    if (! is.matrix(A)) {
      cat(paste0("A must be a binary matrix!!!"))
    }
    if ( any(rowSums(A) != 1) |
         any(A != 0 & A != 1)) {
      cat(paste0("All rows of A must have exactly one 1!!!"))
    }
    if (! is.matrix(X_I) |
        any(X_I != 0 & X_I != 1) |
        ncol(X_I) != nrow(X_I)) {
      cat(paste0("X_I must be a square binary matrix!!!"))
    }
    if (! is.matrix(X_O) |
        any(X_O != 0 & X_I != 1)|
        ncol(X_O) != nrow(X_O)) {
      cat(paste0("X_O must be a square binary matrix!!!"))
    }
    if (ncol(X_I) != nrow(A) |
        ncol(X_O) != ncol(A)) {
      cat(paste0("A, X_I and X_O's dimensions are not compatible!!!"))
    }
    new_mlvsbm <-
      MLVSBM$new(
        n = list("I" = nrow(A),
                 "O" = ncol(A)),
        X = list("I" = X_I,
                 "O" = X_O),
        A = A,
        directed = list("I" = ! isSymmetric(X_I),
                        "O" = ! isSymmetric(X_O))
                            )
    new_mlvsbm$min_clusters <- list("I" = 1, "O" = 1)
    new_mlvsbm$max_clusters <- list("I" = floor(sqrt(nrow(A))),
                                    "O" = floor(sqrt(ncol(A))))
    return (new_mlvsbm)
  }

#' Create a simulated multilevel network (MLVSBM object)
#'
#' @param n_I A positive integer, the number of individuals
#' @param n_O A positive integer, the number of organisations
#' @param Q_I A positive integer, the number of clusters of individuals
#' @param Q_O A positive intger, the number of clusters of organisations
#' @param pi_O A probability vetor of length Q_O,
#' the mixture parameter for the organisations
#' @param gamma A Q_IxQ_O matrix with each column suming to one,
#' the mixture parameters for the individuals
#' @param alpha_I A Q_IxQ_I matrix giving the connectivity probabilities
#' of the individuals
#' @param alpha_O A Q_OxQ_O matrix giving the connectivity probabilities
#' of the organisations
#' @param directed_I A boolean, is level I a directed network?
#' @param directed_O A boolean, is level O a directed network?
#' @param affiliation The distribution under which the affiliation matrix is
#' simulated in c("uniform", "preferential")
#' @param no_empty_org A boolean with FALSE as default, should
#' every organisation have at least one affiliated individuals?
#' Needs to have n_I >= n_O
#' @return An MLVSBM object, a simulated multilevel network with levels,
#' affiliations and memberships
#' @export
#'
#' @examples
mlvsbm_simulate_network <-
  function (n_I, n_O, Q_I, Q_O, pi_O, gamma, alpha_I, alpha_O,
            directed_I, directed_O, affiliation = "uniform",
            no_empty_org = FALSE) {
    if (n_I < 1 | n_O < 1 | n_I%%1 != 0 | n_O %% 1 != 0) {
      cat(paste0("n_I and n_O must be positive integers!!!"))
    }
    if (! is.matrix(gamma) | ! is.matrix(alpha_I) | ! is.matrix(alpha_O)) {
      cat(paste0("gamma, alpha_I and alpha_O must be matrices!!!"))
    }
    if (any(pi_O > 1) | any(pi_O < 0) | sum(pi_O) != 1) {
      cat(paste0("pi_O is a probability vector,
                 its coefficients must sum to one!!!"))
    }
    if ( Q_O != length(pi_O) |
         Q_O != ncol(gamma) |
         Q_O != nrow(alpha_O) |
         Q_O != ncol(alpha_O) |
         Q_I != nrow(alpha_I) |
         Q_I != ncol(alpha_I) |
         Q_I != nrow(gamma)) {
      cat(paste0("Number of clusters and parameters dimension are
                 not compatible!!!"))
    }
    if (any(gamma > 1) | any(gamma < 0) | any(colSums(gamma) != 1)) {
      cat(paste0("Any column of gamma must be a probability vector!!!"))
    }
    if (any(alpha_I > 1) | any(alpha_I < 0)) {
      cat(paste0("Any coefficient of alpha_I must be between 0 and 1!!!"))
    }
    if (any(alpha_O > 1) | any(alpha_O < 0)) {
      cat(paste0("Any coefficient of alpha_I must be between 0 and 1!!!"))
    }
    if (directed_I & !isSymmetric(alpha_I)) {
      cat(paste0("alpha_I is not symmetric but level is directed!!!"))
    }
    if (directed_O & !isSymmetric(alpha_O)) {
      cat(paste0("alpha_O is not symmetric but level is directed!!!"))
    }
    new_mlvsbm <- MLVSBM$new(n = list("I" = n_I, "O" = n_O),
                             directed = list("I" = directed_I, "O" = directed_O),
                             sim_param = list(alpha = list("I" = alpha_I,
                                                           "O" = alpha_O),
                                              gamma = gamma,
                                              pi = list("O" = pi_O),
                                              Q = list("I" = Q_I, "O" = Q_O)))
    new_mlvsbm$min_clusters <- list("I" = 1, "O" = 1)
    new_mlvsbm$max_clusters <- list("I" = floor(sqrt(n_I)),
                                    "O" = floor(sqrt(n_O)))
    new_mlvsbm$simulate(affiliation = affiliation,
                        no_empty_org = no_empty_org)
    return(new_mlvsbm)
  }
