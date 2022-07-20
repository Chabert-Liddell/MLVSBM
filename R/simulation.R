MLVSBM$set(
  "public", "simulate",
# Simulate a multilevel network from a MLVSBM object
  function() {
    if (is.null(private$A))
      private$A <- simulate_affiliation(private$n$I,
                                        private$n$O,
                                        affiliation = private$sim_param$affiliation,
                                        no_empty_org = private$sim_param$no_empty_org)
    private$Z    <-  list(I = NULL, O = NULL)
    private$Z$O  <-  sample(
      x = seq(private$sim_param$Q$O),
      size = private$n$O,
      replace = TRUE,
      prob = private$sim_param$pi$O) # Variables latentes L
    private$Z$I  <-  numeric(private$n$I)
    ind          <- private$A %*% private$Z$O
    for (i in seq(private$n$I)) {
      private$Z$I[i] <- sample(seq(private$sim_param$Q$I),
                               size = 1,
                               prob = private$sim_param$gamma[, ind[i, ]])
    }
    private$X <- list(I = NULL, O = NULL)
    private$X$I <- simulate_adjacency(Z = private$Z$I,
                                     n = private$n$I,
                                     alpha = private$sim_param$alpha$I,
                                     directed = private$directed_$I,
                                     distribution = private$distribution_$I,
                                     no_isolated_node = private$sim_param$no_isolated_node)
    private$X$O <- simulate_adjacency(Z = private$Z$O,
                                     n = private$n$O,
                                     alpha = private$sim_param$alpha$O,
                                     directed = private$directed_$O,
                                     distribution = private$distribution_$O,
                                     no_isolated_node = private$sim_param$no_isolated_node)
    names(private$Z$I) <- paste0("I", seq(length(private$Z$I)))
    names(private$Z$O) <- paste0("O", seq(length(private$Z$O)))
    rownames(private$X$I) <- colnames(private$X$I) <- names(private$Z$I)
    rownames(private$X$O) <- colnames(private$X$O) <- names(private$Z$O)
    rownames(private$A) <- names(private$Z$I)
    colnames(private$A) <- names(private$Z$O)
    }
)


GenMLVSBM$set(
  "public", "simulate",
  # Simulate a multilevel network from a MLVSBM object
  function() {
    if (is.null(private$A)) {
      if(private$sim_param$affiliation == "diagonal") {
        private$A <- lapply(seq(private$L), function(l) diag(1, private$n[l]))
      } else {
        private$A <-
          lapply(seq(private$L -1),
            function(l) {
              simulate_affiliation(n = private$n[l+1],
                                   m = private$n[l],
                                   affiliation = private$sim_param$affiliation,
                                   no_empty_org = private$sim_param$no_empty_org)
            }
          )
      }
    }
    private$Z    <-  vector("list", private$L)
    private$Z[[1]] <-  sample(
      x = seq(private$sim_param$Q[1]),
      size = private$n[1],
      replace = TRUE,
      prob = private$sim_param$pi[[1]]) # Variables latentes L
    for (l in seq(2, private$L)) {
      private$Z[[l]] <-  numeric(private$n[l])
      ind          <- private$A[[l-1]] %*% private$Z[[l-1]]
      for (i in seq(private$n[l])) {
        if (ind[i,] > 0) {
          private$Z[[l]][i] <- sample(seq(private$sim_param$Q[l]),
                                      size = 1,
                                      prob = private$sim_param$gamma[[l-1]][, ind[i, ]])
        } else {
          private$Z[[l]][i] <- sample(seq(private$sim_param$Q[l]),
                                      size = 1,
                                      prob = private$sim_param$pi[[l]])
        }
      }
    }
    private$X <- vector("list", private$L)
    private$X <- lapply(
      X = seq(private$L),
      FUN = function(l) {
        simulate_adjacency(Z = private$Z[[l]],
                           n = private$n[l],
                           alpha = private$sim_param$alpha[[l]],
                           directed = private$directed_[l],
                           distribution = private$distribution_[l],
                           no_isolated_node = private$sim_param$no_isolated_node)
      }
    )
  }
)



#' Simulation an adjacency matrix
#'
#' @importFrom stats rbinom
#'
#' @param Z A vector of integer of size n, the label
#' @param n An integer, the number of rows or columns of the matrix
#' @param alpha A \eqn{\max(Z)\times \max(Z)} matrix, the connectivity parameters
#' @param directed A boolean, Is the network directed or not ?
#' @param distribution The distribution of the indices: only "bernoulli"
#' @param no_isolated_node A boolean, may row and column of adjacency matrices sum to 0
#' @return A nxn adjacency matrix
simulate_adjacency <- function(Z, n, alpha, directed,
                               distribution = "bernoulli", no_isolated_node = FALSE) {
  X <- matrix(0, n, n)
  condition <- TRUE
  it <- 1
  while (condition) {
    X[] <- 0
    for (i in 1:(n-1)){
      X[i, (i+1):n] <- stats::rbinom(n-i, 1, alpha[Z[i], Z[(i+1):n]])
    }
    if (! directed) {
      X <- X + t(X)
    } else {
      for (i in 2:n) {
        X[i, 1:(i-1)] <-
          stats::rbinom(i-1, 1, alpha[Z[i], Z[1:(i-1)]])
      }
    }
    condition <- no_isolated_node & it < 100 & any(colSums(X) + rowSums(X) == 0 )
    it <- it + 1
    if (it == 100 & no_isolated_node) warning("Could not generate a fully connected network!")
  }
  return(X)
}



#' Simulate of matrix of affiliation
#'
#' @param n An integer, the number of individuals
#' @param m An integer, the number of organizations
#' @param affiliation The type of affiliation between c("uniform", "preferential")
#' @param no_empty_org A Boolean. Force all columns to have at least a 1.
#' Must have \eqn{n>m}.
#'
#' @return A \eqn{n \times m} affiliation matrix, with a unique 1 on each rows
simulate_affiliation <-
  function(n, m, affiliation = "uniform", no_empty_org = FALSE) {
    if(is.null(m)) {
      m = n[2]
      n = n[1]
    }
    A <- matrix(0, n, m)
    if (no_empty_org) {
      A[1:m, 1:m] <- diag(1, m, m)
      if (affiliation == "uniform") {
        for (i in seq(m+1, n)) {
          j       <-  sample(x = seq(m), size = 1)
          A[i, j] <- 1
        }
      }
      if (affiliation == "preferential") {
        for (i in seq(m+1, n)) {
          j <- sample( x = seq(m), size = 1, prob = c(colSums(A)))
          A[i, j] <- 1
        }
      }
    } else {
      if (affiliation == "uniform") {
        ind <- sample(x = seq(m), size = n, replace = TRUE)
        A[matrix(c(seq(n), ind), ncol = 2)] <- 1
      }
      if (affiliation == "preferential") {
        for (i in seq(n)) {
          j <- sample( x = seq(m), size = 1, prob = c(colSums(A)) + 1)
          A[i, j] <- 1
        }
      }
    }
    return(A)
  }
