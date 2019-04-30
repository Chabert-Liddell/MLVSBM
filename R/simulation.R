MLVSBM$set(
  "public", "simulate",
  function(affiliation = "uniform", no_empty_org = FALSE) {
    if (is.null(private$A))
      private$A <- simulate_affiliation(private$n$I,
                                        private$n$O,
                                        affiliation = affiliation,
                                        no_empty_org = no_empty_org)
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
    private$X$I = simulate_adjacency(Z = private$Z$I,
                                     n = private$n$I,
                                     alpha = private$sim_param$alpha$I,
                                     directed = private$directed$I,
                                     distribution = private$distribution$I)
    private$X$O = simulate_adjacency(Z = private$Z$O,
                                     n = private$n$O,
                                     alpha = private$sim_param$alpha$O,
                                     directed = private$directed$O,
                                     distribution = private$distribution$O)
    })

#' Simulation an adjacency matrix
#'
#' @importFrom stats rbinom
#'
#' @param Z A vector of integer of size n, the label
#' @param n An intger, the number of rows or columbns of the matrix
#' @param alpha A max(Z)xmax(Z) matrix, the connectivity parameters
#' @param directed A boolean, Is the network directed or not ?
#' @param distribution The distribution of the indices: only "bernouilli"
#'
#' @return A nxn adjacency matrix
#' @export
#'
#' @examples
simulate_adjacency <- function(Z, n, alpha, directed, distribution) {
  X = matrix(0, n, n)
  for (i in 1:(n-1)){
    X[i, (i+1):n] = stats::rbinom(n-i, 1, alpha[Z[i], Z[(i+1):n]])
  }
  if (! directed) {
    X = X + t(X)
  } else {
    for (i in 2:n) {
      X[i, 1:(i-1)] =
        stats::rbinom(i-1, 1, alpha[Z[i], Z[1:(i-1)]])
    }
  }
  return(X)
}



#' Simulate of matrix of affiliation
#'
#' @param n An integer, the number of individuals
#' @param m An integer, the number of organisations
#' @param affiliation The type of affiliation between c("uniform", "preferential")
#' @param no_empty_org A Boolean. Force all column to have at least a 1. Need n>m
#'
#' @return A nxm affiliation matrix, with a unique 1 on each rows
#' @export
#'
#' @examples
simulate_affiliation <-
  function(n, m, affiliation = "uniform", no_empty_org = FALSE) {
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