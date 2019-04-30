#' Perform a spectral clustering
#'
#' @importFrom stats kmeans
#'
#' @param X an Adjacency matrix
#' @param K the number of clusters
#'
#' @return A vector : The clusters labels
#'
#'
#' @examples
spcClust <- function(X, K){
  n <- nrow(X)
  D_moins1_2 <- diag(1/sqrt(colSums(X, na.rm = TRUE) + 1e-3))
  X[is.na(X)] <- mean(X, na.rm=TRUE)
  Labs <- D_moins1_2 %*% X %*% D_moins1_2
  specabs <- eigen(Labs)
  index <- order(abs(specabs$values), decreasing = FALSE)[(n-K+1):n]
  U <- specabs$vectors[,index]
  U <- U / rowSums(U**2)**(1/2)
  U[is.na(U)] <- 0
  clustering <- stats::kmeans(U, K, nstart=100)$cluster
  return(clustering)
}

#' Perform a Hierarchical Clustering
#' @importFrom stats cutree dist hclust
#' @param X An Adjacency Matrix
#' @param K the number of wanted clusters
#'
#' @return A vector : The clusters labels
#'
#'
#' @examples
hierarClust <- function(X, K){
  distance <- stats::dist(x = X, method = "manhattan")
  # distance[which(A == 1)] <- distance[which(A == 1)] - 2
  # distance <- as.dist(ape::additive(distance))
  clust    <- stats::hclust(d = distance , method = "ward.D")
  return(stats::cutree(tree = clust, k = K))
}

#' Merge a list of clusters
#'
#' @param X an adjacency matrix
#' @param Z a vector of cluster memberships
#' @param Q The number of maximal clusters
#'
#' @return A list of Q clustering of Q+1 clusters
#'
#'
#' @examples
split_clust <- function(X, Z, Q) {
  Z_split <-  lapply(
    seq(Q),
    q <- function(q) {
      if (sum(Z==q) < 2) return()
      Z_new        <-  Z
      Z_new[Z==q]  <-  hierarClust(X[Z==q, Z==q], 2)+Q
      Z_new[Z_new == Q + 2]  <-  q
      return(Z_new)
    })
  Z_split  <-  Z_split[which(sapply(Z_split, x <- function(x) ! is.null(x)))]
  return(Z_split)
}

#' Merge a list of clusters
#'
#' @importFrom utils combn
#'
#' @param Z a vector of cluster memberships
#' @param Q the number of original clusters
#'
#' @return A list of Q(Q-1)/2 clustering of Q-1 clusters
#'
#'
#' @examples
merge_clust <- function(Z, Q) {
  Z_merge = lapply(X = 1:choose(Q,2),
                   FUN = function(q) {
                     Z[Z == utils::combn(Q, 2)[2, q]] = utils::combn(Q, 2)[1, q]
                     if(utils::combn(Q, 2)[2, q] < Q){
                       Z[Z > utils::combn(Q, 2)[2, q]] = Z[Z > utils::combn(Q, 2)[2, q]] - 1
                     }
                     return(Z)
                   })
  return(Z_merge)
}

#
#Fonction interne a VEM
#

F.bern <- function(X, alpha, tau){
  return(X %*% tau %*% log(alpha) +
           (1 - X - diag(1, nrow(X))) %*% tau %*% log(1-alpha))
}

rotate <- function(x) t(apply(x, 2, rev))

#
dist_param <- function(param, param_old) {
  sqrt(sum((param-param_old)**2))
}

#' Title
#'
#' @param X An adjacency matrix
#' @param K An integer, the number of folds
#'
#' @return A matrix of the same size than X with class integer as coefficient
#'
#' @examples
build_fold_matrix <- function(X, K) {
  n <- ncol(X)
  arrange     <- sample(x = seq(n))
  labels      <- cut(seq(n), breaks = K, labels = FALSE)
  fold_matrix <- diag(n)
  for (i in seq(n)) {
    fold_matrix[i, ] <- (labels + labels[i]) %% K
  }
  fold_matrix <- fold_matrix[arrange, arrange] + 1
  diag(fold_matrix) <- 0
  return(fold_matrix)
}


xlogx      <- function(x) ifelse(x < 2*.Machine$double.eps, 0, x*log(x))
quad_form  <- function(X, tau) tau %*% tcrossprod(X, tau)
logistic   <- function(x) 1/(1 + exp(-x))
logit      <- function(x) log(x/(1 - x))
