#' Perform a spectral clustering
#'
#' @importFrom stats kmeans
#'
#' @param X an Adjacency matrix
#' @param K the number of clusters
#'
#' @return A vector : The clusters labels
spcClust <- function(X, K){
  if (K == 1) return (rep(1L, nrow(X)))
  n <- nrow(X)
  X[X == -1] <- NA
  isolated <- which(rowSums(X, na.rm = TRUE) == 0)
  connected <- setdiff(seq(n), isolated)
  X <- X[connected, connected]
  D_moins1_2 <- diag(1/sqrt(rowSums(X, na.rm = TRUE) + 1e-3))
  X[is.na(X)] <- mean(X, na.rm=TRUE)
  Labs <- D_moins1_2 %*% X %*% D_moins1_2
  specabs <- eigen(Labs, symmetric = TRUE)
  index <- rev(order(abs(specabs$values)))[1:K]
  U <- specabs$vectors[,index]
  U <- U / rowSums(U**2)**(1/2)
  U[is.na(U)] <- 0
  cl <- stats::kmeans(U, K, iter.max = 100, nstart=100)$cluster
  clustering <- rep(1, n)
  clustering[connected] <- cl
  clustering[isolated] <-   which.min(rowsum(rowSums(X, na.rm = TRUE),cl))
  return(clustering)
}

#' Perform a Hierarchical Clustering
#' @importFrom stats cutree dist hclust
#' @importFrom ape additive
#' @param X An Adjacency Matrix
#' @param K the number of wanted clusters
#'
#' @return A vector : The clusters labels
hierarClust <- function(X, K){
  if (K == 1) return (rep(1L, nrow(X)))
  # distance <- stats::dist(x = X, method = "manhattan")
  # X[X == -1] <- NA
  # distance[which(A == 1)] <- distance[which(A == 1)] - 2
  # distance <- stats::as.dist(ape::additive(distance))
  clust <- cluster::agnes(x = X, metric = "manhattan", method = "ward")
  # clust    <- stats::hclust(d = distance , method = "ward.D2")
  return(stats::cutree(tree = clust, k = K))
}

#' Merge a list of clusters
#'
#' @param X an adjacency matrix
#' @param Z a vector of cluster memberships
#' @param Q The number of maximal clusters
#'
#' @return A list of Q clustering of Q+1 clusters
split_clust <- function(X, Z, Q) {
  Z_split <-  lapply(
    X = seq(Q),
    FUN =  function(q) {
      if (sum(Z==q) < 2) return()
      Z_new        <-  Z
      Z_new[Z==q]  <-  hierarClust(X[Z==q, Z==q], 2) + Q
      Z_new[Z_new == Q + 2]  <-  q
      return(Z_new)
    })
  Z_split  <-  Z_split[which(sapply(X = Z_split, FUN = function(x) ! is.null(x)))]
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
merge_clust <- function(Z, Q) {
  Z_merge <- lapply(
    X = 1:choose(Q,2),
    FUN = function(q) {
      Z[Z == utils::combn(Q, 2)[2, q]] <- utils::combn(Q, 2)[1, q]
      if(utils::combn(Q, 2)[2, q] < Q){
        Z[Z > utils::combn(Q, 2)[2, q]] <- Z[Z > utils::combn(Q, 2)[2, q]] - 1
        }
      return(Z)
      }
    )
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

#' Compare two clustering with the Adjusted Rand Index
#'
#' @param x A vector of integers, the clusters labels
#' @param y A vector of integers of the same length as x, the clusters labels
#'
#' @return A number between 0 (random clustering) and 1 (identical clustering)
#' @export
#'
#' @examples ARI(x = c(1, 2, 1), y = c(2, 2, 1))
ARI <- function (x, y)
{
  x <- as.vector(x)
  y <- as.vector(y)
  if (length(x) != length(y))
    stop("arguments must be vectors of the same length")
  tab <- table(x, y)
  if (all(dim(tab) == c(1, 1)))
    return(1)
  a <- sum(choose(tab, 2))
  b <- sum(choose(rowSums(tab), 2)) - a
  c <- sum(choose(colSums(tab), 2)) - a
  d <- choose(sum(tab), 2) - a - b - c
  ARI <-
    (a - (a + b) * (a + c)/(a + b + c + d))/
    ((a + b + a + c)/2 - (a + b) * (a + c)/(a + b + c + d))
  return(ARI)
}


plot_multilevel_matrix <- function(X, X_hat, A, Z) {
  nodes <- group <- lvl <- edges <- weight <- NULL
  Z_sup <- c(Z$I, Z$O + max(Z$I))
  QI <- max(Z$I)
  QO <- max(Z$O)
  g_ind <- tidygraph::as_tbl_graph(t(X_hat$I * X$I ))  %>%
    tidygraph::activate(nodes) %>%
    tidygraph::mutate(group = Z$I) %>%
    tidygraph::activate(edges) %>%
    tidygraph::mutate(lvl = "ind")
  g_org <-
    tidygraph::as_tbl_graph(X_hat$O * X$O )  %>%
    tidygraph::activate(nodes) %>%
    tidygraph::mutate(group = Z$O + QI) %>%
    tidygraph::activate(edges) %>%
    tidygraph::mutate(lvl = "org")
  g_aff <-
    tidygraph::as_tbl_graph(A) %>%
    tidygraph::activate(edges) %>%
    tidygraph::mutate(lvl = "aff")

  p_mat <-  tidygraph::graph_join(g_ind, g_org) %>%
    tidygraph::graph_join(g_aff) %>%
    ggraph::ggraph('matrix', sort.by = group)+
    ggraph::geom_edge_point(ggplot2::aes(filter = (lvl == "aff")),
                            edge_colour = "black",  edge_size = 1.2)+
    ggraph::geom_edge_point(ggplot2::aes(filter = (lvl == "ind"),
                                         edge_colour = weight),edge_size = 1.2)+
    ggraph::geom_edge_point(ggplot2::aes(filter = (lvl == "org"), edge_fill = weight),
                    edge_size = 2, edge_shape = 22, stroke = 0)+
    ggplot2::geom_hline(yintercept = c( cumsum(table(Z_sup))[QI+seq(QO)]+.5)) +
    ggplot2::geom_hline(yintercept = cumsum(table(Z_sup))[QI]+.5, size = 1.1) +
    ggplot2::geom_vline(xintercept = cumsum(table(Z_sup))[QI]+.5, size = 1.1) +
    ggplot2::geom_vline(xintercept = c(0, cumsum(table(Z_sup))[seq(QI)]+.5)) +
    ggplot2::annotate(geom = "segment",
                      x = c(0, cumsum(table(Z_sup))[QI+seq(QO)]+.5),
                      y = cumsum(table(Z_sup))[QI]+.5,
                      xend = c(0,cumsum(table(Z_sup))[QI+seq(QO)]+.5),
                      yend = cumsum(table(Z_sup))[QI+QO]+.5) +
    ggplot2::annotate(geom = "segment",
                      y = c(0, cumsum(table(Z_sup))[seq(QI)]+.5),
                      x = cumsum(table(Z_sup))[QI]+.5,
                      yend = c(0, cumsum(table(Z_sup))[seq(QI)]+.5),
                      xend = 0) +
    ggraph::scale_edge_fill_gradient(
      name = 'Organizations',
      low = "#fcbba1",
      high = "#67000d",
      guide = ggraph::guide_edge_colorbar(order = 2,title.position = "top")) +
    ggraph::scale_edge_color_gradient(
      name = 'Individuals',
      low = "#deebf7",
      high = "#08519c",
      guide = ggraph::guide_edge_colorbar(order = 1,
                                          title.position = "top",
                                          title.hjust = 1)) +
    ggplot2::coord_fixed(xlim = c(-1, sum(dim(A)+1)), ylim = c(-1, sum(dim(A)+1))) +
    ggraph::theme_graph()
    return(p_mat)
}
