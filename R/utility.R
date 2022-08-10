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
  if (! isSymmetric(X)) {
    X <- 1*((X+t(X)) >0)
  }
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
  diss <- cluster::daisy(x = X, metric = "manhattan", warnType = FALSE)
  if (! any(is.na(diss))) {
    clust <- cluster::agnes(x = X, metric = "manhattan", method = "ward")
  } else {
    return (rep(1L, nrow(X)))
  }
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


.xlogx      <- function(x) ifelse(x < 2*.Machine$double.eps, 0, x*log(x))

quad_form  <- function(X, tau) tau %*% tcrossprod(X, tau)
logistic   <- function(x) 1/(1 + exp(-x))
logit      <- function(x) log(x/(1 - x))



.xlogy     <- function(x, y, eps = NULL) {
  ifelse(x < 2*.Machine$double.eps, 0, x*.log(y, eps = eps))
}
.quadform  <- function(x, y)  tcrossprod(x %*% y, x)
.tquadform  <- function(x, y)  crossprod(x, y %*% x)
logistic   <- function(x) 1/(1 + exp(-x))
logit      <- function(x) log(x/(1 - x))
.logit <- function(x, eps = NULL) {
  if(is.null(eps)) {
    res <- log(x/(1-x))
  } else {
    res <- log(pmax(pmin(x, 1-eps), eps)/pmax(pmin(1-x, 1-eps), eps))
  }
  return (res)
}


.threshold <- function(x, eps = 1e-9) {
  #  x <- .softmax(x)
  x[x < eps] <- eps
  x[x > 1-eps] <- 1-eps
  x <- x/.rowSums(x, nrow(x), ncol(x))
  x
}

.softmax <- function(x) {
  x_max <- apply(x, 1, max)
  x <- exp(x - x_max)
  x <- x/.rowSums(x, nrow(x), ncol(x))
  x
}
.log <- function(x, eps = NULL) {
  if(is.null(eps)) {
    res <- log(x)
  } else {
    res <- log(pmax(pmin(x, 1-eps), eps))
  }
  return (res)
}
.one_hot <- function(x, Q) {
  O <- matrix(0, length(x),Q)
  O[cbind(seq.int(length(x)), x)] <- 1
  return(O)
}
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


# plot_multilevel_matrix <- function(X, X_hat, A, Z) {
#   nodes <- group <- lvl <- edges <- weight <- NULL
#   Z_sup <- c(Z$I, Z$O + max(Z$I))
#   QI <- max(Z$I)
#   QO <- max(Z$O)
#   g_ind <- tidygraph::as_tbl_graph(t(X_hat$I * X$I ))  %>%
#     tidygraph::activate(nodes) %>%
#     tidygraph::mutate(group = Z$I) %>%
#     tidygraph::activate(edges) %>%
#     tidygraph::mutate(lvl = "ind")
#   g_org <-
#     tidygraph::as_tbl_graph(X_hat$O * X$O )  %>%
#     tidygraph::activate(nodes) %>%
#     tidygraph::mutate(group = Z$O + QI) %>%
#     tidygraph::activate(edges) %>%
#     tidygraph::mutate(lvl = "org")
#   g_aff <-
#     tidygraph::as_tbl_graph(A) %>%
#     tidygraph::activate(edges) %>%
#     tidygraph::mutate(lvl = "aff")
#
#   p_mat <-  tidygraph::graph_join(g_ind, g_org) %>%
#     tidygraph::graph_join(g_aff) %>%
#     ggraph::ggraph('matrix', sort.by = group)+
#     ggraph::geom_edge_point(ggplot2::aes(filter = (lvl == "aff")),
#                             edge_colour = "black",  edge_size = 1.2)+
#     ggraph::geom_edge_point(ggplot2::aes(filter = (lvl == "ind"),
#                                          edge_colour = weight),edge_size = 1.2)+
#     ggraph::geom_edge_point(ggplot2::aes(filter = (lvl == "org"), edge_fill = weight),
#                     edge_size = 2, edge_shape = 22, stroke = 0)+
#     ggplot2::geom_hline(yintercept = c( cumsum(table(Z_sup))[QI+seq(QO)]+.5)) +
#     ggplot2::geom_hline(yintercept = cumsum(table(Z_sup))[QI]+.5, size = 1.1) +
#     ggplot2::geom_vline(xintercept = cumsum(table(Z_sup))[QI]+.5, size = 1.1) +
#     ggplot2::geom_vline(xintercept = c(0, cumsum(table(Z_sup))[seq(QI)]+.5)) +
#     ggplot2::annotate(geom = "segment",
#                       x = c(0, cumsum(table(Z_sup))[QI+seq(QO)]+.5),
#                       y = cumsum(table(Z_sup))[QI]+.5,
#                       xend = c(0,cumsum(table(Z_sup))[QI+seq(QO)]+.5),
#                       yend = cumsum(table(Z_sup))[QI+QO]+.5) +
#     ggplot2::annotate(geom = "segment",
#                       y = c(0, cumsum(table(Z_sup))[seq(QI)]+.5),
#                       x = cumsum(table(Z_sup))[QI]+.5,
#                       yend = c(0, cumsum(table(Z_sup))[seq(QI)]+.5),
#                       xend = 0) +
#     ggraph::scale_edge_fill_gradient(
#       name = 'Organizations',
#       low = "#fcbba1",
#       high = "#67000d",
#       guide = ggraph::guide_edge_colorbar(order = 2,title.position = "top")) +
#     ggraph::scale_edge_color_gradient(
#       name = 'Individuals',
#       low = "#deebf7",
#       high = "#08519c",
#       guide = ggraph::guide_edge_colorbar(order = 1,
#                                           title.position = "top",
#                                           title.hjust = 1)) +
#     ggplot2::coord_fixed(xlim = c(-1, sum(dim(A)+1)), ylim = c(-1, sum(dim(A)+1))) +
#     ggraph::theme_graph()
#     return(p_mat)
# }


plot_multilevel_graphon <- function(fit, order = "degree") {
  xmin <- xmax <- ymin <- ymax <- value <- NULL
  ord <- list()
  ord$O <- seq(fit$nb_clusters$O)
  ord$I <- seq(fit$nb_clusters$I)
  if (order == "degree") {
    if(fit$nb_clusters$O > 1) {
      ord$O <- order(fit$parameters$pi$O %*% fit$parameters$alpha$O, decreasing=TRUE)
    }
    if(fit$nb_clusters$I > 1) {
      ord$I <- order( t(fit$parameters$gamma %*% fit$parameters$pi$O) %*% fit$parameters$alpha$I, decreasing=TRUE)
    }
  }
  if (order == "affiliation") {
    if(fit$nb_clusters$O > 1) {
      ord$O <- order(fit$parameters$pi$O %*% fit$parameters$alpha$O, decreasing=TRUE)
    }
    if(fit$nb_clusters$I > 1) {
    ord$I <- order(apply(fit$parameters$gamma, 1, "which.max"), decreasing = TRUE)
    }
  }

  p <- list()
  p$O <- fit$parameters$alpha$O[ord$O, ord$O] %>% t() %>%
    reshape2::melt() %>%
    dplyr::mutate(xmax = rep(c(0,cumsum(fit$parameters$pi$O[ord$O][1:(fit$nb_clusters$O-1)])), fit$nb_clusters$O),
                  xmin = rep(cumsum(fit$parameters$pi$O[ord$O]), fit$nb_clusters$O),
                  ymax = rep(c(0,cumsum(fit$parameters$pi$O[ord$O][1:(fit$nb_clusters$O-1)])), each = fit$nb_clusters$O),
                  ymin = rep(cumsum(fit$parameters$pi$O[ord$O]), each = fit$nb_clusters$O)) %>%
    ggplot2::ggplot(ggplot2::aes(xmin = xmin, ymin = ymin, xmax = xmax, ymax = ymax, fill = value)) +
    ggplot2::geom_rect() +
    ggplot2::scale_fill_gradient2("Org", low = "white", mid = "red", midpoint = 1) +
    ggplot2::geom_hline(yintercept = cumsum(fit$parameters$pi$O[ord$O][1:(fit$nb_clusters$O-1)]), size = .2) +
    ggplot2::geom_vline(xintercept = cumsum(fit$parameters$pi$O[ord$O][1:(fit$nb_clusters$O-1)]), size = .2) +
    ggplot2::scale_y_reverse() +
    ggplot2::theme_void(base_size = 15, base_rect_size = 1, base_line_size  = 1) +
    ggplot2::xlab("") + ggplot2::ylab("") + ggplot2::coord_equal(expand = FALSE)

  fit$parameters$pi$I <-
    as.vector(fit$parameters$gamma %*% fit$parameters$pi$O)

  p$I <- fit$parameters$alpha$I[ord$I, ord$I] %>% t() %>%
    reshape2::melt() %>%
    dplyr::mutate(xmax = rep(c(0,cumsum(fit$parameters$pi$I[ord$I][1:(fit$nb_clusters$I-1)])), fit$nb_clusters$I),
                  xmin = rep(cumsum(fit$parameters$pi$I[ord$I]), fit$nb_clusters$I),
                  ymax = rep(c(0,cumsum(fit$parameters$pi$I[ord$I][1:(fit$nb_clusters$I-1)])), each = fit$nb_clusters$I),
                  ymin = rep(cumsum(fit$parameters$pi$I[ord$I]), each = fit$nb_clusters$I)) %>%
    ggplot2::ggplot(ggplot2::aes(xmin = xmin, ymin = ymin, xmax = xmax, ymax = ymax, fill = value)) +
    ggplot2::geom_rect() +
    ggplot2::scale_fill_gradient2("Ind", low = "white", mid = "blue", midpoint = 1) +
    ggplot2::geom_hline(yintercept = cumsum(fit$parameters$pi$I[ord$I][1:(fit$nb_clusters$I-1)]), size = .2) +
    ggplot2::geom_vline(xintercept = cumsum(fit$parameters$pi$I[ord$I][1:(fit$nb_clusters$I-1)]), size = .2) +
    ggplot2::scale_y_reverse() +
    ggplot2::theme_void(base_size = 15, base_rect_size = 1, base_line_size  = 1) +
    ggplot2::xlab("") + ggplot2::ylab("") + ggplot2::coord_equal(expand = FALSE)

  p$A <- fit$parameters$gamma[ord$I, ord$O] %>% t() %>%
    reshape2::melt() %>%
    dplyr::mutate(xmax = rep(c(0,cumsum(fit$parameters$pi$O[ord$O][1:(fit$nb_clusters$O-1)])), fit$nb_clusters$I),
                  xmin = rep(cumsum(fit$parameters$pi$O[ord$O]), fit$nb_clusters$I),
                  ymax = rep(c(0,cumsum(fit$parameters$pi$I[ord$I][1:(fit$nb_clusters$I-1)])), each = fit$nb_clusters$O),
                  ymin = rep(cumsum(fit$parameters$pi$I[ord$I]), each = fit$nb_clusters$O)) %>%
    ggplot2::ggplot(ggplot2::aes(xmin = xmin, ymin = ymin, xmax = xmax, ymax = ymax, fill = value)) +
    ggplot2::geom_rect() +
    ggplot2::scale_fill_gradient2("Aff", low = "white", mid = "black", midpoint = 1) +
    ggplot2::geom_hline(yintercept = cumsum(fit$parameters$pi$I[ord$I][1:(fit$nb_clusters$I-1)]), size = .2) +
    ggplot2::geom_vline(xintercept = cumsum(fit$parameters$pi$O[ord$O][1:(fit$nb_clusters$O-1)]), size = .2) +
    ggplot2::scale_y_reverse() +
    ggplot2::theme_void(base_size = 15, base_rect_size = 1, base_line_size  = 1) +
    ggplot2::xlab("") + ggplot2::ylab("") + ggplot2::coord_equal(expand = FALSE)

  leg <- list()
  leg$O <- cowplot::get_legend(p$O)
  p$O <- p$O + ggplot2::theme(legend.position = "none")
  leg$I <- cowplot::get_legend(p$I)
  p$I <- p$I + ggplot2::theme(legend.position = "none")
  leg$A <- cowplot::get_legend(p$A)
  p$A <- p$A + ggplot2::theme(legend.position = "none")

  p_mat <- cowplot::ggdraw()+
    cowplot::draw_plot(p$I, x = 0, y = 0, width = .5, height = .5) +
    cowplot::draw_plot(p$O, x = 0.5, y = .5, width = .5, height = .5) +
    cowplot::draw_plot(p$A, x = .5, y = 0, width = .5, height = .5) +
    ggplot2::geom_hline(yintercept = .5) +
    ggplot2::geom_vline(xintercept = .5) +
    cowplot::draw_plot(leg$I, x = .0, y = .7, width = .2, height = .2 ) +
    cowplot::draw_plot(leg$A, x = .15, y = .7, width = .2, height = .2 ) +
    cowplot::draw_plot(leg$O, x = .3, y = .7, width = .2, height = .2 ) +
    ggplot2::coord_equal(xlim = c(0,1), ylim = c(0,1), )
  return(p_mat)
}





plot_generalized_multilevel_graphon <- function(fit, order = "affiliation") {
  # browser()
  xmin <- xmax <- ymin <- ymax <- value <- NULL
  color <- RColorBrewer::brewer.pal(min(9, max(3, fit$nb_levels)), name = "Set1")
  if (length(fit$nb_levels) > length(color)) {
    color <- rep(fit$nb_levels, "blue")
  }
  fit$reorder(order = order)
  p <- list()
  for(l in seq(fit$nb_levels)) {
    if (fit$nb_clusters[l] > 1) {
      p[[l]] <- fit$parameters$alpha[[l]] %>% t() %>%
        reshape2::melt() %>%
        dplyr::mutate(xmax = rep(c(0,cumsum(fit$block_proportions[[l]][1:(fit$nb_clusters[l]-1)])),
                                 fit$nb_clusters[l]),
                      xmin = rep(cumsum(fit$block_proportions[[l]]), fit$nb_clusters[l]),
                      ymax = rep(c(0,cumsum(fit$block_proportions[[l]][1:(fit$nb_clusters[l]-1)])),
                                 each = fit$nb_clusters[l]),
                      ymin = rep(cumsum(fit$block_proportions[[l]]), each = fit$nb_clusters[l])) %>%
        ggplot2::ggplot(ggplot2::aes(xmin = xmin, ymin = ymin, xmax = xmax, ymax = ymax, fill = value)) +
        ggplot2::geom_rect() +
        ggplot2::scale_fill_gradient2(low = "white", mid = color[l], midpoint = 1) +
        ggplot2::geom_hline(yintercept = cumsum(fit$block_proportions[[l]][1:(fit$nb_clusters[l]-1)]), size = .2) +
        ggplot2::geom_vline(xintercept = cumsum(fit$block_proportions[[l]][1:(fit$nb_clusters[l]-1)]), size = .2) +
        ggplot2::geom_hline(yintercept = c(0, 1), size = 1.1) +
        ggplot2::geom_vline(xintercept = c(0, 1), size = 1.1) +
        ggplot2::scale_y_reverse() +
        ggplot2::theme_void(base_size = 15, base_rect_size = 1, base_line_size  = 1) +
        ggplot2::xlab("") + ggplot2::ylab("") + ggplot2::coord_equal(expand = FALSE)
    } else {
      p[[l]] <- fit$parameters$alpha[[l]] %>% t() %>%
        reshape2::melt() %>%
        ggplot2::ggplot(ggplot2::aes(xmin = 0, ymin = 0, xmax = 1, ymax = 1, fill = value)) +
        ggplot2::geom_rect() +
        ggplot2::scale_fill_gradient2(low = "white", mid = color[l], midpoint = 1) +
        ggplot2::geom_hline(yintercept = cumsum(fit$block_proportions[[l]][1:(fit$nb_clusters[l]-1)]), size = .2) +
        ggplot2::geom_vline(xintercept = cumsum(fit$block_proportions[[l]][1:(fit$nb_clusters[l]-1)]), size = .2) +
        ggplot2::geom_hline(yintercept = c(0, 1), size = 1.1) +
        ggplot2::geom_vline(xintercept = c(0, 1), size = 1.1) +
        ggplot2::scale_y_reverse() +
        ggplot2::theme_void(base_size = 15, base_rect_size = 1, base_line_size  = 1) +
        ggplot2::xlab("") + ggplot2::ylab("") + ggplot2::coord_equal(expand = FALSE)
    }

  }

  pA <- list()

  for (l in seq(fit$nb_levels -1)) {
    if(fit$nb_clusters[l] == 1) {
      xmin <- rep(1, fit$nb_clusters[l+1])
      xmax <- rep(0, fit$nb_clusters[l+1])
    } else {
      xmax <-  rep(c(0,cumsum(fit$block_proportions[[l]][1:(fit$nb_clusters[l]-1)])),
                 fit$nb_clusters[l+1])
      xmin <-  rep(cumsum(fit$block_proportions[[l]]), fit$nb_clusters[l+1])
    }
    if(fit$nb_clusters[l+1] == 1) {
      ymin <- rep(1, fit$nb_clusters[l])
      ymax <- rep(0, fit$nb_clusters[l])
    } else {
      ymax <-  rep(c(0,cumsum(fit$block_proportions[[l+1]][1:(fit$nb_clusters[l+1]-1)])),
                 each = fit$nb_clusters[l])
      ymin <-  rep(cumsum(fit$block_proportions[[l+1]]), each = fit$nb_clusters[l])
    }
    pA[[l]] <- fit$parameters$gamma[[l]] %>% t() %>%
      reshape2::melt() %>%
      dplyr::mutate(xmax = xmax,
                    xmin = xmin,
                    ymax = ymax,
                    ymin = ymin) %>%
      ggplot2::ggplot(ggplot2::aes(xmin = xmin, ymin = ymin, xmax = xmax, ymax = ymax, fill = value)) +
      ggplot2::geom_rect() +
      ggplot2::scale_fill_gradient2("Aff", low = "white", mid = "black", midpoint = 1) +
      ggplot2::geom_hline(yintercept = cumsum(fit$block_proportions[[l+1]][1:(fit$nb_clusters[l+1]-1)]),
                          size = .2) +
      ggplot2::geom_vline(xintercept = cumsum(fit$block_proportions[[l]][1:(fit$nb_clusters[l]-1)]),
                          size = .2) +
      ggplot2::geom_hline(yintercept = c(0, 1), size = 1.1) +
      ggplot2::geom_vline(xintercept = c(0, 1), size = 1.1) +
      ggplot2::scale_y_reverse() +
      ggplot2::theme_void(base_size = 15, base_rect_size = 1, base_line_size  = 1) +
      ggplot2::xlab("") + ggplot2::ylab("") + ggplot2::coord_equal(expand = FALSE)
  }

  leg <- list()
  for (l in seq(fit$nb_levels)) {
    leg[[l]] <- cowplot::get_legend(p[[l]])
    p[[l]] <- p[[l]] + ggplot2::theme(legend.position = "none") +
      ggplot2::labs(tag = l)
  }
  legA <- list()
  for (l in seq(fit$nb_levels-1)) {
    legA[[l]] <- cowplot::get_legend(pA[[l]])
    pA[[l]] <- pA[[l]] + ggplot2::theme(legend.position = "none") +
      ggplot2::labs(tag = paste0(l, "-", l+1))
  }
  pl <- vector("list", 2*fit$nb_levels)
  if (fit$nb_levels%%2 == 0) {
    idl <- sort(c(seq(1, 2*fit$nb_levels-1-2, by = 4),
                  seq(4, 2*fit$nb_levels-1-2, by = 4)))
  } else {
    idl <- sort(c(seq(1, 2*fit$nb_levels-1-1, by = 4),
                  seq(4, 2*fit$nb_levels-1-1, by = 4)))
  }

  pbl <- ggplot2::ggplot() + ggplot2::theme_void()


  pl[setdiff(seq(2*fit$nb_levels-1), idl)] <- p
  pl[idl] <- pA
  if (fit$nb_levels%%2 != 0) {
    pl[[2*fit$nb_levels-1]] <- pbl
    pl[[2*fit$nb_levels]] <- p[[fit$nb_levels]]
  } else {
    pl[[2*fit$nb_levels]] <- pbl
  }
  p_mat <- do.call(eval(parse(text="patchwork::wrap_plots")),
                   c(pl, ncol = fit$nb_levels, byrow=FALSE))
  # p_mat <- do.call(eval(parse(text="gridExtra::grid.arrange")),
  #         c(pl, ncol = fit$nb_levels, as.table=FALSE))
  return(p_mat)
}
