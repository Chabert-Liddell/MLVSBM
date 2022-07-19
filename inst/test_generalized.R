n <- 100
L <- 4
alpha <- list(
  diag(.4, 3, 3) + .1,
  -diag(.3, 3, 3) + .5,
  matrix(c(.8, .2, .1,
        .4, .4, .1,
        .2, .1, .1), 3, 3),
  matrix(c(.3, .5, .5,
         .1, .4, .5,
         .1, .3, .1), 3, 3)
)
gamma <- lapply(seq(3),
                function(m) matrix(c(.7, .2, .1,
                                     .1, .7, .2,
                                     .2, .1, .7), 3, 3))
pi <- list(rep(1, 3)/3, NULL, c(.1, .3, .6), NULL)

directed = c(FALSE, FALSE, FALSE, TRUE)

X <- list()
Z <- list()

Z[[1]] <- c(rep(1, 40), rep(2, 30), rep(3, 30))
X[[1]] <- MLVSBM:::simulate_adjacency(Z = Z[[1]],
                                      n = 100,
                                      alpha = alpha[[1]], directed = FALSE,
                                      distribution = "bernoulli")
A <- lapply(seq(3), function(m) diag(1, 100))
A[[2]][81:100,] <- 0

Z[[4]] <- Z[[3]] <- Z[[2]] <- rep(0, 100)
#ind          <- A[[1]] %*% Z[[1]]
for (m in seq(2,4)) {
  for (i in seq(n)) {
    ind          <- A[[m-1]] %*% Z[[m-1]]
    if (ind[i,] > 0) {
      Z[[m]][i] <- sample(seq(3),
                          size = 1,
                          prob = gamma[[m-1]][, ind[i, ]])
    } else {
      Z[[m]][i] <- sample(seq(3),
                          size = 1,
                          prob = pi[[m]])
    }

  }
  X[[m]] <- MLVSBM:::simulate_adjacency(Z = Z[[m]],
                              n = 100,
                              alpha = alpha[[m]],
                              directed = directed[m],
                              distribution = "bernoulli")
}

#A[[2]] <- diag(1, 100)

fit <- MLVSBM:::FitGenMLVSBM$new(Q = rep(3, 4), A = A, X = X, directed = directed,
                                 no_affiliation = c(FALSE, FALSE, TRUE, FALSE),
                                 distribution = rep("bernoulli", 4))

fit$init_clustering()
fit$do_vem()
#fit$initialize()

genmlv <- MLVSBM:::mlvsbm_create_generalized_network(X, A, directed,
                                                     rep("bernoulli", 4))

fit <- MLVSBM:::mlvsbm_estimate_generalized_network(genmlv)

myGenMLVSBM <- MLVSBM:::GenMLVSBM$new(n = n, X = X, A = A, L = 4,
                                      directed = directed,
                                      distribution = rep("bernoulli", 4))


myGenMLVSBM


fitone <- myGenMLVSBM$estimate_one(Q = rep(3, 4))
fitall <- myGenMLVSBM$estimate_all_bm(Q = rep(1,4), Z = lapply(seq(4), function(m) rep(1, 100)))

fitsbm <- list()
for (m in seq(4)) {
  fitsbm[[m]] <- myGenMLVSBM$estimate_level(level = m)
}
fitcomplete <- myGenMLVSBM$estimate_all_bm(Q = rep(3,4),
                            Z =lapply(seq(4), function(m) myGenMLVSBM$fittedmodels_sbm[[m]][[3]]$Z))

(sapply(seq(4), function(m) ARI(Z[[m]], fit$Z[[m]])))
(sapply(seq(4), function(m) ARI(Z[[m]], fitall$Z[[m]])))
(sapply(seq(4), function(m) ARI(Z[[m]], fitone$Z[[m]])))
(sapply(seq(4), function(m) ARI(Z[[m]], fitcomplete$Z[[m]])))

#===============================================================================
# Testing 2 levels de la vignette
#===============================================================================

X2 <- list(upper_level, lower_level)
A <- list(affiliation)


myGenMLVSBM <- MLVSBM:::GenMLVSBM$new(n = c(40, 60), X = X2, A = A,
                                      L = 2, directed = c(FALSE, TRUE),
                                      distribution = rep("bernoulli", 2))
fitall <- myGenMLVSBM$estimate_all_bm(Q = rep(1,2),
                                      Z = list(rep(1, 40), rep(1, 60)))


#===============================================================================
# hard to infer vignette
#===============================================================================

myGenMLVSBM <- MLVSBM:::GenMLVSBM$new(n = unlist(n), X = mlvl$adjacency_matrix[c(2,1)],
                                      A = list(mlvl$affiliation_matrix),
                                      L = 2, directed = c(FALSE, FALSE),
                                      distribution = rep("bernoulli", 2))
fitall <- myGenMLVSBM$estimate_all_bm(Q = rep(1,2),
                                      Z = list(rep(1, 40), rep(1, 120)))
for (m in seq(2)) {
  myGenMLVSBM$estimate_level(level = m)
}
fitcomplete <- myGenMLVSBM$estimate_all_bm( Q = c(1, 2),
  Z =list(myGenMLVSBM$fittedmodels_sbm[[1]][[1]]$Z,
            myGenMLVSBM$fittedmodels_sbm[[2]][[2]]$Z))
fitcheat <- myGenMLVSBM$estimate_all_bm( Q = c(2, 2),
                                         Z = mlvl$fittedmodels[[2]]$Z[c(2,1)])
