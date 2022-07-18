n <- 100
L <- 4
alpha <- list(
  diag(.4, 3, 3) + .1,
  -diag(.3, 3, 3) + .5,
  matrix(c(.7, .4, .1,
        .4, .4, .1,
        .4, .1, .1), 3, 3),
  matrix(c(.3, .5, .5,
         .1, .4, .5,
         .1, .3, .1), 3, 3)
)
gamma <- lapply(seq(3),
                function(m) matrix(c(.8, .1, .1,
                                     .1, .8, .1,
                                     .1, .1, .8), 3, 3))

directed = c(FALSE, FALSE, FALSE, TRUE)

X <- list()
Z <- list()
Z[[1]] <- c(rep(1, 40), rep(2, 30), rep(3, 30))
X[[1]] <- MLVSBM:::simulate_adjacency(Z = Z[[1]],
                                      n = 100,
                                      alpha = alpha[[1]], directed = FALSE,
                                      distribution = "bernoulli")
A <- lapply(seq(3), function(m) diag(1, 100))
Z[[4]] <- Z[[3]] <- Z[[2]] <- rep(0, 100)
#ind          <- A[[1]] %*% Z[[1]]
for (m in seq(2, 4)) {
  for (i in seq(n)) {
    ind          <- A[[m-1]] %*% Z[[m-1]]
    Z[[m]][i] <- sample(seq(3),
                        size = 1,
                        prob = gamma[[m-1]][, ind[i, ]])
  }
  X[[m]] <- MLVSBM:::simulate_adjacency(Z = Z[[m]],
                              n = 100,
                              alpha = alpha[[m]],
                              directed = directed[m],
                              distribution = "bernoulli")
}

fit <- MLVSBM:::FitGenMLVSBM$new(Q = rep(3, 4), A = A, X = X, directed = directed,
                                 no_affiliation = rep(FALSE,4),
                                 distribution = rep("bernoulli", 4))

fit$init_clustering()
fit$do_vem()
fit$initialize()
