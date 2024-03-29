n <- 100
L <- 4
alpha <- list(
  diag(.4, 3, 3) + .1,
  -diag(.2, 3, 3) + .3,
  matrix(c(.8, .2, .1,
        .4, .4, .1,
        .2, .1, .1), 3, 3),
  matrix(c(.3, .5, .5,
         .1, .4, .5,
         .1, .3, .1), 3, 3)
)
alpha[[1]][1,1] <- .8
alpha[[1]][3,3] <- .2
gamma <- lapply(seq(3),
                function(m) matrix(c(.8, .1, .1,
                                     .1, .8, .1,
                                     .1, .1, .8), 3, 3, byrow = TRUE))
pi <- list(rep(1, 3)/3, NULL, c(.1, .3, .6), NULL)

directed = c(FALSE, FALSE, FALSE, TRUE)

gmlv <- mlvsbm_simulate_generalized_network(n = rep(n, 4),
                                    Q = rep(3, 4),
                                    pi = pi,
                                    gamma = gamma,
                                    alpha = alpha,
                                    directed = directed,
                                    distribution = rep("bernoulli", 4))

fit <- mlvsbm_estimate_generalized_network(gmlv, fit_options = list(ve = "joint"))
fit2 <- mlvsbm_estimate_generalized_network(gmlv, fit_options = list(ve = "sequential"))
fitone <- mlvsbm_estimate_generalized_network(gmlv2, nb_clusters = rep(3, 4),
                                              fit_options = list(ve = "sequential"))
fit_from_scratch <-  mlvsbm_estimate_generalized_network(gmlv,
                                    init_clustering = lapply(seq(4),
                                                             function(x)rep(1, n)),
                                    init_method = "merge_split",
                                    fit_options = list(ve = "joint"))
plot(fit)

mlvsbm_estimate_generalized_network(gmlv,
                                    init_clustering = lapply(seq(4),
                                                             function(x)rep(1, n)),
                                    init_method = "merge_split",
                                    fit_options = list(ve = "sequential"))
mlvsbm_estimate_generalized_network(gmlv,
                                    init_clustering = lapply(seq(4),
                                                             function(x)rep(1, n)),
                                    init_method = "merge_split",
                                    fit_options = list(ve = "joint"))



gmlv2 <- mlvsbm_create_generalized_network(X = gmlv$adjacency_matrix,
                                           A = gmlv$affiliation_matrix,
                                           directed = gmlv$directed,
                                           distribution = gmlv$distribution)



#==============================================================================
# Compare with dynsbm
#
library(dynsbm)
data(simdataT5Q4N40binary)
data("simdataT5Q4N40discrete")

dim(simdataT5Q4N40binary)
Xdyn <- lapply(seq(5), function(l) simdataT5Q4N40binary[l,,])
Adyn <- lapply(seq(4), function(l) diag(1, 40))

gdyn <- mlvsbm_create_generalized_network(X = Xdyn,
                                          A = Adyn,
                                          directed = rep(FALSE, 5),
                                          distribution = rep("bernoulli", 5))
fitdyn <- mlvsbm_estimate_generalized_network(gdyn)
fitdyn$parameters
plot(fitdyn)

Xpois <- lapply(seq(5), function(l) simdataT5Q4N40discrete[l,,])

gpois <- mlvsbm_create_generalized_network(X = Xpois,
                                          A = Adyn,
                                          directed = rep(FALSE, 5),
                                          distribution = rep("poisson", 5))
fitpois <- mlvsbm_estimate_generalized_network(gpois)
fitpois$parameters
plot(fitpois)

mlvsbm_estimate_generalized_network(gpois, init_clustering = lapply(seq(5), function(x)rep(1, 40)), init_method = "merge_split")

list.dynsbm <- select.dynsbm(simdataT5Q4N40binary,
                             Qmin=1, Qmax=5, edge.type="binary", nstart=25)
dynsbm <- list.dynsbm[[4]]

plot(fi)


## plotting intra/inter connectivity patterns
connectivity.plot(dynsbm, simdataT5Q4N40binary)


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
