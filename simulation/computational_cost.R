library(tidyverse)
estimate_computation_cost <-
  function(n_iter = 30, n = NA, Q = list(I = 3, O = 3),
           alpha = NA, gamma = NA, pi = NA,
           directed = list(I = FALSE, O = FALSE),
           affiliation = "preferential", nb_cores = 1L) {
  result <- pbmcapply::pbmclapply(
    X = seq(n_iter),
    FUN = function(i) {
      mlvl <-
        MLVSBM::mlvsbm_simulate_network(
          n = n, Q = Q, pi = pi, gamma = gamma, alpha = alpha,
          directed = directed, affiliation = affiliation)
      t0 <- Sys.time()
      fit <- MLVSBM::mlvsbm_estimate_network(mlvl, nb_cores = nb_cores)
      t1 <- Sys.time()
      return(difftime(t1, t0, units = "secs"))
    }, mc.cores = floor(12/nb_cores))
  return(result)
}

Q     <- list(I = 2, O = 2)
n     <- list(I = 50*Q$I, O = 50*Q$O)
alpha <- list()
alpha$O <- .1 * (diag(x = 1, nrow = Q$O, ncol = Q$O) + 1)
pi      <- rep(1/Q$O, Q$O)
gamma   <- .1 * (diag(8, Q$I, Q$O) + 2/Q$I)
alpha$I <- .1 * (diag(x = -2, nrow = Q$I, ncol = Q$I) + 2)
# alpha$I <- .1 * matrix(c(5, 3, 1,
#                        3, 1, 1,
#                        1, 1, 1), ncol = Q$I, nrow = Q$I)






mlvl <-
  MLVSBM::mlvsbm_simulate_network(
    n = n, Q = Q, pi = pi, gamma = gamma, alpha = alpha,
    directed = list(I = FALSE, O = FALSE))
sum(colSums(mlvl$adjacency_matrix$O)==0) +sum(colSums(mlvl$adjacency_matrix$I)==0)
t0 <- Sys.time()
fit <- MLVSBM::mlvsbm_estimate_network(mlvl, nb_cores = 2L)
t1 <- Sys.time()
difftime(t1, t0)

res_cost <- tibble()

for (nb_cores in c(2L)) {
  for(nO in c(500)) {
    for(nI in c(3*nO)) {
      for(K in c(2, 4,  8)) {
        Q     <- list(I = K, O = K)
        n     <- list(I = nI, O = nO)
        alpha <- list()
        alpha$O <- Q$O*sqrt(Q$O) * .01 * (diag(x = 1, nrow = Q$O, ncol = Q$O) + 1)
        pi      <- rep(1/Q$O, Q$O)
        gamma   <- .1 * (diag(8, Q$I, Q$O) + 2/Q$I)
        alpha$I <- Q$I*sqrt(Q$I) * .01 * (diag(x = -2, nrow = Q$I, ncol = Q$I) + 3)
        set.seed(1)
        time_cost <- estimate_computation_cost(n_iter = 30,
                                                 n = n,
                                                 Q = Q,
                                                 alpha = alpha,
                                                 gamma = gamma,
                                                 pi = pi, nb_cores = nb_cores)
        res_cost <-
          bind_rows(res_cost, tibble(nI = nI,
                           nO = nO,
                           Q = K,
                           nb_cores = nb_cores,
                           time = unlist(time_cost)))
        saveRDS(res_cost, file = "~/Documents/Multilevel/Simulations/article/comp_cost.rds")
      }
    }
  }
}
res_cost <- readRDS(file = "~/Documents/Multilevel/Simulations/article/comp_cost.rds")
res_cost %>%
  ggplot(aes(x = as.factor(Q), y = time)) +
  geom_boxplot() +
  facet_grid(nO + nI ~ .)

res_cost %>% group_by(nO, nI, Q) %>%
  summarise(q10 =  quantile(time, probs = .1),
            q50 = quantile(time, probs = .5),
            q90 = quantile(time, probs = .9))



################################################################################
# Near detectability threshold
################################################################################

Q     <- list(I = 2, O = 2)
n     <- list(I = 60*Q$I, O = 20*Q$O)
alpha <- list()
alpha$I <- matrix(c(.1, .21, .21, .1), 2, 2)
alpha$O <- matrix(c(.27, .1, .1, .27), 2, 2)
pi      <- rep(1/Q$O, Q$O)
gamma   <- .1 * (diag(8, Q$I, Q$O) + 2/Q$I)



detect_threshold <- function(n, q) {
  polyroot(c(n*q*(n*q-2),-2*n*(n*q+1),n**2 ))
}
detect_threshold(40, .1)
detect_threshold(120, .1)


res_detect_thresh <- tibble()
for (g in seq(.75, 1, .05)) {
  gamma <-   matrix(c(g, 1-g, 1-g, g), 2, 2)
  res <- pbmcapply::pbmclapply(
    X = seq(5),
    FUN = function(i) {
      mlvl <-
        MLVSBM::mlvsbm_simulate_network(
          n = n, Q = Q, pi = pi, gamma = gamma, alpha = alpha,
          directed = list(I = FALSE, O = FALSE))
      sbm_lower <- blockmodels::BM_bernoulli("SBM_sym",
                                             mlvl$adjacency_matrix$I,
                                             verbosity = 0, plotting = "")
      sbm_lower$estimate()
      sbm_upper <- blockmodels::BM_bernoulli("SBM_sym",
                                             mlvl$adjacency_matrix$O,
                                             verbosity = 0, plotting = "")
      sbm_upper$estimate()

      fit <- MLVSBM::mlvsbm_estimate_network(mlvl, nb_cores = 2L)
      return(tibble("SBM_ARI_I" =
                      aricode::ARI(
                        c1 = mlvl$memberships$I,
                        c2 = apply(sbm_lower$memberships[[which.max(sbm_lower$ICL)]]$Z,
                   1, which.max)),
                    "SBM_ARI_O" =
                     aricode::ARI(
                       c1 = mlvl$memberships$O,
                       c2 = apply(sbm_upper$memberships[[which.max(sbm_upper$ICL)]]$Z,
                                  1, which.max)),
                    "MLVL_ARI_I" = aricode::ARI(
                      c1 = mlvl$memberships$I,
                      c2 = fit$Z$I),
                    "MLVL_ARI_O" = aricode::ARI(
                      c1 = mlvl$memberships$O,
                      c2 = fit$Z$O)))
    },
    mc.cores = 6)
  res <- bind_rows(res)
  res$d <- as.factor(g)
  res_detect_thresh <- bind_rows(res_detect_thresh, res)
}

res_detect_thresh %>%
  pivot_longer(-d, names_to = "Model", values_to = "ARI") %>%
  ggplot(aes(x = d, y = ARI)) +
  geom_boxplot() +
  facet_wrap(~ Model)


saveRDS("~/Documents/r_project/mlvsbm_analysis/")


