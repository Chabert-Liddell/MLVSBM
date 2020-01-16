# Simulation pour la prédiction de liens

library(tidyverse)
library(MLVSBM)
library(blockmodels)
library(aricode)
setwd(dir = "Documents/Multilevel/Simulations/article/")
#-------------------------------------------------------------------------------
# Function for simulating data
#-------------------------------------------------------------------------------
simulation <- function(n_iter = 30, n = NA, Q = list(I = 3, O = 3),
                       alpha = NA, gamma = NA, pi = NA,
                       directed = list(I = FALSE, O = FALSE),
                       affiliation = "preferential") {
  result <- pbmcapply::pbmclapply(
    X = seq(n_iter),
    FUN = function(i) {
      mlvl <-
        MLVSBM::mlvsbm_simulate_network(
          n = n, Q = Q, pi = pi, gamma = gamma, alpha = alpha,
          directed = directed, affiliation = affiliation)
      sbm_upper <- blockmodels::BM_bernoulli("SBM_sym",
                                             mlvl$adjacency_matrix$O,
                                             verbosity = 0, plotting = "")
      sbm_upper$estimate()
      sbm_lower <- blockmodels::BM_bernoulli("SBM_sym",
                                             mlvl$adjacency_matrix$I,
                                             verbosity = 0, plotting = "")
      sbm_lower$estimate()
      ZR <- apply(sbm_lower$memberships[[which.max(sbm_lower$ICL)]]$Z,
                  1, which.max)
      ZL <- apply(sbm_upper$memberships[[which.max(sbm_upper$ICL)]]$Z,
                  1, which.max)
#      fit <- MLVSBM::mlvsbm_estimate_network(mlv = mlvl)
      fit_bm_init <-
        mlvl$estimate_all_bm(Q = list(I = which.max(sbm_lower$ICL),
                                           O = which.max(sbm_upper$ICL)),
                                  Z = list(I = ZR, O = ZL))
      # arir_mlvl <- aricode::ARI(
      #   c1 = mlvl$memberships$I,
      #   c2 = fit$Z$I)
      # aril_mlvl <- aricode::ARI(
      #   c1 = mlvl$memberships$O,
      #   c2 = fit$Z$O)
      arir_mlvl_bm_init <- aricode::ARI(
        c1 = mlvl$memberships$I,
        c2 = fit_bm_init$Z$I)
      aril_mlvl_bm_init <- aricode::ARI(
        c1 = mlvl$memberships$O,
        c2 = fit_bm_init$Z$O)
      arir_sbm <- aricode::ARI(
        c1 = mlvl$memberships$I,
        c2 = ZR)
      aril_sbm <- aricode::ARI(
        c1 = mlvl$memberships$O,
        c2 = ZL)
      nb_rclust_sbm <- which.max(sbm_lower$ICL)
      nb_lclust_sbm <- which.max(sbm_upper$ICL)
      # nb_rclust_mlvl <- fit$nb_clusters$I
      # nb_lclust_mlvl <- fit$nb_clusters$O
      nb_rclust_mlvl_bm_init <- fit_bm_init$nb_clusters$I
      nb_lclust_mlvl_bm_init <- fit_bm_init$nb_clusters$O
      return(tibble("SBM_ARI_R" = arir_sbm,
                    "SBM_ARI_L" = aril_sbm,
                    # "MLVL_ARI_R" = arir_mlvl,
                    # "MLVL_ARI_L" = aril_mlvl,
                    "MLVL_INIT_ARI_R" = arir_mlvl_bm_init,
                    "MLVL_INIT_ARI_L" = aril_mlvl_bm_init,
                    "SBM_R_CLUST" = nb_rclust_sbm,
                    "SBM_L_CLUST" = nb_lclust_sbm,
                    # "MLVL_R_CLUST" = nb_rclust_mlvl,
                    # "MLVL_L_CLUST" = nb_lclust_mlvl,
                    "MLVL_INIT_R_CLUST" = nb_rclust_mlvl_bm_init,
                    "MLVL_INIT_L_CLUST" = nb_lclust_mlvl_bm_init)
             )
    })
  return(result)
}



#-------------------------------------------------------------------------------
# Set level 1
# Community structure for macro level (planted partition) (varying)
# Preferential attachment interlevel
# Community structure for micro level (varying)
# Gamma is near diagonal
#-------------------------------------------------------------------------------
description_1 <- "#-------------------------------------------------------------------------------
# Set level 1
# Community structure for macro level (planted partition) (varying)
# Preferential attachment interlevel
# Community structure for micro level (varying)
# Gamma is near diagonal
#-------------------------------------------------------------------------------"

n <- list(I = 60*3, O = 20*3)
Q <- list(I = 3, O = 3)
delta <- c(.01, .05, .1)
alpha <- list(I = NULL, O = NULL)
alpha$O <- .1*(diag(x = 4, nrow = 3, ncol = 3) + 1)
pi    <- rep(1/3, 3)
gamma   <- .1*(diag(7, 3, 3) + 1)

epsilon <- seq(0, 5, 1)
sim_1 <- tibble()

set.seed(1)

for (d in seq(delta)) {
  result_tmp <- tibble()
  for(i in seq(epsilon)) {
    alpha$I <- delta[d]*(diag(x = epsilon[i], nrow = 3, ncol = 3) + 1)
    result <- simulation(n = n, Q = Q,
                         alpha = alpha,
                         gamma = gamma, pi = pi, n_iter = 30)
    result <- unnest(tibble(result))
    result$eps <- epsilon[i]
    result_tmp <- bind_rows(result_tmp, result)
  }
  result_tmp$d <- as.factor(delta[d])
  sim_1 <- bind_rows(sim_1, result_tmp)
}

#save(list = c("sim_1", "description_1"),
 #    file = "./sim_1.Rdata")

ggplot(gather(sim_1, "model", "ARI", "MLVL_INIT_ARI_R", "SBM_ARI_R")) +
  stat_summary(mapping = aes(x = eps, y = ARI, col = model, linetype = as.factor(d)),
               geom = "line", fun.y = mean)

#
#
# sim_1 %>% group_by(eps) %>% group_by(d, add = TRUE) %>%
#   summarise_at(c("SBM_R_CLUST", "MLVL_R_CLUST"), mean)


#-------------------------------------------------------------------------------
# Set level 2
# Community structure for macro level (planted partition) (varying)
# Preferential attachment interlevel
# Community structure for micro level (planted partition)
# Gamma is near diagonal
# Same as sim 1 but we invert alpha R and L
#-------------------------------------------------------------------------------
description_2 <- "#-------------------------------------------------------------------------------
# Set level 2
# Community structure for macro level (planted partition) (varying)
# Preferential attachment interlevel
# Community structure for micro level (planted partition)
# Gamma is near diagonal
# Same as sim 1 but we invert alpha R and L
#-------------------------------------------------------------------------------"
n_L     <- 20*3
n_R     <- 60*3

delta <- c(.01, .05, .1)
alpha_R <- .1*(diag(x = 4, nrow = 3, ncol = 3) + 1)
pi_L    <- rep(1/3, 3)
gamma   <- .1*(diag(7, 3, 3) + 1)

epsilon <- seq(0, 8, by = .5)
sim_2 <- tibble()

set.seed(1)

for (d in seq(delta)) {
  result_tmp <- tibble()
  for(i in seq(epsilon)) {
    alpha_L <- delta[d]*(diag(x = epsilon[i], nrow = 3, ncol = 3) + 1)
    result <- simulation(n_L = n_L, n_R = n_R, Q_R = 3, Q_L = 3,
                         alpha_R = alpha_R, alpha_L = alpha_L,
                         gamma = gamma, pi_L = pi_L, n_iter = 50)
    result <- unnest(tibble(result))
    result$eps <- epsilon[i]
    result_tmp <- bind_rows(result_tmp, result)
  }
  result_tmp$d <- delta[d]
  sim_2 <- bind_rows(sim_2, result_tmp)
}
save(list = c("sim_2", "description_2"),
     file = "./sim_2.Rdata")
# ggplot(gather(sim_2, "model", "ARI_L", "SBM_ARI_L_CLUST", "MLVL_ARI_L_CLUST")) +
#   stat_summary(mapping = aes(x = eps, y = ARI_L, linetype = model, col = as.factor(d)),
#                geom = "line", fun.y = mean)
# sim_2 %>% group_by(eps) %>% group_by(d, add = TRUE) %>%
#   summarise_at(c("SBM_L_CLUST", "MLVL_L_CLUST"), mean) %>% print()

#-------------------------------------------------------------------------------
# Set level 3
# Community structure for macro level (planted partition)
# Preferential attachment interlevel
# dissortative structure for micro level (varying)
# Gamma is near diagonal
#-------------------------------------------------------------------------------
description_3 <- "#-------------------------------------------------------------------------------
# Set level 3
# Community structure for macro level (planted partition)
# Preferential attachment interlevel
# dissortative structure for micro level (varying)
# Gamma is near diagonal
#-------------------------------------------------------------------------------"
n_L     <- 20*3
n_R     <- 60*3
delta <- c(.01, .05, .1)
alpha_L <- .1*(diag(x = 4, nrow = 3, ncol = 3) + 1)
pi_L    <- rep(1/3, 3)
gamma   <- .1*(diag(7, 3, 3) + 1)

epsilon <- seq(0, 8, by = .5)
sim_3 <- tibble()

set.seed(1)

for(d in seq(delta)) {
  result_tmp <- tibble()
  for(i in seq(epsilon)) {
    print(i)
    alpha_R <- delta[d]*(diag(x = - epsilon[i], nrow = 3, ncol = 3) + 1 + epsilon[i])
    result <- simulation(n_L = n_L, n_R = n_R, Q_R = 3, Q_L = 3,
                         alpha_R = alpha_R, alpha_L = alpha_L,
                         gamma = gamma, pi_L = pi_L, n_iter = 50)
    result <- unnest(tibble(result))
    result$eps <- epsilon[i]
    result_tmp <- bind_rows(result_tmp, result)
  }
  result_tmp$d <- delta[d]
  sim_3 <- bind_rows(sim_3, result_tmp)
}
save(list = c("sim_3", "description_3"),
     file = "./sim_3.Rdata")
ggplot(gather(sim_3, "model", "ARI_R", "SBM_ARI_R_CLUST", "MLVL_ARI_R_CLUST")) +
  stat_summary(mapping = aes(x = eps, y = ARI_R, linetype = model, col = as.factor(d)),
               geom = "line", fun.y = mean)
sim_3 %>% group_by(eps) %>% group_by(d, add = TRUE) %>%
  summarise_at(c("SBM_R_CLUST", "MLVL_R_CLUST"), mean) %>% print()
ggplot(sim_3 %>% group_by(eps) %>% group_by(d, add = TRUE) %>%
         summarise_at(c("SBM_R_CLUST", "MLVL_R_CLUST"), mean) %>%
         gather("model", "CLUST", "SBM_R_CLUST", "MLVL_R_CLUST")) +
  facet_wrap(cols = vars(model)) +
  geom_raster(mapping = aes(x = as.integer(eps), y = as.factor(d), fill = CLUST))

#-------------------------------------------------------------------------------
# Set level 4
# Community structure for macro level (planted partition)
# Preferential attachment interlevel
# dissortative structure for micro level (varying)
# Gamma is near diagonal
#-------------------------------------------------------------------------------
description_4 <- "#-------------------------------------------------------------------------------
# Set level 4
# Community structure for macro level (planted partition)
# Preferential attachment interlevel
# dissortative structure for micro level (varying)
# Gamma is near diagonal
# Same as Set 3 with exchanged alpha
#-------------------------------------------------------------------------------"
n_L     <- 20*3
n_R     <- 60*3
delta <- c(.01, .05, .1)
alpha_R <- .1*(diag(x = 4, nrow = 3, ncol = 3) + 1)
pi_L    <- rep(1/3, 3)
gamma   <- .1*(diag(7, 3, 3) + 1)

epsilon <- seq(0, 8, by = .5)
sim_4 <- tibble()

set.seed(1)

for(d in seq(delta)) {
  result_tmp <- tibble()
  for(i in seq(epsilon)) {
    print(i)
    alpha_L <- delta[d]*(diag(x = - epsilon[i], nrow = 3, ncol = 3) + 1 + epsilon[i])
    result <- simulation(n_L = n_L, n_R = n_R, Q_R = 3, Q_L = 3,
                         alpha_R = alpha_R, alpha_L = alpha_L,
                         gamma = gamma, pi_L = pi_L, n_iter = 50)
    result <- unnest(tibble(result))
    result$eps <- epsilon[i]
    result_tmp <- bind_rows(result_tmp, result)
  }
  result_tmp$d <- delta[d]
  sim_4 <- bind_rows(sim_4, result_tmp)
}
save(list = c("sim_4", "description_4"),
     file = "./sim_4.Rdata")

#-------------------------------------------------------------------------------
# Set level 5
# Community structure for macro level (planted partition)
# Preferential attachment interlevel
# Core periphery for micro level (varying)
# Gamma is near diagonal
#-------------------------------------------------------------------------------
description_5 <- "#-------------------------------------------------------------------------------
# Set level 5
# Community structure for macro level (planted partition)
# Preferential attachment interlevel
# Core periphery for micro level (varying)
# Gamma is near diagonal
#-------------------------------------------------------------------------------"
n_L     <- 20*3
n_R     <- 60*3
delta <- c(.01, .05, .1)
alpha_L <- .1*(diag(x = 4, nrow = 3, ncol = 3) + 1)
pi_L    <- rep(1/3, 3)
gamma   <- .1*(diag(7, 3, 3) + 1)

epsilon <- seq(0, 8, by = .5)
sim_5 <- tibble()

set.seed(1)

for(d in seq(delta)) {
  result_tmp <- tibble()
  for(i in seq(epsilon)) {
    print(i)
    alpha_R <- delta[d]*matrix(c(1+epsilon[i],   1+epsilon[i], 1,
                                 1+epsilon[i], 1, 1,
                                 1, 1, 1), ncol = 3, nrow = 3)
    result <- simulation(n_L = n_L, n_R = n_R, Q_R = 3, Q_L = 3,
                         alpha_R = alpha_R, alpha_L = alpha_L,
                         gamma = gamma, pi_L = pi_L, n_iter = 50)
    result <- unnest(tibble(result))
    result$eps <- epsilon[i]
    result_tmp <- bind_rows(result_tmp, result)
  }
  result_tmp$d <- delta[d]
  sim_5 <- bind_rows(sim_5, result_tmp)
}
save(list = c("sim_5", "description_5"),
     file = "./sim_5.Rdata")

#-------------------------------------------------------------------------------
# Set level 6
# Community structure for micro level (planted partition)
# Preferential attachment interlevel
# Core periphery for macro level (varying)
# Gamma is near diagonal
# Same as Set 5 with inversed alpha
#-------------------------------------------------------------------------------
description_6 <- "#-------------------------------------------------------------------------------
# Set level 6
# Community structure for micro level (planted partition)
# Preferential attachment interlevel
# Core periphery for macro level (varying)
# Gamma is near diagonal
# Same as set 5 with inversed alpha
#-------------------------------------------------------------------------------"
n_L     <- 20*3
n_R     <- 60*3
delta <- c(.01, .05, .1)
alpha_R <- .1*(diag(x = 4, nrow = 3, ncol = 3) + 1)
pi_L    <- rep(1/3, 3)
gamma   <- .1*(diag(7, 3, 3) + 1)

epsilon <- seq(0, 8, by = .5)
sim_6 <- tibble()

set.seed(1)

for(d in seq(delta)) {
  result_tmp <- tibble()
  for(i in seq(epsilon)) {
    print(i)
    alpha_L <- delta[d]*matrix(c(1+epsilon[i],   1+epsilon[i], 1,
                                 1+epsilon[i], 1, 1,
                                 1, 1, 1), ncol = 3, nrow = 3)
    result <- simulation(n_L = n_L, n_R = n_R, Q_R = 3, Q_L = 3,
                         alpha_R = alpha_R, alpha_L = alpha_L,
                         gamma = gamma, pi_L = pi_L, n_iter = 50)
    result <- unnest(tibble(result))
    result$eps <- epsilon[i]
    result_tmp <- bind_rows(result_tmp, result)
  }
  result_tmp$d <- delta[d]
  sim_6 <- bind_rows(sim_6, result_tmp)
}
save(list = c("sim_6", "description_6"),
     file = "./sim_6.Rdata")

