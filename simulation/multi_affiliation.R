set.seed(1)
my_mlvsbm <- MLVSBM::mlvsbm_simulate_network(
  n = list(100, 50),
  Q = list(2, 2),
  pi = c(.5, .5),
  gamma = matrix(c(.8, .2, .1, .9), ncol = 2),
  alpha = list(matrix(c(.25, .1, .1, .25), nrow = 2, ncol = 2),
               matrix(c(.7, .4, .4, .1), nrow = 2, ncol = 2)),
  directed = list(FALSE, FALSE),
  affiliation = "preferential",
  no_empty_org = FALSE)


aff <- matrix(0, 100, 50)
nb_aff <- rpois(n = 100, lambda = .5) + 1
for(i in seq(100)) {
  aff[i, sample(x = seq(50), size = nb_aff[i], replace = FALSE)] <- 1/nb_aff[i]
}
my_mlvsbm$affiliation_matrix <- aff
fit <- MLVSBM:::mlvsbm_estimate_network(my_mlvsbm)

fit$Z
