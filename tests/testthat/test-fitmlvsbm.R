set.seed(1)
my_mlvsbm <- mlvsbm_simulate_network(
  n = list(60, 40),
  Q = list(2, 2),
  pi = c(.5, .5),
  gamma = matrix(c(.8, .2, .1, .9), ncol = 2),
  alpha = list(matrix(c(.5, .1, .1, .5), nrow = 2, ncol = 2),
               matrix(c(.7, .4, .4, .1), nrow = 2, ncol = 2)),
  directed = list(FALSE, FALSE),
  affiliation = "preferential",
  no_empty_org = FALSE)

fit_hier <- FitMLVSBM$new(Q = my_mlvsbm$simulation_parameters$Q,
                            A = my_mlvsbm$affiliation_matrix,
                            X = my_mlvsbm$adjacency_matrix,
                            directed = my_mlvsbm$directed,
                            distribution = my_mlvsbm$distribution,
                            independent = FALSE)
fit_hier$do_vem(init = "hierarchical")
fit_spec <- FitMLVSBM$new(Q = my_mlvsbm$simulation_parameters$Q,
                            A = my_mlvsbm$affiliation_matrix,
                            X = my_mlvsbm$adjacency_matrix,
                            directed = my_mlvsbm$directed,
                            distribution = my_mlvsbm$distribution,
                            independent = FALSE)
fit_spec$do_vem(init = "spectral")
test_that("spec_hierarchical", {
  expect_equal(
    c(MLVSBM:::ARI(fit_hier$Z$I, fit_spec$Z$I),
      MLVSBM:::ARI(fit_hier$Z$O, fit_spec$Z$O)),
    c(1, 1))
})
set.seed(1)
my_mlvsbm <- MLVSBM::mlvsbm_simulate_network(
  n = list(100, 50),
  Q = list(2, 2),
  pi = c(.5, .5),
  gamma = matrix(c(.8, .2, .1, .9), ncol = 2),
  alpha = list(matrix(c(.15, .1, .1, .15), nrow = 2, ncol = 2),
               matrix(c(.7, .4, .4, .1), nrow = 2, ncol = 2)),
  directed = list(FALSE, FALSE),
  affiliation = "preferential",
  no_empty_org = FALSE)
fit <- MLVSBM:::mlvsbm_estimate_network(my_mlvsbm)
test_that(
  "is_independent",
  {
    expect_equal(
      TRUE,
      max(my_mlvsbm$ICL_sbm$lower) + max(my_mlvsbm$ICL_sbm$upper) <= my_mlvsbm$ICL
    )
  })

