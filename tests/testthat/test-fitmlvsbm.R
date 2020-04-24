set.seed(1)
my_mlvsbm <- mlvsbm_simulate_network(
  n = list(I = 60, O = 40),
  Q = list(I = 2, O = 2),
  pi = c(.2, .8),
  gamma = matrix(c(.8, .2, .1, .9), ncol = 2),
  alpha = list(I = matrix(c(.5, .1, .1, .5), nrow = 2, ncol = 2),
               O = matrix(c(.9, .4, .4, .1), nrow = 2, ncol = 2)),
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
fit <- mlvsbm_estimate_network(my_mlvsbm, nb_cores = 1L)
test_that(
  "is_independent",
  {
    expect_equal(
      TRUE,
      max(my_mlvsbm$ICL_sbm$lower) + max(my_mlvsbm$ICL_sbm$upper) <= max(my_mlvsbm$ICL$ICL)
    )
  })
test_that(
  "good_model_size",
  {
    expect_equal(
      (mlvsbm_estimate_network(my_mlvsbm,
                               nb_clusters = list(I = 2, O = 2),
                               nb_cores = 1L))$nb_clusters,
      list(I = 2., O = 2.)
    )
  }
)
test_that(
  "fit_object",
  {
    expect_equal(
      TRUE,
      inherits(x = mlvsbm_estimate_network(my_mlvsbm, nb_cores = 1L),
               what = "FitMLVSBM")
    )
  }
)
# set.seed(1)
# my_mlvsbm <- MLVSBM::mlvsbm_simulate_network(
#   n = list(100, 50),
#   Q = list(2, 2),
#   pi = c(.5, .5),
#   gamma = matrix(c(.8, .2, .1, .9), ncol = 2),
#   alpha = list(matrix(c(.3, .1, .1, .3), nrow = 2, ncol = 2),
#                matrix(c(.7, .4, .4, .1), nrow = 2, ncol = 2)),
#   directed = list(FALSE, FALSE),
#   affiliation = "preferential",
#   no_empty_org = FALSE)
# fit <- MLVSBM:::mlvsbm_estimate_network(my_mlvsbm)

