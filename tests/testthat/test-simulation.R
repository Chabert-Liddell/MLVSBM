test_that("simulate_adjacency", {
  expect_identical(
    rowSums(mlvsbm_simulate_network(
        n = list(60, 40),
        Q = list(2, 2),
        pi = c(.5, .5),
        gamma = matrix(c(.8, .2, .1, .9), ncol = 2),
        alpha = list(matrix(c(.5, .1, .1, .5), nrow = 2, ncol = 2),
                     matrix(c(.7, .4, .4, .1), nrow = 2, ncol = 2)),
        directed = list(FALSE, FALSE),
        affiliation = "preferential",
        no_empty_org = FALSE)$affiliation_matrix)
      ,
    rep(1, 60))
})
