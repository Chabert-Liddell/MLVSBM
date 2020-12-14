test_that("simulate_adjacency", {
  expect_equal(
    unname(rowSums(mlvsbm_simulate_network(
        n = list(60, 40),
        Q = list(2, 2),
        pi = c(.5, .5),
        gamma = matrix(c(.8, .2, .1, .9), ncol = 2),
        alpha = list(matrix(c(.5, .1, .1, .5), nrow = 2, ncol = 2),
                     matrix(c(.7, .4, .4, .1), nrow = 2, ncol = 2)),
        directed = list(FALSE, FALSE),
        affiliation = "preferential",
        no_empty_org = FALSE)$affiliation_matrix))
      ,
    rep(1, 60))
})
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
MR <- my_mlvsbm$adjacency_matrix$I
ML <- my_mlvsbm$adjacency_matrix$O
MR[sample(seq(my_mlvsbm$nb_nodes$I), size = 12, replace = FALSE),
   sample(seq(my_mlvsbm$nb_nodes$I), size = 12, replace = FALSE)] <- NA
ML[sample(seq(my_mlvsbm$nb_nodes$O), size = 8, replace = FALSE),
   sample(seq(my_mlvsbm$nb_nodes$O), size = 8, replace = FALSE)] <- NA
my_na <- mlvsbm_create_network(X = list(I = MR, O = ML),
                               A = my_mlvsbm$affiliation_matrix,
                               directed = list(FALSE, FALSE))
test_that("dealing_na", {
  expect_equal(MR, my_na$adjacency_matrix$I)
  expect_equal(ML, my_na$adjacency_matrix$O)
})
