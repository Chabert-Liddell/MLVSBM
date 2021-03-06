% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/user_function.R
\name{mlvsbm_log_likelihood}
\alias{mlvsbm_log_likelihood}
\title{Compute the complete log likelihood of a multilevel network for a given
clustering of the nodes.}
\usage{
mlvsbm_log_likelihood(mlv, clustering)
}
\arguments{
\item{mlv}{A MLVSBM object, the network data}

\item{clustering}{A list of 2 vectors of integers of the same length as
the number of node of each level.}
}
\value{
A numeric, the log likelihood of the multilevel network
for the given clustering.
}
\description{
This function is useful to compute the likelihood for clusters
obtained by different methods.
}
\examples{
my_mlvsbm <- MLVSBM::mlvsbm_simulate_network(
  n = list(I = 40, O = 20), # Number of nodes for the lower level and the upper level
  Q = list(I = 2, O = 2), # Number of blocks for the lower level and the upper level
  pi = c(.3, .7), # Block proportion for the upper level, must sum to one
  gamma = matrix(c(.9, .2,   # Block proportion for the lower level,
                   .1, .8), # each column must sum to one
                 nrow = 2, ncol = 2, byrow = TRUE),
  alpha = list(I = matrix(c(.8, .2,
                            .2, .1),
                          nrow = 2, ncol = 2, byrow = TRUE), # Connection matrix
               O = matrix(c(.99, .3,
                            .3, .1),
                          nrow = 2, ncol = 2, byrow = TRUE)),# between blocks
  directed = list(I = FALSE, O = FALSE), # Are the upper and lower level directed or not ?
  affiliation = "preferential") # How the affiliation matrix is generated
mlvsbm_log_likelihood(mlv = my_mlvsbm, clustering = my_mlvsbm$memberships)
}
