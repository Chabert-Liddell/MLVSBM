#' An R6Class for multilevel object
#' @import R6
#'
#' @export


MLVSBM <-
  R6::R6Class(
    "MLVSBM",
    ## fields for internal use (refering to mathematical notations)
    private = list(
      n            = NULL, # Number of nodes
      sim_param    = NULL, # parameters used for simulation
      X            = NULL, # List of adjacency matrices
      A            = NULL, # Affiliation matrix
      Z            = NULL, # List of latent variables vectors
      fitted       = NULL, # List of fitted model for this network
      ICLtab       = NULL, #
      tmp_fitted   = NULL, # List of all fitted model
      min_Q        = NULL, # List of minimum clusters for inference
      max_Q        = NULL, # List of maximum clusters for inference
      directed_     = NULL, # Are levels directed
      M            = NULL,  # list of NA masks for CV and missig data for X
      distribution_ = NULL,
      fitted_sbm    = NULL,
      ICLtab_sbm    = NULL
      ),
    public = list(
      ## constructor
      initialize = function(n = NULL, X = NULL, A = NULL,
                            Z = NULL, directed = NULL, sim_param = NULL,
                            distribution = list("bernoulli", "bernoulli")) {
        private$n            = n
        private$X            = X
        private$Z            = Z
        private$A            = A
        private$directed_     = directed
        private$sim_param    = sim_param
        private$distribution_ = distribution
        private$min_Q        = list(I = 1,
                                    O = 1)
        private$max_Q        = list(I = floor(sqrt(n[[1]])),
                                    O = floor(sqrt(n[[2]])))
        private$fitted               = list()
        private$tmp_fitted           = list()
        private$fitted_sbm           = list("lower" = list(),
                                            "upper" = list())
        private$ICLtab_sbm   = list("lower" = list(),
                                    "upper" = list())
      }
      ),
    active = list(
      ## active binding to access fields outside the class
      nb_nodes              = function(value) private$n,
      simulation_parameters = function(value) private$sim_param,
      affiliation_matrix    = function(value)
        if (missing(value)) private$A else private$A = value,
      adjacency_matrix      = function(value)
        if (missing(value)) private$X else private$X = value,
      memberships           = function(value)
        if (missing(value)) private$Z else private$Z = value,
      fittedmodels          = function(value) private$fitted,
      ICL                   = function(value) private$ICLtab,
      ICL_sbm               = function(value) private$ICLtab_sbm,
      tmp_fittedmodels      = function(value) private$tmp_fitted,
      fittedmodels_sbm      = function(value) private$fitted_sbm,
      max_clusters          = function(value)
        if (missing(value))  private$max_Q else private$max_Q = value,
      min_clusters          = function(value)
        if (missing(value))  private$min_Q else private$min_Q = value,
      directed              = function(value)
        if (missing(value)) private$directed_ else private$directed_ = value,
      distribution          = function(value)
        if (missing(value)) private$distribution_ else private$distribution_ = value
      )
    )


MLVSBM$set(
  "public", "findmodel",
  function (nb_clusters = NA, fit = NA) {
    if (fit == "best") {
      return(self$ICL %>%
               filter(ICL == max(ICL)) %>%
               filter(index == max(index)) %>%
               select(index) %>%
               as.numeric() %>%
               self$fittedmodels[[.]])
      } else {
        return(self$ICL %>%
                 filter(Q_I == nb_clusters$I, Q_O == nb_clusters$O) %>%
                 filter(ICL == max(ICL)) %>%
                 filter(index == max(index)) %>%
                 select(index) %>%
                 as.numeric() %>%
                 self$fittedmodels[[.]])}
            }
            )
MLVSBM$set(
  "public", "clearmodels",
  function () {
    private$fitted     <-  list()
    private$tmp_fitted <-  list()
    private$ICLtab     <-  NULL
    }
  )
MLVSBM$set(
  "public", "addmodel",
  function (fit) {
    private$fitted = c(private$fitted, list(fit))
    if (is.null(private$ICLtab)) {
      private$ICLtab <-
        dplyr::tibble(index  = as.integer(1),
                      Q_I      = fit$nb_clusters$I,
                      Q_O      = fit$nb_clusters$O,
                      ICL       = fit$ICL
                  )
      } else {
        private$ICLtab <-
          dplyr::bind_rows(
            private$ICLtab,
            dplyr::tibble(index  = as.integer(nrow(private$ICLtab) +1),
                          Q_I    = fit$nb_clusters$I,
                          Q_O    = fit$nb_clusters$O,
                          ICL    = fit$ICL
                          )
          )
        }
    }
  )
