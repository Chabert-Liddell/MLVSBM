#' @import R6
#' @export
MLVSBM <-
  R6::R6Class(
    "MLVSBM",
    ## fields for internal use (refering to mathematical notations)
    private = list(
      n          = NULL, # Number of node
      sim_param  = NULL, # parameters used for simulation
      X          = NULL, # List of adjacency matrices
      A          = NULL, # Affiliation matrix
      Z          = NULL, # List of latent variables vectors
      fitted     = NULL, # List of fitted model for this network
      ICLtab     = NULL, #
      tmp_fitted = NULL, # List of all fitted model
      min_Q      = NULL, # List of minimum clusters for inference
      max_Q      = NULL, # List of maximum clusters for inference
      directed   = NULL, # Are levels directed
      M          = NULL  # list of NA masks for CV and missig data for X
      ),
    public = list(
      ## constructor
      initialize = function(n = NULL, X = NULL, A = NULL,
                            Z = NULL, directed = NULL, sim_param = NULL) {
        n         = n
        X         = X
        Z         = Z
        A         = A
        directed  = directed
        sim_param = sim_param
      }
      ),
    active = list(
      ## active binding to access fields outside the class
      nb_nodes              = function(value) private$n,
      simulation_parameters = function(value) sim_param,
      affiliation_matrox    = function(value) {
        if (missing(value)) private$A else private$A = value},
      adjacency_matrix      = function(value) {
        if (missing(value)) private$X else private$X = value},
      memberships           = function(value) {
        if (missing(value)) private$Z else private$Z = value},
      fittedmodels          = function(value) private$fitted,
      ICL                   = function(value) private$ICLtab,
      tmp_fittedmodels      = function(value) private$tmp_fitted,
      max_clusters          = function(value)
        if (missing(value))  private$max_Q else private$max_Q = value,
      min_clusters          = function(value)
        if (missing(value))  private$min_Q else private$min_Q = value
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
