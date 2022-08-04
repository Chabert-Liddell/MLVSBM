#' R6Class for multilevel object
#'
#' @description Store all simulation parameters and list of fittedmodels.
#' Methods for global inference and model selection are included.
#' @import R6
#'
#' @export
MLVSBM <-
  R6::R6Class(
    "MLVSBM",
    ## fields for internal use (referring to mathematical notations)
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
      M            = NULL,  # list of NA masks for CV and missing data for X
      distribution_ = NULL,
      fitted_sbm    = NULL,
      ICLtab_sbm    = NULL
      ),
    public = list(
      #' @param n A list of size 2, the number of nodes
      #' @param X A list of 2 adjacency matrices
      #' @param A The affiliation matrix
      #' @param Z A list of 2 vectors, the blocks membership
      #' @param directed A list of 2 booleans
      #' @param sim_param A list of MLVSBM parameters for simulating networks
      #' @param distribution The distributions of the interactions ("bernoulli")
      #'
      #' @return A MLVSBM object
      #' @description
      #' Constructor for R6 class MLVSBM
      initialize = function(n = NULL, X = NULL, A = NULL,
                            Z = NULL, directed = NULL, sim_param = NULL,
                            distribution = list("bernoulli", "bernoulli")) {
        private$n            = n
        if(! is.null(X)) {
          if(is.null(rownames(X[[1]])) & is.null(colnames(X[[1]]))) {
            rownames(X[[1]]) <- colnames(X[[1]]) <- paste0("I", seq(nrow(X$I)))
          }
          if(is.null(rownames(X[[2]])) & is.null(colnames(X[[2]]))) {
            rownames(X[[2]]) <- colnames(X[[2]]) <- paste0("O", seq(nrow(X$O)))
          }
          rownames(A) <- rownames(X$I)
          colnames(A) <- colnames(X$O)
        }
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
      },
      #' @description
      #' Find a fitted model of a given size
      #' @param nb_clusters A list of the size of the model
      #' @param fit if fit = "best" return the best model according to the ICL
      #' @return A FitMLVSBM object
      findmodel =
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
        },
      #' @description delete all fitted models
      clearmodels =
        function () {
          private$fitted     <-  list()
          private$tmp_fitted <-  list()
          private$ICLtab     <-  NULL
        },
      #' @description Added a FitMLVSBM object to the list of fitted model
      #' @param fit The FitMLVSBM object to be added
      addmodel =
        function (fit) {
          private$fitted <- c(private$fitted, list(fit))
          if (is.null(private$ICLtab)) {
            private$ICLtab <-
              data.frame(index  = as.integer(1),
                         Q_I      = fit$nb_clusters$I,
                         Q_O      = fit$nb_clusters$O,
                         ICL       = fit$ICL)
          } else {
            private$ICLtab <-
              rbind(
                private$ICLtab,
                data.frame(index  = as.integer(nrow(private$ICLtab) +1),
                           Q_I    = fit$nb_clusters$I,
                           Q_O    = fit$nb_clusters$O,
                           ICL    = fit$ICL
                )
              )
          }
        }
      ),
    active = list(
      ## active binding to access fields outside the class
      #' @field nb_nodes List of the umber of nodes for each levels
      nb_nodes              = function(value) private$n,
      #' @field simulation_parameters List of parameters of the MLVSBM
      simulation_parameters = function(value) private$sim_param,
      #' @field affiliation_matrix Access the affiliation matrix
      affiliation_matrix    = function(value)
        if (missing(value)) private$A else private$A = value,
      #' @field adjacency_matrix Access the list of adjacency_matrix
      adjacency_matrix      = function(value)
        if (missing(value)) private$X else private$X = value,
      #' @field memberships Access the list of the clusterings
      memberships           = function(value)
        if (missing(value)) private$Z else private$Z = value,
      #' @field fittedmodels Get the list of selected fitted FitMLVSBM objects
      fittedmodels          = function(value) private$fitted,
      #' @field ICL A summary table of selected fitted models and ICL model selection criterion
      ICL                   = function(value) private$ICLtab,
      #' @field ICL_sbm Summary table of ICL by levels
      ICL_sbm               = function(value) private$ICLtab_sbm,
      #' @field tmp_fittedmodels A list of all fitted FitMLVSBM objects
      tmp_fittedmodels      = function(value) private$tmp_fitted,
      #' @field fittedmodels_sbm A list of selected fitted FitSBM objects of each levels
      fittedmodels_sbm      = function(value) private$fitted_sbm,
      #' @field max_clusters Access the list of maximum model size
      max_clusters          = function(value)
        if (missing(value))  private$max_Q else private$max_Q = value,
      #' @field min_clusters Access the list of minimum model size
      min_clusters          = function(value)
        if (missing(value))  private$min_Q else private$min_Q = value,
      #' @field directed Access the list of boolean for levels  direction
      directed              = function(value)
        if (missing(value)) private$directed_ else private$directed_ = value,
      #' @field directed Access the list of the distribution used for each levels
      distribution          = function(value)
        if (missing(value)) private$distribution_ else private$distribution_ = value
      )
    )

