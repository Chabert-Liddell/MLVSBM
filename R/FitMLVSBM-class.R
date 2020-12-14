#' An R6 Class object, a fitted multilevel network once $dovem() is done
#'
#' @import R6
#' @export
FitMLVSBM <-
  R6::R6Class(
    "FitMLVSBM",
    private = list(
      n            = NULL, # number of nodes in I and O
      Q            = NULL, # Number of clusters in I and O
      param        = NULL, # List of fitted parameters of the model
      tau          = NULL, # List of variational parameters of the fitted model
      A            = NULL, # Affiliation Matrix
      X            = NULL, # List of adjacency Matrices
      directed_    = NULL,  # is XI directed ? Is XO directed ?
      M            = NULL,  # List of mask matrix
      distribution_ = NULL, # List of the distribution of X
      independent_ = NULL
    ),
    ##
    ##  Public Methods
    ##
    public = list(
      #' @description Constructor for the FitMLVSBM class
      #' @param Q List of number of blocks
      #' @param A Affiliation matrix
      #' @param X List of adjacency matrices
      #' @param M List of Mask matrices
      #' @param directed List of boolean
      #' @param distribution List of string
      #' @param independent Boolean
      #'
      #' @return A FitMLVSBM object
      initialize = function(Q = list(I = 1, O = 1),
                            A = NA, X = NA,
                            M = list(I = NA, O = NA),
                            directed = NA,
                            distribution = list("bernoulli", "bernoulli"),
                            independent = FALSE) {
        private$A     = A
        private$X     = X
        private$n     = list(I = nrow(A), O = ncol(A))
        private$Q     = Q
        private$directed_ = directed
        private$distribution_ = distribution
        private$independent_    = independent
        if (is.na(M$I)) {
          private$M$I <- diag(-1, private$n$I)
        } else {
          private$M$I          <- M$I
          private$X$I[M$I == 0] <- -1
          diag(private$X$I)    <- 0
        }
        if (is.na(M$O)) {
          private$M$O <- diag(-1, private$n$O)
        } else {
          private$M$O           <- M$O
          private$X$O[M$O == 0] <- -1
          diag(private$X$O)     <- 0
        }
        private$M$I    = if (is.na(M$I)) diag(-1, private$n$I) + 1 else M$I
        private$M$O    = if (is.na(M$O)) diag(-1, private$n$O) + 1 else M$O
        if (sum(is.na(private$X$I)) + sum(is.na(private$X$O)) > 0) {
          private$M$I[is.na(private$X$I)] <- 0
          private$X$I[is.na(private$X$I)] <- -1
          private$M$O[is.na(private$X$O)] <- 0
          private$X$O[is.na(private$X$O)] <- -1
        }
      },
      #' @field vbound The vector of variational bound for monitoring convergence
      vbound = NULL,
      ##------------------------------------------------------------------------
      ##      ## Parameters update
      ##------------------------------------------------------------------------
      #' @description  Update the connection parameters for the M step
      #' @param safeguard Parameter live in a compact [safeguard, 1-safeguard]
      update_alpha =
        function(safeguard = 2*.Machine$double.eps) {
          ## alpha
          if(private$Q$I == 1) {
            private$param$alpha$I <-
              as.matrix( sum(private$M$I * private$X$I)/ sum(private$M$I))
          } else {
            alpha <-
              crossprod(private$tau$I,
                        (private$M$I*private$X$I) %*% private$tau$I) /
              crossprod(private$tau$I,
                        private$M$I %*% private$tau$I)
            private$param$alpha$I <- alpha
          }
          if(private$Q$O == 1) {
            private$param$alpha$O <-
              as.matrix( sum(private$X$O) / sum(private$M$O))
          } else {
            alpha <-
              crossprod(private$tau$O,
                        (private$M$O * private$X$O) %*% private$tau$O) /
              crossprod(private$tau$O,
                        private$M$O %*% private$tau$O)
            private$param$alpha$O <- alpha
          }
          return (private$param$alpha)
        },
      #' @description Update the upper level mixture parameter for the M step
      #' @param safeguard Parameter live in a compact [safeguard, 1-safeguard]
      update_pi =
        function(safeguard = 1e-3) {
          ## rho
          if (private$Q$O == 1) {
            private$param$pi$O = 1
          } else {
            pi <- colMeans(private$tau$O)
            pi[pi < safeguard] <- safeguard
            private$param$pi$O <- pi/sum(pi)
          }
          return(private$param$pi$O)
        },
      #' @description Update the lower level mixture parameter for the M step
      #' @param safeguard Parameter live in a compact [safeguard, 1-safeguard]
      update_gamma =
        function(safeguard = 1e-6){
          ## gamma
          if (private$independent_ == FALSE) {
            if (private$Q$I == 1) {
              private$param$gamma <-
                matrix(1, nrow = 1, ncol = private$Q$O)
            } else {
              gamma <-
                crossprod(private$tau$I,
                          private$A) %*% private$tau$O
              gamma <-
                t(t(gamma)/colSums(private$A %*% private$tau$O))
              gamma[gamma < safeguard] <- safeguard
              private$param$gamma <- t(t(gamma)/colSums(gamma))
            }
            return(private$param$gamma)
          }
          if (private$independent_ == TRUE) {
            if (private$Q$I == 1) {
              private$param$gamma <-
                matrix(1, nrow = 1, ncol = private$Q$O)
            } else {
              gamma <- matrix(colMeans(private$tau$I),
                              nrow = private$Q$I,
                              ncol = private$Q$O)
              gamma[gamma < safeguard] <- safeguard
              private$param$gamma <- t(t(gamma)/colSums(gamma))
            }
            return(private$param$gamma)
          }
        },

      #-------------------------------------------------------------------------
      #  Inference
      #-------------------------------------------------------------------------
      #' @description   init_clustering Initial clustering for VEM algorithm
      #' @param safeguard Parameter live in a compact [safeguard, 1-safeguard]
      #'
      #' @param method Algorithm used to initiate the clustering, either
      #' "spectral", "hierarchical" or "merge_split" (if \code{Z} is provided)
      #' @param Z Initial clustering if provided
      #'
      init_clustering =
        function(safeguard = 2*.Machine$double.eps,
                 method = "hierarchical",
                 Z = NULL) {
          if (private$Q$I == 1) {
            private$tau$I <-
              matrix(1, nrow = private$n$I, ncol = 1)
          } else {
            init_clust <-
              switch(method,
                     "spectral"     = spcClust(private$X$I, private$Q$I),
                     "hierarchical" = hierarClust(private$X$I, private$Q$I),
                     "merge_split"  = Z$I)
            private$tau$I <-
              1 * sapply(X = seq(private$Q$I), FUN = function(x) init_clust %in% x)
            private$tau$I[private$tau$I < safeguard] <- safeguard
            private$tau$I <-
              private$tau$I / rowSums(private$tau$I)
          }
          if(private$Q$O == 1){
            private$tau$O <-
              matrix(1, nrow = private$n$O, ncol = 1)
          } else {
            init_clust <-
              switch(method,
                     "spectral"     = spcClust(private$X$O, private$Q$O),
                     "hierarchical" = hierarClust(private$X$O, private$Q$O),
                     "merge_split"  = Z$O)
            private$tau$O <-
              1 * sapply(X = seq(private$Q$O), FUN = function(x) init_clust %in% x)
            private$tau$O[private$tau$O < safeguard] <-  safeguard
            private$tau$O <-
              private$tau$O/rowSums(private$tau$O)
          }
        },
      #' @description Reset all parameters
      clear = function(){
        private$param = NULL
        private$tau = NULL
      },
      #-------------------------------------------------------------------------
      # Varational EM algorithm
      #-------------------------------------------------------------------------
      #' @description m_step Compute the M step of the VEM algorithm
      #' @param safeguard Parameter live in a compact [safeguard, 1-safeguard]
      m_step =
        function(safeguard = 1e-6){
          self$update_alpha(safeguard = safeguard)
          self$update_pi(safeguard = safeguard)
          self$update_gamma(safeguard = safeguard)
        },
        #' @description Compute the VE step of the VEM algorithm
        #' @param threshold The convergence threshold
        #' @param fixPointIter The maximum number of fixed point iterations
        #' @param safeguard Parameter live in a compact [safeguard, 1-safeguard]
        ve_step =
          function(threshold = 1e-6, fixPointIter = 10, safeguard = 1e-6){
            condition <- TRUE
            it        <- 0
            tau_old   <- private$tau
            tau       <- private$tau
            while (condition) {
              ## sigma
              tau$I <-
                (private$M$I * private$X$I) %*%
                tcrossprod(tau_old$I, log(private$param$alpha$I)) +
                (private$M$I * (1 - private$X$I)) %*%
                tcrossprod(tau_old$I, log(1 - private$param$alpha$I)) +
                private$A %*% tcrossprod(tau_old$O, log(private$param$gamma))
              if (private$Q$I == 1) {
                tau$I  <-
                  as.matrix(exp( apply(X = tau$I,
                                       MARGIN = 1,
                                       FUN = function(x) x - max(x))),
                            nrow = private$n$I, ncol = 1 )
              } else {
                tau$I <- exp( t(apply(X = tau$I,
                                      MARGIN = 1,
                                      FUN = function(x) x - max(x))) )
              }
              tau$I[tau$I < safeguard] <-  safeguard
              tau$I <- tau$I/rowSums(tau$I)
              ## tau
              tau$O <-
                matrix(log(private$param$pi$O), private$n$O, private$Q$O, byrow = TRUE) +
                (private$M$O * private$X$O) %*%
                tcrossprod(tau_old$O, log(private$param$alpha$O)) +
                (private$M$O * (1 - private$X$O)) %*%
                tcrossprod(tau_old$O,  log(1 - private$param$alpha$O)) +
                crossprod(private$A, tau_old$I) %*% log(private$param$gamma)
              if (private$Q$O == 1) {
                tau$O <- as.matrix(exp(apply(X = tau$O,
                                             MARGIN = 1,
                                             FUN = function(x) x - max(x))),
                                   nrow = private$n$O, ncol = 1 )
              } else {
                tau$O  <- exp(t(apply(X = tau$O,
                                      MARGIN = 1,
                                      FUN = function(x) x - max(x))) )
              }
              tau$O[tau$O < safeguard] = safeguard
              tau$O  <-  tau$O/rowSums(tau$O)
              it  <-  it + 1
              condition  <-  (max(dist_param(tau$O, tau_old$O),
                                  dist_param(tau$I, tau_old$I)) > threshold &&
                                it <= fixPointIter)
              tau_old   <- tau
            }
            private$tau <-  tau
          },
        #' @param init The method for \code{self$init_clustering}
        #' @param threshold The convergence threshold
        #' @param maxIter The max number of VEM iterations
        #' @param fixPointIter The max number of fixed point iterations for VE step
        #' @param safeguard Parameter live in a compact [safeguard, 1-safeguard]
        #' @param Z Initial clustering if provided
        #'
        #' @description Launch a Variational EM algorithm
        do_vem =
          function(init = "hierarchical", threshold = 1e-6,
                   maxIter = 1000, fixPointIter = 100,
                   safeguard = 1e-6, Z = NULL) {
            self$init_clustering(method = init, safeguard = safeguard, Z = Z)
            self$m_step(safeguard = safeguard)
            self$vbound <-  c(self$vbound, self$bound)
            condition   <-  TRUE
            it          <-  0
            if (private$Q$I != 1 | private$Q$O != 1) {
              while (condition) {
                param_old <- private$param
                tau_old   <- private$tau
                bound_old <- self$vbound[length(self$vbound)]
                self$ve_step(safeguard = safeguard)
                self$m_step(safeguard = safeguard)
                if (bound_old > self$bound) {
                  private$tau <- tau_old
                  private$param     <- param_old
                  condition         <- FALSE
                } else {
                  it          <-  it + 1
                  self$vbound <-  c(self$vbound, self$bound)
                  #            cat(it, " : ", self$bound, "\r" )
                  condition <-
                    (max(c(dist_param(private$param$alpha$I, param_old$alpha$I),
                           dist_param(private$param$alpha$O, param_old$alpha$O),
                           dist_param(private$param$gamma, param_old$gamma)))
                     > threshold &
                       it <= maxIter)
                }
              }
              self$permute_empty_class()
            }
          },
        #' @description permute_empty_class Put empty blocks numbers at the end
        permute_empty_class =
          function() {
            if (length(unique(self$Z$I)) < private$Q$I) {
              perm  <-  c(unique(self$Z$I), setdiff(seq(private$Q$I), self$Z$I))
              private$tau$I <-  private$tau$I[, perm]
              private$param$alpha$I <-  private$param$alpha$I[perm, perm]
              private$param$gamma <- matrix(data = private$param$gamma[perm,],
                                            nrow = private$Q$I, ncol = private$Q$O)
            }
            if (length(unique(self$Z$O)) < private$Q$O) {
              perm <-  c(unique(self$Z$O),
                         setdiff( seq(private$Q$O), self$Z$O))
              private$tau$O <-  private$tau$O[, perm]
              private$param$alpha$O   <-  private$param$alpha$O[perm, perm]
              private$param$pi$O      <-  private$param$pi$O[perm]
              private$param$gamma     <-  matrix(data = private$param$gamma[,perm],
                                                 nrow = private$Q$I, ncol = private$Q$O)
            }
          },
      #' @description Plot of FitMLVSBM objects
      #' @param type A string for the type of plot, just "matrix" for now
      #' @return a ggplot2 object
      plot = function(type = c('matrix')) {
        if(type == "matrix") {
          if (! requireNamespace("ggplot2", quietly = TRUE)) {
            stop("Please install ggplot2: install.packages('ggplot2')")
          }
          if (! requireNamespace("ggraph", quietly = TRUE)) {
            stop("Please install ggraph: install.packages('ggraph')")
          }
          if (! requireNamespace("tidygraph", quietly = TRUE)) {
            stop("Please install tidygraph: install.packages('tidygraph')")
          }
          p <- plot_multilevel_matrix(private$X, self$X_hat, private$A, self$Z)
        }
        p
      },
      #' @description print method
      #' @param type character to tune the displayed name
      show = function(type = "Multilevel Stochastic Block Model") {
        cat(type, "--", self$distribution[[1]], "variant\n")
        cat("=====================================================================\n")
        cat("Dimension = (", dim(private$A), ") - (",
            self$nb_clusters$I, self$nb_clusters$O,  ") blocks.\n")
        cat("=====================================================================\n")
        cat("* Useful fields \n")
        cat("  $independent, $distribution, $nb_nodes, $nb_clusters, $Z \n")
        cat("  $membership, $parameters, $ICL, $vbound, $X_hat \n")
      },
      #' @description print method
      print = function() self$show()
    ),
    ##
    ## Active fields
    ##
    active = list(
      ## accessor and
      #' @field affiliation_matrix Get the affiliation matrix
      affiliation_matrix = function(value) private$A,
      #' @field adjacency_matrix Get the list of adjacency matrices
      adjacency_matrix   = function(value) private$X,
      #' @field nb_nodes Get the list of the number of nodes
      nb_nodes           = function(value) private$n,
      #' @field nb_clusters Get the list of the number of blocks
      nb_clusters        = function(value) private$Q,
      #' @field parameters Get the list of the model parameters
      parameters = function(value) {
        if(missing(value)) return(private$param)
        else private$param <- value
        },
      #' @field membership Get the list of the variational parameters
      membership = function(value) {
        if(missing(value)) return(private$tau)
        else private$tau <- value
        },
      #' @field independent Are the levels independent?
      independent    = function(value) private$independent_,
      #' @field distribution Emission distribution of each level
      distribution   = function(value) private$distribution_,
      #' @field directed Are the levels directed?
      directed       = function(value) private$directed_,
      ## other functions
      #' @field entropy Get the entropy of the model
      entropy     = function(value) {
        - sum(private$tau$O * log(private$tau$O)) -
          sum(private$tau$I * log(private$tau$I))
      },
      #' @field bound Get the variational bound of the model
      bound      = function(value) self$complete_likelihood + self$entropy,
      #' @field df_mixture Get the degrees of freedom of the
      #' mixture parameters
      df_mixture = function(value) list(I = private$Q$I -1,
                                        O = private$Q$O -1),
      #' @field df_connect Get the degrees of freedom of the
      #' connection parameters
      df_connect = function(value) {
        list(I = if (private$directed_$I) private$Q$I**2
             else choose(private$Q$I + 1, 2),
             O = if (private$directed_$O) private$Q$O**2
             else choose(private$Q$O + 1, 2))
      },
      #' @field connect Get the number of possible observed connections
      connect    = function(value) {
        list(I = ifelse(private$directed_$I, 1, .5) * sum(private$M$I),
             O = ifelse(private$directed_$O, 1, .5) * sum(private$M$O))
      },
      #' @field ICL Get the ICL model selection criterion of the model
      ICL        = function(value) {
        self$complete_likelihood - self$full_penalty #+ self$entropy
        },
      #' @field full_penalty Get the penalty used to compute the ICL
      full_penalty    = function(value) {
        self$penalty$O + self$penalty$I
        },
      #' @field Z Get the list of block memberships (vector form)
      Z          = function(value) {
        Z = list()
        if (private$Q$I == 1) {
          Z$I = rep(1, private$n$I)
        } else {
          Z$I = apply(private$tau$I, 1, which.max)
        }
        if (private$Q$O == 1) {
          Z$O = rep(1, private$n$O)
        } else {
          Z$O = apply(private$tau$O, 1, which.max)
        }
        return(Z)
      },
      #' @field X_hat Get the list of the matrices of probability conection
      #' predictions
      X_hat = function(value) {
        list(
          I = quad_form(private$param$alpha$I, private$tau$I),
          O = quad_form(private$param$alpha$O, private$tau$O)
          )
        },
      #' @field map Get the list of block memberships (matrix form)
      map = function(value) {
        tmp_map <-
          list(I = 1 * sapply(1:private$Q$O, function(x) self$Z$I %in% x),
               O = 1 * sapply(1:private$Q$I, function(x) self$Z$O %in% x))
        if (private$Q$I == 1)
          tmp_map$I <- matrix(1, nrow = private$n$I, ncol = 1)
        if (private$Q$O == 1)
          tmp_map$O <- matrix(1, nrow = private$n$O, ncol = 1)
        return(tmp_map)
      },
      ##------------------------------------------------------------------------
      ## Computing penalties
      ##------------------------------------------------------------------------
      #' @field penalty Get the ICL penalty
      penalty = function(value) {
          penalty <- list()
          penalty$I <- .5 * private$Q$O * self$df_mixture$I * log(private$n$I) +
            .5 * self$df_connect$I * log(self$connect$I)
          penalty$O <- .5 * self$df_mixture$O * log(private$n$O) +
            .5 * self$df_connect$O * log(self$connect$O)
          return(penalty)
        },
      #' @field likelihood Compute the likelihood of both levels
      likelihood =
        function(value) {
          likelihood <- list()
          factor <-  if (private$directed_$I) 1 else .5
          likelihood$I <-
            factor * (
              sum((private$M$I * private$X$I) *
                    quad_form(log(private$param$alpha$I), private$tau$I)) +
                sum((private$M$I * (1 - private$X$I)) *
                      quad_form(log(1 - private$param$alpha$I), private$tau$I))
            ) +
            sum(private$A * private$tau$I %*%
                  tcrossprod(log(private$param$gamma), private$tau$O))
          factor = if (private$directed_$O) 1 else .5
          likelihood$O <-
            factor * (
              sum((private$M$O * private$X$O) *
                    quad_form(log(private$param$alpha$O),
                              private$tau$O)) +
                sum((private$M$O * (1 - private$X$O)) *
                      quad_form(log(1 - private$param$alpha$O),
                                private$tau$O))) +
            sum(private$tau$O%*%log(private$param$pi$O))
          return(likelihood)
        },
      #' @field  complete_likelihood Get the complete likelihood of the model
      complete_likelihood =
        function(value) self$likelihood$I + self$likelihood$O
      )
  )

#' Extract model coefficients
#'
#' Extracts model coefficients from objects with class \code{\link[=FitMLVSBM]{FitMLVSBM}}
#'
#' @param object an R6 object of class FitMLVSBM
#' @param ... additional parameters for S3 compatibility. Not used
#' @return List of parameters.
#' @export
coef.FitMLVSBM <- function(object, ...) {
  stopifnot(inherits(object, "FitMLVSBM"))
  object$parameters
}

#' Model Predictions
#'
#' Make predictions from an SBM.
#'
#' @param object an R6 object of class \code{\link[=FitMLVSBM]{FitMLVSBM}}
#' @param ... additional parameters for S3 compatibility. Not used
#' @return A list with the following entries:
#' \describe{
#'   \item{dyads}{A list of matrix with the probability of each dyads}
#'   \item{nodes}{A list of vectors with the clustering of each nodes}
#' }
#' @importFrom stats predict
#' @export
predict.FitMLVSBM <- function(object, ...) {
  stopifnot(inherits(object, "FitMLVSBM"))
  list(dyads = object$X_hat,
       nodes = object$Z)
}


#' Multilevel SBM Plot
#'
#' Basic matrix plot method for a FitMLVSBM object
#' @description basic matrix plot method for a FitMLVSBM object
#' @param x an R6 object of class \code{\link[=FitMLVSBM]{FitMLVSBM}}
#' @param type A string for the type of plot, just "matrix" for now
#' @param ... additional parameters for S3 compatibility. Not used
#' @return a ggplot2 object
#' @export
plot.FitMLVSBM = function(x, type = c('matrix'), ...){
  stopifnot(inherits(x, "FitMLVSBM"))
  p <- x$plot(type)
  p
}

