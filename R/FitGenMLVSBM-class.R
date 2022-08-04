#' An R6 Class object, a fitted generalized multilevel network once $dovem() is done
#'
#' @import R6
#' @export
FitGenMLVSBM <-
  R6::R6Class(
    "FitGenMLVSBM",
    private = list(
      n            = NULL, # Number of nodes, a vector of size L
      Q            = NULL, # Number of clusters, a vector of size L
      L            = NULL, # An integer, the number of layers
      param        = NULL, # List of fitted parameters of the model
      tau          = NULL, # List of variational parameters of the fitted model
      A            = NULL, # List of affiliation matrices
      X            = NULL, # List of adjacency matrices
      directed_    = NULL,  # Boolean vector of size L, are the levels directed?
      M            = NULL,  # List of mask matrices
      distribution_ = NULL, # Characters vector of length L, the distribution of X
      independent_ = NULL,
      no_aff       = NULL,   # Boolean vector of size L
      emqr         = NULL, # List of QxQ matrices, number of edges per block
      nmqr         = NULL # List of QxQ matrices, number of nodes per block
    ),
    ##
    ##  Public Methods
    ##
    public = list(
      fit_options = list(ve = "joint"),
      #' @description Constructor for the FitMLVSBM class
      #' @param Q Vector with the number of blocks
      #' @param A List of affiliation matrice
      #' @param X List of adjacency matrices
      #' @param M List of Mask matrices
      #' @param directed Vector of boolean
      #' @param distribution Vector of string
      #' @param independent Boolean
      #' @param no_affiliation A vector of boolean.
      #' For each level, are there any nodes with no affiliations?
      #'
      #' @return A FitGenMLVSBM object
      initialize = function(Q = NULL,
                            A = NULL,
                            X = NULL,
                            M = NULL,
                            directed = NULL,
                            distribution = NULL,
                            independent = FALSE,
                            no_affiliation = NULL) {
        private$A     = A
        private$X     = X
        private$n     = vapply(X, nrow, FUN.VALUE = 1L)
        private$Q     = Q
        private$L     = length(X)
        private$directed_ = directed
        private$distribution_ = distribution
        private$independent_    = independent
        for (m in seq(private$L)) {
          if (is.null(M[[m]])) {
            private$M[[m]] <- diag(-1, private$n[m])
          } else {
            private$M[[m]]          <- M[[m]]
            private$X[[m]][M[[m]] == 0] <- -1
            diag(private$X[[m]])    <- 0
          }
          private$M[[m]] <- if (is.null(M[[m]])) diag(-1, private$n[m]) + 1 else M[[m]]
        }
        if (any(vapply(seq(private$L), function(m) {
          sum(is.na(private$X[[m]])) > 0 } , FUN.VALUE = TRUE))) {
          for (m in seq(private$L)) {
            private$M[[m]][is.na(private$X[[m]])] <- 0
            private$X[[m]][is.na(private$X[[m]])] <- -1
          }
        }
        if (is.null(no_affiliation)) {
          private$no_aff <- c(FALSE, vapply(
            seq_along(A),
            function(m) {
              any(rowSums(A[[m]]) == 0)
            }, FUN.VALUE = TRUE))
        } else {
          private$no_aff <- no_affiliation
        }
        private$emqr <- lapply(seq(private$L), function(m) matrix(NA, Q[m], Q[m]))
        private$nmqr <- lapply(seq(private$L), function(m) matrix(NA, Q[m], Q[m]))
        private$param$pi <- vector("list", private$L)
      },
      #' @field vbound The vector of variational bound for monitoring convergence
      vbound = NULL,
      ##------------------------------------------------------------------------
      ##      ## Parameters update
      ##------------------------------------------------------------------------
      #' @description  Update the connection parameters for the M step
      #' @param safeguard Parameter live in a compact [safeguard, 1-safeguard]
      update_alpha =
        function(m, safeguard = 2*.Machine$double.eps) {
          ## alpha
          if(private$Q[m] == 1) {
            private$param$alpha[[m]] <-
              as.matrix( sum(private$M[[m]] * private$X[[m]])/ sum(private$M[[m]]))
          } else {
            alpha <-
              crossprod(private$tau[[m]],
                        (private$M[[m]]*private$X[[m]]) %*% private$tau[[m]]) /
              crossprod(private$tau[[m]],
                        private$M[[m]] %*% private$tau[[m]])
            private$param$alpha[[m]] <- alpha
          }
          return (private$param$alpha[[m]])
        },
      #' @description Update the  mixture parameter for the M step of level m
      #' @param safeguard Parameter live in a compact [safeguard, 1-safeguard]
      update_pi =
        function(m, safeguard = 1e-3) {
          ## rho
          if (m == 1) {
            if (private$Q[m] == 1) {
              private$param$pi[[m]] <- 1
            } else {
              pi <- colMeans(private$tau[[m]])
              pi[pi < safeguard] <- safeguard
              private$param$pi[[m]] <- pi/sum(pi)
            }
          } else {
            if (private$no_aff[m]) {
              if (private$Q[m] == 1) {
                private$param$pi[[m]] <- 1
              } else {
                pi <- colMeans(private$tau[[m]][rowMeans(private$A[[m-1]]) == 0,])
                pi[pi < safeguard] <- safeguard
                private$param$pi[[m]] <- pi/sum(pi)
              }
            }
          }
          return(private$param$pi[[m]])
        },
      #' @description Update the hierarchical mixture parameter for the M step
      #' @param safeguard Parameter live in a compact [safeguard, 1-safeguard]
      update_gamma =
        function(m, safeguard = 1e-6){
          ## gamma
          if (private$independent_ == FALSE) {
            if (private$Q[m+1] == 1) {
              private$param$gamma[[m]] <-
                matrix(1, nrow = 1, ncol = private$Q[m])
            } else {
              gamma <-
                crossprod(private$tau[[m+1]],
                          private$A[[m]]) %*% private$tau[[m]]
              gamma <-
                t(t(gamma)/colSums(private$A[[m]] %*% private$tau[[m]]))
              gamma[gamma < safeguard] <- safeguard
              private$param$gamma[[m]] <- t(t(gamma)/colSums(gamma))
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
          for (m in seq(private$L)) {
            if (private$Q[m] == 1) {
              private$tau[[m]] <-
                matrix(1, nrow = private$n[m], ncol = 1)
            } else {
              init_clust <-
                switch(method,
                       "spectral"     = spcClust(private$X[[m]], private$Q[m]),
                       "hierarchical" = hierarClust(private$X[[m]], private$Q[m]),
                       "merge_split"  = Z[[m]])
              private$tau[[m]] <-
                1 * sapply(X = seq(private$Q[m]), FUN = function(x) init_clust %in% x)
              private$tau[[m]][private$tau[[m]] < safeguard] <- safeguard
              private$tau[[m]] <-
                private$tau[[m]] / rowSums(private$tau[[m]])
            }
          }
        },
      #' @description Reset all parameters
      clear = function(){
        private$param <- NULL
        private$tau <- NULL
      },
      #-------------------------------------------------------------------------
      # Varational EM algorithm
      #-------------------------------------------------------------------------
      #' @description m_step Compute the M step of the VEM algorithm
      #' @param safeguard Parameter live in a compact [safeguard, 1-safeguard]
      m_step =
        function(m, safeguard = 1e-6){
          self$update_alpha(m, safeguard = safeguard)
          if (m == 1 | private$no_aff[m]) {
            self$update_pi(m, safeguard = safeguard)
          }
          if (m > 1) {
            self$update_gamma(m-1, safeguard = safeguard)
          }
          if (m < private$L) {
            self$update_gamma(m, safeguard = safeguard)
          }
        },
      #' @description Compute the VE step of the VEM algorithm
      #' @param m The level to be updated
      #' @param threshold The convergence threshold
      #' @param fixPointIter The maximum number of fixed point iterations
      #' @param safeguard Parameter live in a compact [safeguard, 1-safeguard]
      ve_step =
        function(m, threshold = 1e-6, fixPointIter = 3, safeguard = 1e-6) {
          condition <- TRUE
          if (private$Q[m] == 1) {
            return (matrix(1, private$n[m], private$Q[m]))
          }
          it        <- 0
          tau_old   <- private$tau[[m]]
          tau       <- private$tau[[m]]
          while (condition) {
            ## sigma
            tau <-
              switch(
                private$distribution_[m],
                "bernoulli" = {
                  tau_new <- (private$M[[m]] * private$X[[m]]) %*%
                    tau_old %*%
                    t(.logit(private$param$alpha[[m]], eps = 1e-9)) +
                    private$M[[m]] %*%
                    tau_old %*%
                    t(.log(1-private$param$alpha[[m]], eps = 1e-9))
                  if (private$directed_[m]) {
                    tau_new <- tau_new +
                      crossprod(private$M[[m]]*private$X[[m]],
                                tau_old %*%
                                  .logit(private$param$alpha[[m]], eps=1e-9)) +
                      crossprod(private$M[[m]],
                                tau_old %*%
                                  .log(1-private$param$alpha[[m]], eps = 1e-9))
                  }
                  invisible(tau_new)
                },
                "poisson" = {
                  tau_new <-
                    (private$M[[m]] * private$X[[m]]) %*%
                    tau_old %*%
                    t(log(private$param$alpha[[m]])) -
                    private$M[[m]] %*%
                    tau_old %*%
                    t(private$param$alpha[[m]])
                  if (private$directed_[m]) {
                    tau_new <- tau_new +
                      crossprod(private$M[[m]]*private$X[[m]],
                                tau_old %*%
                                  log(private$param$alpha[[m]])) -
                      crossprod(private$M[[m]],
                                tau_old %*%
                                  private$param$alpha[[m]])
                  }
                  invisible(tau_new)
                }
              )
            # tau <- density_function(private$X[[m]],
            #                         private$M[[m]],
            #                         tau_old,
            #                         private$direction_[m],
            #                         private$distribution_[m])
            if (private$no_aff[m]) {
              tau <- tau +
                matrix(log(private$param$pi[[m]]),
                       private$n[m],
                       private$Q[m],
                       byrow = TRUE)
            }
            if (m > 1) {
              tau <- tau + private$A[[m-1]] %*%
                tcrossprod(private$tau[[m-1]], log(private$param$gamma[[m-1]]))
            }
            if (m < private$L) {
              tau <- tau + crossprod(private$A[[m]], private$tau[[m+1]]) %*%
                log(private$param$gamma[[m]])
            }
            tau  <- exp(t(apply(X = tau,
                                MARGIN = 1,
                                FUN = function(x) x - max(x))) )
            tau[tau < safeguard] <- safeguard
            tau  <-  tau/rowSums(tau)
            it  <-  it + 1
            condition  <-  (dist_param(tau, tau_old) > threshold &&
                              it <= fixPointIter)
            tau_old   <- tau
          }
          private$tau[[m]] <-  tau
        },

      #' @description Compute the VE step of the VEM algorithm
      #' @param threshold The convergence threshold
      #' @param fixPointIter The maximum number of fixed point iterations
      #' @param safeguard Parameter live in a compact [safeguard, 1-safeguard]
      ve_step2 =
        function(threshold = 1e-6, fixPointIter = 5, safeguard = 1e-6) {
          condition <- TRUE
          tau <- private$tau
          tau_old <- private$tau
          it <- 1
          while (condition) {
          for (m in seq(private$L)) {
            if (private$Q[m] == 1) {
              tau[[m]] <- matrix(1, private$n[m], private$Q[m])
            } else {
              ## sigma
              tau[[m]] <-
                switch(
                  private$distribution_[m],
                  "bernoulli" = {
                    tau_new <- (private$M[[m]] * private$X[[m]]) %*%
                      tau_old[[m]] %*%
                      t(.logit(private$param$alpha[[m]], eps = 1e-9)) +
                      private$M[[m]] %*%
                      tau_old[[m]] %*%
                      t(.log(1-private$param$alpha[[m]], eps = 1e-9))
                    if (private$directed_[m]) {
                      tau_new <- tau_new +
                        crossprod(private$M[[m]]*private$X[[m]],
                                  tau_old[[m]] %*%
                                    .logit(private$param$alpha[[m]], eps=1e-9)) +
                        crossprod(private$M[[m]],
                                  tau_old[[m]] %*%
                                    .log(1-private$param$alpha[[m]], eps = 1e-9))
                    }
                    invisible(tau_new)
                  },
                  "poisson" = {
                    tau_new <-
                      (private$M[[m]] * private$X[[m]]) %*%
                      tau_old[[m]] %*%
                      t(log(private$param$alpha[[m]])) -
                      private$M[[m]] %*%
                      tau_old[[m]] %*%
                      t(private$param$alpha[[m]])
                    if (private$directed_[m]) {
                      tau_new <- tau_new +
                        crossprod(private$M[[m]]*private$X[[m]],
                                  tau_old[[m]] %*%
                                    log(private$param$alpha[[m]])) -
                        crossprod(private$M[[m]],
                                  tau_old[[m]] %*%
                                    private$param$alpha[[m]])
                    }
                    invisible(tau_new)
                  }
                )
            # tau <- density_function(private$X[[m]],
            #                         private$M[[m]],
            #                         tau_old,
            #                         private$direction_[m],
            #                         private$distribution_[m])
            if (private$no_aff[m]) {
              tau[[m]] <- tau[[m]] +
                matrix(log(private$param$pi[[m]]),
                       private$n[m],
                       private$Q[m],
                       byrow = TRUE)
            }
            if (m > 1) {
              tau[[m]] <- tau[[m]] + private$A[[m-1]] %*%
                tcrossprod(tau_old[[m-1]], log(private$param$gamma[[m-1]]))
            }
            if (m < private$L) {
              tau[[m]] <- tau[[m]] + crossprod(private$A[[m]], tau_old[[m+1]]) %*%
                log(private$param$gamma[[m]])
            }
            tau[[m]]  <- exp(t(apply(X = tau[[m]],
                                MARGIN = 1,
                                FUN = function(x) x - max(x))) )
            tau[[m]][tau[[m]] < safeguard] <- safeguard
            tau[[m]]  <-  tau[[m]]/rowSums(tau[[m]])
            }
          }
            it  <-  it + 1
            condition  <-
              (max(vapply(seq_along(tau),
                          function(m) dist_param(tau[[m]], tau_old[[m]]),
                          FUN.VALUE = .1)) > threshold &&
                              it <= fixPointIter)
            tau_old   <- tau
          }
          private$tau <-  tau
        },


      update_mqr = function(m) {
        tau_tmp <- private$tau[[m]]
        private$emqr[[m]] <-
          .tquadform(tau_tmp, private$X[[m]] * private$M[[m]])
        private$nmqr[[m]] <-
          .tquadform(tau_tmp, private$M[[m]])
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
                 maxIter = 1000, fixPointIter = 10,
                 safeguard = 1e-6, Z = NULL) {
          self$init_clustering(method = init, safeguard = safeguard, Z = Z)
          lapply(seq(private$L), function(m) self$m_step(m, safeguard = safeguard))
          lapply(seq(private$L), function(m) self$update_mqr(m))
          self$vbound <-  c(self$vbound, self$bound)
          condition   <-  TRUE
          it          <-  0
          if (any(private$Q >1 )) {
            while (condition) {
              param_old <- private$param
              tau_old   <- private$tau
              bound_old <- self$vbound[length(self$vbound)]
              if (self$fit_options$ve == "sequential") {
                ## Forward pass
                for (m in seq(private$L)) {
                  self$ve_step(m, safeguard = safeguard)
                  self$m_step(m, safeguard = safeguard)
                }
                ## Backward pass
                for (m in seq(private$L-1, 1)) {
                  self$ve_step(m, safeguard = safeguard)
                  self$m_step(m, safeguard = safeguard)
                }
              } else {
                self$ve_step2(safeguard = safeguard)
                for (m in seq(private$L)) self$m_step(m, safeguard = safeguard)
              }

              lapply(seq(private$L), function(m) self$update_mqr(m))
             ## Calculer la vbound par morceau de maniere a ne regarder que la
             ## difference entre les update pour un niveau donné
              if (bound_old > self$bound) {
                private$tau <- tau_old
                private$param     <- param_old
                condition         <- FALSE
              } else {
                it          <-  it + 1
                self$vbound <-  c(self$vbound, self$bound)
                #            cat(it, " : ", self$bound, "\r" )
                condition <-
                  ((max(vapply(
                    seq(private$L),
                    function(m) dist_param(private$param$alpha[[m]],
                                           param_old$alpha[[m]]),
                    FUN.VALUE = .1)))
                   > threshold &
                     it <= maxIter)
                if (it == maxIter) {
                  warning(paste0("VEM for Q = ", private$Q, "did not converge!"))
                }
              }
            }
            invisible(lapply(seq(private$L), self$permute_empty_class))
            #self$show()
          }
        },
      #' @description permute_empty_class Put empty blocks numbers at the end
      permute_empty_class =
        function(m) {
          if (length(unique(self$Z[[m]])) < private$Q[m]) {
            perm  <-  c(unique(self$Z[[m]]), setdiff(seq(private$Q[m]), self$Z[[m]]))
            private$tau[[m]] <-  private$tau[[m]][, perm]
            private$param$alpha[[m]] <-  private$param$alpha[[m]][perm, perm]
            if (private$no_aff[m])
              private$param$pi[[m]]      <-  private$param$pi[[m]][perm]
            if (m > 1) {
              private$param$gamma[[m-1]] <- matrix(data = private$param$gamma[[m-1]][perm,],
                                            nrow = private$Q[m], ncol = private$Q[m-1])
            }
            if (m < private$L) {
              private$param$gamma[[m]]  <-  matrix(data = private$param$gamma[[m]][,perm],
                                                 nrow = private$Q[m+1], ncol = private$Q[m])
            }
          }
        },

      xz_loglikelihood = function(m) {
        factor <-  if (private$directed_[m]) 1 else .5
        switch(
          private$distribution_[m],
          "bernoulli" = {
            alpha <- private$param$alpha[[m]]
            factor * sum(
              .xlogy(private$emqr[[m]], alpha, eps = 1e-12) +
                .xlogy(private$nmqr[[m]] - private$emqr[[m]], 1 - alpha, eps = 1e-12))
          },
          "poisson" = {
            alpha <- private$param$alpha[[m]]
            factor * sum(private$emqr[[m]] * log( alpha ) -
                              private$nmqr[[m]]  * alpha )
          }
        )
      },

      za_loglikelihood = function(m) {
        sum(private$A[[m]] * private$tau[[m+1]] %*%
              tcrossprod(log(private$param$gamma[[m]]), private$tau[[m]]))
      },


    #' @param order One of c("affiliation", "degree")
    #' @description Reorder the block memberships and parameters of the networks
      reorder = function(order = "affiliation") {
      #  browser()
        ord <- list()
        for(l in seq(private$L)) {
          ord[[l]] <- seq(private$Q[l])
        }
        if (order == "degree") {
          for(l in seq(private$L)) {
            if(private$Q[l] > 1) {
              ord[[l]] <- order(self$block_proportions[[l]] %*%
                                  private$param$alpha[[l]], decreasing=TRUE)
            }
          }
          if(private$Q[1] > 1) {
            private$tau[[1]] <- private$tau[[1]][,ord[[1]]]
            private$param$alpha[[1]] <- private$param$alpha[[1]][ord[[1]],ord[[1]]]
            private$param$pi[[1]] <- private$param$pi[[l]][ord[[1]]]
          }
          for(l in seq(2, private$L)) {
            if(private$Q[l] >1) {
              private$tau[[l]] <- private$tau[[l]][,ord[[l]]]
              private$param$alpha[[l]] <- private$param$alpha[[l]][ord[[l]],ord[[l]]]
              if(! is.null(private$param$pi[[l]])) {
                private$param$pi[[l]] <- private$param$pi[[l]][ord[[l]]]
              }
              private$param$gamma[[l-1]] <- private$param$gamma[[l-1]][ord[[l]],ord[[l-1]]]
            }
          }
        }
        if (order == "affiliation") {
          if(private$Q[1] > 1) {
            ord[[1]] <- order(private$param$pi[[1]] %*% private$param$alpha[[1]],
                              decreasing=TRUE)
            private$tau[[1]] <- private$tau[[1]][,ord[[1]]]
            private$param$alpha[[1]] <- private$param$alpha[[1]][ord[[1]],ord[[1]]]
            private$param$pi[[1]] <- private$param$pi[[1]][ord[[1]]]
          }
          for(l in seq(2, private$L)) {
            if(private$Q[l] >1) {
              ord[[l]] <- order(apply(private$param$gamma[[l-1]][,ord[[l-1]]], 1, "which.max"),
                                decreasing = FALSE)
              private$tau[[l]] <- private$tau[[l]][,ord[[l]]]
              private$param$alpha[[l]] <- private$param$alpha[[l]][ord[[l]],ord[[l]]]
              if(! is.null(private$param$pi[[l]])) {
                private$param$pi[[l]] <- private$param$pi[[l]][ord[[l]]]
              }
              private$param$gamma[[l-1]] <- private$param$gamma[[l-1]][ord[[l]],ord[[l-1]]]
            }
          }
        }
        return(ord)
      },
      #' @description Plot of FitMLVSBM objects
      #' @param type A string for the type of plot, just "matrix" for now
      #' @return a ggplot2 object
      plot = function(type = c('matrix'), ...) {
        #if(type == "matrix") {
          if (! requireNamespace("ggplot2", quietly = TRUE)) {
            stop("Please install ggplot2: install.packages('ggplot2')")
          }
          if (! requireNamespace("patchwork", quietly = TRUE)) {
            stop("Please install ggplot2: install.packages('patchwork')")
          }
          if (! requireNamespace("RColorBrewer", quietly = TRUE)) {
            stop("Please install ggplot2: install.packages('RColorBrewer')")
          }
          # if (! requireNamespace("ggraph", quietly = TRUE)) {
          #   stop("Please install ggraph: install.packages('ggraph')")
          # }
          # if (! requireNamespace("tidygraph", quietly = TRUE)) {
          #   stop("Please install tidygraph: install.packages('tidygraph')")
          # }
          p <- plot_generalized_multilevel_graphon(self, ...)
     #   }
        p
      },
      #' @description print method
      #' @param type character to tune the displayed name
      show = function(type = "Multilevel Stochastic Block Model") {
        cat(type, "--", self$distribution, "variant\n")
        cat("=====================================================================\n")
        cat("Dimension = (", self$nb_nodes, ") - (",
            self$nb_clusters,  ") blocks.\n")
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
      #' @field nb_levels Get the number of levels
      nb_levels          = function(value)
        if (missing(value)) private$L else private$L = value,
      #' @field block_proportions Get the block proportions of each level
      block_proportions          = function(value)
        lapply(private$tau, function(tau) colMeans(tau)),
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
        - sum(vapply(
          X = seq(private$L),
          FUN = function(m) sum(.xlogx(private$tau[[m]])),
              FUN.VALUE = .1))
      },
      #' @field bound Get the variational bound of the model
      bound      = function(value) self$complete_likelihood + self$entropy,
      #' @field df_mixture Get the degrees of freedom of the
      #' mixture parameters
      df_mixture = function(value) {
        vapply(seq(private$L),
               function(m) private$Q[m] -1,
               FUN.VALUE = 1)
      },
      #' @field df_connect Get the degrees of freedom of the
      #' connection parameters
      df_connect = function(value) {
        vapply(seq(private$L),
               function(m) {
                 if (private$directed_[m]) private$Q[m]**2
                 else choose(private$Q[m] + 1, 2)
                 },
               FUN.VALUE = 1)
      },
      #' @field connect Get the number of possible observed connections
      connect    = function(value) {
        vapply(seq(private$L),
               function(m) {
                 ifelse(private$directed_[m], 1, .5) * sum(private$M[[m]])
               },
               FUN.VALUE = 1)
      },
      #' @field ICL Get the ICL model selection criterion of the model
      ICL        = function(value) {
        self$complete_likelihood - self$full_penalty #+ self$entropy
      },
      #' @field full_penalty Get the penalty used to compute the ICL
      full_penalty    = function(value) {
        sum(self$penalty)
      },
      #' @field Z Get the list of block memberships (vector form)
      Z          = function(value) {
        Z <- lapply(seq(private$L),
               function(m) {
                 if (private$Q[m] == 1) {
                   return(rep(1, private$n[m]))
                 } else {
                   apply(private$tau[[m]], 1, which.max)
                 }
               }
        )
        return(Z)
      },
      #' @field X_hat Get the list of the matrices of probability connection
      #' predictions
      X_hat = function(value) {
        lapply(seq(private$L),
               quad_form(private$param$alpha[[m]], private$tau[[m]]) *
                 (1-diag(1, private$n[m])))
      },
      #' @field map Get the list of block memberships (matrix form)
      map = function(value) {
        tmp_map <-
          lapply(seq(private$L),
                 function(m) {
                   if (private$Q[m] == 1) {
                     return (matrix(1, nrow = private$n[m], ncol = 1))
                   } else {
                     return (1 * vapply(seq(private$Q[m]),
                                        function(x) self$Z[[m]] %in% x,
                                        FUN.VALUE = TRUE))
                   }
                 }
                 )
        return(tmp_map)
      },
      ##------------------------------------------------------------------------
      ## Computing penalties
      ##------------------------------------------------------------------------
      #' @field penalty Get the ICL penalty
      penalty = function(value) {
        penalty <- vapply(
          seq(private$L),
          function(m) {
            if (m == 1) {
              pen <- .5 * self$df_mixture[m] * log(private$n[m]) +
                .5 * self$df_connect[m] * log(self$connect[m])
            } else {
              pen <-  .5 * self$df_connect[m] * log(self$connect[m])
              pen <-
                if (private$no_aff[m]) {
                  nb_noaff <- sum(rowSums(private$A[[m-1]]) == 0)
                  pen <- pen + .5 * self$df_mixture[m] * log(nb_noaff)
                  + .5 * private$Q[m-1] * self$df_mixture[m] *
                    log(private$n[m] - nb_noaff)
                } else {
                  pen <- pen + .5 * private$Q[m-1] * self$df_mixture[m] *
                    log(private$n[m])
                }
            }
            return(pen)
          },
          FUN.VALUE = .1)
        return(penalty)
      },
      #' @field likelihood Compute the likelihood of both levels
      likelihood =
        function(value) {
          likelihood <-
            vapply(
              seq(private$L),
              function(m) {
                ll <-  self$xz_loglikelihood(m)
                if (m == 1) {
                  ll <- ll + sum(private$tau[[m]]%*%log(private$param$pi[[m]]))
                } else {
                  ll <-  ll + self$za_loglikelihood(m-1)
                  if (private$no_aff[m]) {
                    ll <- ll + sum(private$tau[[m]]%*%log(private$param$pi[[m]]))
                  }
                }
                return(ll)
              },
              FUN.VALUE = .1
            )
          return(likelihood)
        },
      #' @field  complete_likelihood Get the complete likelihood of the model
      complete_likelihood =
        function(value) sum(self$likelihood)
    )
  )

#' Extract model coefficients
#'
#' Extracts model coefficients from objects with class
#' \code{\link[=FitGenMLVSBM]{FitGenMLVSBM}}
#'
#' @param object an R6 object of class FitMLVSBM
#' @param ... additional parameters for S3 compatibility. Not used
#' @return List of parameters.
#' @export
coef.FitGenMLVSBM <- function(object, ...) {
  stopifnot(inherits(object, "FitGenMLVSBM"))
  object$parameters
}

#' Model Predictions
#'
#' Make predictions of dyads or give the nodes memberships for all levels
#' with a generalized multilevel SBM.
#'
#' @param object an R6 object of class \code{\link[=FitGenMLVSBM]{FitGenMLVSBM}}
#' @param ... additional parameters for S3 compatibility. Not used
#' @return A list with the following entries:
#' \describe{
#'   \item{dyads}{A list of matrix with the probability of each dyads}
#'   \item{nodes}{A list of vectors with the clustering of each nodes}
#' }
#' @importFrom stats predict
#' @export
predict.FitGenMLVSBM <- function(object, ...) {
  stopifnot(inherits(object, "FitGenMLVSBM"))
  list(dyads = object$X_hat,
       nodes = object$Z)
}


#' Generalized Multilevel SBM Plot
#'
#' Basic matrix plot method for a FitMLVSBM object
#' @description basic matrix plot method for a FitMLVSBM object
#' @param x an R6 object of class \code{\link[=FitGenMLVSBM]{FitGenMLVSBM}}
#' @param type A string for the type of plot, just "matrix" for now
#' @param ... additional parameters for S3 compatibility. Not used
#' @return a ggplot2 object
#' @export
plot.FitGenMLVSBM <- function(x, type = c('graphon'), ...){
  stopifnot(inherits(x, "FitGenMLVSBM"))
  p <- x$plot(type, ...)
  p
}

