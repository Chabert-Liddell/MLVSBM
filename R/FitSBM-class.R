#' An R6 Class ocject, a fitted level of a multilelvel network once $do_vem() is done
#'
#' @import R6
#' @export


#-------------------------------------------------------------------------------
# Class declaration
#-------------------------------------------------------------------------------

FitSBM <-
  R6::R6Class(
    "FitSBM",
    private = list(
      Q         = NULL, # Number of clusters
      pi        = NULL, # Mixture parameter
      alpha     = NULL, # Connectivity parameter
      tau       = NULL, # Variational parameter
      X         = NULL, # Adjacency Matrix
      n         = NULL, # number of nodes
      directed_  = NULL, # is X directed ?
      M         = NULL,  # mask matrix,
      distribution_ = NULL
    ),
    public = list(
      ## constructor
      initialize = function(Q = 1,
                            X = NULL,
                            M = NULL,
                            directed = FALSE,
                            distribution = "bernoulli") {
        private$X        = X
        private$n        = nrow(X)
        private$Q        = Q
        private$directed_ = directed
        private$distribution_ = distribution
        if (is.null(M)) {
          private$M  <-  diag(-1, private$n) + 1
        } else {
          private$M          <-  M
          private$X[M == 0]  <-  -1
          diag(private$X)    <-  0
        }
        if (sum(is.na(X)) > 0) {
          private$M[is.na(X)] <- 0
          private$X[is.na(X)] <- -1
        }
      },
      vbound = NULL
    ),
    active = list(
      ## accessor and mutator
      #' @field adjacency Get the adjacency matrix
      adjacency  = function(value) private$X,
      #' @field mask Get the mask matrix for dealing with NA
      mask = function(value) private$M,
      #' @field nb_nodes Get the number of nodes of the level
      nb_nodes   = function(value) private$n,
      #' @field nb_clusters Get the number of blocks
      nb_clusters = function(value) private$Q,
      #' @field distribution Get the distribution used for the connections
      distribution = function(value) private$distribution_,
      #' @field directed Get if the level is directed or not
      directed = function(value) private$directed_,
      #' @field mixture_parameter Access the block proportions
      mixture_parameter = function(value) {
        if (missing(value)) return(private$pi) else private$pi <- value
      },
      #' @field connectivity_parameter Access the connectivity matrix
      connectivity_parameter = function(value) {
        if (missing(value)) return(private$alpha) else private$alpha <- value
      },
      #' @field membership Access the variational parameters
      membership = function(value) {
        if (missing(value)) return(private$tau) else private$tau <- value
      },
      ## other functions
      #' @field entropy Get the entropy of the model
      entropy    = function(value) - sum( xlogx(private$tau)),
      #' @field bound Get the variational bound of the model
      bound     = function(value) self$likelihood + self$entropy,
      #' @field df_mixtures Get the degree of freedom of the block proportion
      df_mixture = function(value) private$Q -1,
      #' @field df_connect Get the degree of freedom of the connection parameters
      df_connect = function(value) {
        if (private$directed_) private$Q**2 else private$Q*(private$Q+1)/2
      },
      connect = function(value) ifelse (private$directed_, 1, .5)*sum(private$M),
      #' @field ICL Get the ICL model selection criterion
      ICL        = function(value) self$likelihood - self$penalty,
      #' @field penalty Get the penalty used for computing the ICL
      penalty    = function(value) {
        .5*self$df_connect*log(self$connect) + .5*self$df_mixture*log(private$n)
      },
      #' @field Z Access the vector of block membership (clustering)
      Z          = function(value){
        if (private$Q == 1) rep(1, private$n) else apply(private$tau, 1, which.max)
      },
      #' @field X_hat Get the connection probability matrix
      X_hat = function(value) {
        quad_form(private$alpha,private$tau)
      }
    )
  )
#-------------------------------------------------------------------------------
# Parameters update
#-------------------------------------------------------------------------------
FitSBM$set(
  "public",
  "update_alpha",
  function(safeguard = 1e-6) {
    ## alpha
    if (private$Q == 1) {
      private$alpha <-
        as.matrix( sum(private$M * private$X) / sum(private$M))
      } else {
        alpha <- crossprod(private$tau, (private$M*private$X) %*% private$tau) /
          crossprod(private$tau, private$M %*% private$tau)
        private$alpha <- alpha
        }
    return(private$alpha)
    }
  )
FitSBM$set("public", "update_pi",
        function(safeguard = 1e-6){
          ## pi
          if (private$Q == 1) {
            private$pi = 1
          } else {
            pi                 <- colMeans(private$tau)
            pi[pi < safeguard] <- safeguard
            private$pi         <- pi/sum(pi)
          }
          return(private$pi)
        })
FitSBM$set("active", "likelihood",
        function(value) self$X_likelihood + self$Z_likelihood
)
#-------------------------------------------------------------------------------
#  Inference
#-------------------------------------------------------------------------------
FitSBM$set(
  "public",
  "init_clustering",
  function(safeguard = 1e-6,
           method = "hierarchical",
           Z = NULL) {
    if (private$Q == 1) {
      private$tau <-
        as.matrix(rep(1, private$n), nrow = private$n, ncol = 1)
    } else {
      init_clustering <-
        switch(method,
               "spectral"     = spcClust(private$X, private$Q),
               "hierarchical" = hierarClust(private$X, private$Q),
               "merge_split"  = Z)
      private$tau <-
        1 * sapply(1:private$Q, function(x) init_clustering %in% x)
      private$tau[private$tau < safeguard] <- safeguard
      private$tau <- private$tau / rowSums(private$tau)
    }
  })




#-------------------------------------------------------------------------------
# Likelihood computation
#-------------------------------------------------------------------------------
FitSBM$set("active", "X_likelihood",
        function(value){
          facteur <-  if (private$directed_) 1 else .5
          return(
            facteur * (
              sum(private$M * private$X *
                    quad_form(log(private$alpha), private$tau)) +
                sum(private$M * (1 - private$X) *
                      quad_form(log(1 - private$alpha), private$tau))
            )
          )
        }
)
FitSBM$set("active", "Z_likelihood",
        function(value) sum(private$tau%*%log(private$pi))
)


#-------------------------------------------------------------------------------
# Varational EM algorithm
#-------------------------------------------------------------------------------
FitSBM$set("public", "m_step",
        function(safeguard = 1e-6){
          self$update_alpha()
          self$update_pi()
        })

FitSBM$set(
  "public",
  "ve_step",
  function(threshold = 1e-6, fixPointIter = 100, safeguard = 1e-6){
    condition <-  TRUE
    it        <-  0
    tau_old   <-  private$tau
    while(condition){
      ## tau
      tau  <-
        matrix(log(private$pi), private$n, private$Q, byrow = TRUE) +
        (private$M * private$X) %*% tcrossprod(tau_old,
                                               log(private$alpha)) +
        (private$M * (1 - private$X)) %*% tcrossprod(tau_old,
                                                     log(1 - private$alpha))
      if (private$Q == 1) {
        tau  <- as.matrix(exp(apply(tau, 1, x <- function(x) x - max(x))), ncol = 1 )
        } else {
          tau  <-  exp(t(apply(tau, 1, x <- function(x) x - max(x))) )
          }
      tau[tau < safeguard] = safeguard
      tau <-  tau/rowSums(tau)
      it <-  it + 1
      condition  <- dist_param(tau, tau_old) > threshold && it <= fixPointIter
      tau_old   <- tau
      }
    private$tau <- tau
    }
  )

FitSBM$set(
  "public",
  "do_vem",
  function(init = "hierarchical", threshold = 1e-6, maxIter = 1000,
           fixPointIter = 100, safeguard = 1e-6, Z = NULL,
           bound = NA) {
    self$init_clustering(method = init, safeguard = safeguard, Z = Z)
    self$m_step(safeguard = safeguard)
    self$vbound <-  c(self$vbound, self$bound)
    condition   <-  TRUE
    it          <-  0
    if (private$Q != 1) {
      while (condition) {
        alpha_old <- private$alpha
        pi_old    <- private$pi
        tau_old   <- private$tau
        bound_old <- self$vbound[length(self$vbound)]
        self$ve_step(safeguard = safeguard)
        self$m_step(safeguard = safeguard)
        if (bound_old > self$bound) {
          private$tau   <- tau_old
          private$alpha <- alpha_old
          private$pi    <- pi_old
          condition     <- FALSE
          } else {
            it          <-  it + 1
            self$vbound <-  c(self$vbound, self$bound)
#            cat(it, " : ", self$bound, "\r" )
            condition <- dist_param(private$alpha, alpha_old) >
              (threshold && it <= maxIter)
          }
        }
      self$permute_empty_class()
      }
    }
  )

FitSBM$set("public", "permute_empty_class",
        function(){
          if(length(unique(self$Z)) < private$Q){
            perm  <-  c(unique(self$Z), setdiff( 1:private$Q, self$Z))
            private$tau     <-  private$tau[, perm]
            private$alpha         = private$alpha[perm, perm]
            private$pi            = private$pi[perm]
          }
        }
)

FitSBM$set("public", "clear",
        function(){
          private$pi     <-  NULL
          private$alpha  <-  NULL
          private$tau    <-  NULL
        })

