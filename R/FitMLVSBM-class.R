#-------------------------------------------------------------------------------
# Class declaration
#-------------------------------------------------------------------------------

#' @import R6
#' @export
FitMLVSBM <-
  R6::R6Class(
    "FitMLVSBM",
    private = list(
      n         = NULL, # number of nodes in I and O
      Q         = NULL, # Number of clusters in I and O
      param     = NULL, # List of fitted parameters of the model
      var_param = NULL, # List of variational parameters of the fitted model
      A         = NULL, # Affiliation Matrix
      X         = NULL, # List of adjacency Matrices
      directed  = NULL,  # is XI directed ? Is XO directed ?
      M         = NULL,  # List of mask matrix
      distribution = NULL # List of the distribution of X
    ),
    public = list(
      ## constructor
      initialize = function(Q = list(I = 1, O = 1),
                            A = NA, X = NA,
                            M = NA,
                            directed = NA,
                            distribution = list("bernoulli", "bernoulli")) {
        private$A     = A
        private$X     = X
        private$n     = list(I = nrow(A), O = ncol(A))
        private$Q     = Q
        private$directed = directed
        distribution = distribution
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
      vbound = NULL
    ),
    active = list(
      ## accessor and mutator
      affiliation_matrix = function(value) private$A,
      adjacency_matrix   = function(value) private$X,
      nb_nodes           = function(value) private$n,
      nb_clusters        = function(value) private$Q,
      parameters = function(value) {
        if(missing(value)) return(private$param)
        else private$param <- value
        },
      membership = function(value) {
        if(missing(value)) return(private$var_param)
        else private$var_param <- value
        },
      ## other functions
      entropy     = function(value) {
        - sum(private$var_param$tau$O * log(private$var_param$tau$O)) -
          sum(private$var_param$tau$I * log(private$var_param$tau$I))
      },
      bound      = function(value) self$complete_likelihood + self$entropy,
      df_mixture = function(value) list(I = private$Q$I -1,
                                        O = private$Q$O -1),
      df_connect = function(value) {
        list(I = if (private$directed$I) private$Q$I**2
             else choose(private$Q$I + 1, 2),
             O = if (private$directed$O) private$Q$O**2
             else choose(private$Q$O + 1, 2))
      },
      connect    = function(value) {
        list(I = dplyr::if_else(private$directed$I, 1, .5) * sum(private$M$I),
             O = dplyr::if_else(private$directed$O, 1, .5) * sum(private$M$O))
      },
      ICL        = function(value) {
        self$complete_likelihood + self$entropy - self$full_penalty
        },
      full_penalty    = function(value) {
        self$penalty$O + self$penalty$I
        },
      Z          = function(value) {
        Z = list()
        if (private$Q$I == 1) {
          Z$I = rep(1, private$n$I)
        } else {
          Z$R = apply(private$var_param$tau$I, 1, which.max)
        }
        if (private$Q$O == 1) {
          Z$L = rep(1, private$n$O)
        } else {
          Z$L = apply(private$var_param$tau$O, 1, which.max)
        }
        return(Z)
      }
      X_hat = function(value) {
        list(
          I = quad_form(private$param$alpha$I, private$var_param$tau$I),
          O = quad_form(private$param$alpha$O, private$var_param$tau$O)
          )
        },
      map = function(value) {
        tmp_map <-
          list(I = 1 * sapply(1:private$Q$O, function(x) self$Z$I %in% x),
               O = 1 * sapply(1:private$Q$I, function(x) self$Z$O %in% x))
        if (private$Q$I == 1)
          tmp_map$I <- matrix(1, nrow = private$n$I, ncol = 1)
        if (private$Q$O == 1)
          tmp_map$O <- matrix(1, nrow = private$n$O, ncol = 1)
        return(tmp_map)
      }
    )
  )

##------------------------------------------------------------------------------
## Computing penalties
##------------------------------------------------------------------------------
FitMLVSBM$set(
  "active",
  "penalty",
  function(value) {
  penalty <- list()
  penalty$I <- .5 * private$Q$O * self$df_mixture$I * log(private$n$I) +
    .5 * self$df_connect$I * log(self$connect$I)
  penalty$O <- .5 * self$df_mixture$O * log(private$n$O) +
    .5 * self$df_connect$O * log(self$connect$O)
  }
  )


#-------------------------------------------------------------------------------
# Likelihood computation
#-------------------------------------------------------------------------------
FitMLVSBM$set(
  "active",
  "likelihood",
  function(value) {
    likelihood <- list()
    factor <-  if (private$directed$I) 1 else .5
    likelihood$I <-
      factor * (
        sum((private$M$I * private$X$I) *
              quad_form(log(private$param$alpha$I), private$var_param$tau$I)) +
          sum((private$M$I * (1 - private$X$I)) *
                quad_form(log(1 - private$param$alpha), private$var_param$tau$I))
        ) +
      sum(private$A * private$var_param$tau$I %*%
            tcrossprod(log(private$param$gamma), private$var_param$tau$O))
    factor = if (private$directed$O) 1 else .5
    likelihood$O <-
      factor * (
        sum((private$M$O * private$X$O) *
              quad_form(log(private$param$alpha$O),
                        private$var_param$tau$O)) +
          sum((private$M$O * (1 - private$X$O)) *
                quad_form(log(1 - private$param$alpha$O),
                          private$var_param$tau$O))) +
      sum(private$var_param$tau%*%log(private$param$rho))
    }
  )

FitMLVSBM$set("active", "complete_likelihood",
             function(value) self$likelihood$I + self$likelihood$O
)


##------------------------------------------------------------------------------
## Parameters update
##------------------------------------------------------------------------------
FitMLVSBM$set(
  "public",
  "update_alpha_I",
  function(safeguard = 2*.Machine$double.eps) {
    ## alpha
    if(private$Q$I == 1) {
      private$param$alpha$I <-
        as.matrix( sum(private$M$I * private$X$I)/ sum(private$M$I))
      } else {
        alpha <-
          crossprod(private$var_param$tau$I,
                    (private$M$I*private$X$I) %*% private$var_param$tau$I) /
          crossprod(private$var_param$tau$I,
                    private$M$I %*% private$var_param$tau$I)
        private$param$alpha$I <- alpha
      }

    return (private$param$alpha$I)
    }
  )
FitMLVSBM$set(
  "public",
  "update_alpha_O",
  function(safeguard = 2*.Machine$double.eps) {
    if(private$Q$O == 1) {
      private$param$alpha$O <-
        as.matrix( sum(private$X$O) / sum(private$M$O))
      } else {
        alpha <-
          crossprod(private$var_param$tau$O,
                    (private$M$O * private$X$O) %*% private$var_param$tau$O) /
          crossprod(private$var_param$tau$O,
                    private$M$O %*% private$var_param$tau$O)
        private$param$alpha$O <- alpha
               }
    return (private$param$alpha$O)
    }
  )
FitMLVSBM$set(
  "public",
  "update_pi_O",
  function(safeguard = 1e-2) {
    ## rho
    if (private$Q$O == 1) {
      private$param$pi$O = 1
      } else {
        pi <- colMeans(private$var_param$tau$O)
        pi[pi < safeguard] <- safeguard
        private$param$pi$O <- pi/sum(pi)
        }
    return(private$param$pi$O)
    }
  )
FitMLVSBM$set(
  "public",
  "update_gamma",
  function(safeguard = 1e-6){
    ## gamma
    if (private$name == "classic") {
      if (private$Q$I == 1) {
        private$param$gamma <-
          matrix(rep(1, private$Q$O), nrow = 1, ncol = private$Q$O)
        } else {
          gamma <-
            crossprod(private$var_param$tau$I,
                      private$A) %*% private$var_param$tau
          gamma <-
            t(t(gamma)/colSums(private$A %*% private$var_param$tau))
          gamma[gamma < safeguard] <- safeguard
          private$param$gamma <- t(t(gamma)/colSums(gamma))
          }
      return(private$param$gamma)
      }
    if (private$name == "independent") {
      if (private$Q$I == 1) {
        private$param$gamma <-
          matrix(1, nrow = 1, ncol = private$Q$O)
        } else {
          gamma <- matrix(colMeans(private$var_param$tau$I),
                          nrow = private$Q$I,
                          ncol = private$Q$O)
                       gamma[gamma < safeguard] <- safeguard
                       private$param$gamma <- t(t(gamma)/colSums(gamma))
                     }
                     return(private$param$gamma)
                   }
                 })


#-------------------------------------------------------------------------------
#  Inference
#-------------------------------------------------------------------------------
SBMModel$set(
  "public",
  "initClustering",
  function(safeguard = 2*.Machine$double.eps,
           method = "hierarchical",
           Z = NULL) {
    if (private$Q == 1) {
      private$var_param$tau$I <-
        as.matrix(rep(1, private$n$I), nrow = private$n$I, ncol = 1)
    } else {
      RinitClust <-
        switch(method,
               "spectral"     = spcClust(private$X$I, private$Q),
               "hierarchical" = hierarClust(private$X$I, private$Q),
               "merge_split"  = Z$R)
      private$var_param$tau$I <-
        1 * sapply(1:private$Q, function(x) RinitClust %in% x)
      private$var_param$tau$I[private$var_param$tau$I < safeguard] <-
        safeguard
      private$var_param$tau$I <-
        private$var_param$tau$I / rowSums(private$var_param$tau$I)
    }
    if(private$S == 1){
      private$var_param$tau$O <-
        as.matrix(rep(1, private$n$O), nrow = private$n$O, ncol = 1)
    }else{
      LinitClust <-
        switch(method,
               "spectral"     = spcClust(private$X$O, private$S),
               "hierarchical" = hierarClust(private$X$O, private$S),
               "merge_split"  = Z$L)
      private$var_param$tau$O = 1 * sapply(1:private$S, function(x) LinitClust %in% x)
      private$var_param$tau$O[private$var_param$tau$O < safeguard] = safeguard
      private$var_param$tau$O = private$var_param$tau$O/rowSums(private$var_param$tau$O)
    }
  })



SBMModel$set("public", "clear",
             function(){
               private$param = NULL
               private$var_param = NULL
             })





#-------------------------------------------------------------------------------
# Varational EM algorithm
#-------------------------------------------------------------------------------
SBMModelStoc$set("public", "MStep",
                 function(safeguard = 1e-6){
                   self$update_alpha(safeguard = safeguard)
                   self$update_beta(safeguard = safeguard)
                   self$update_rho(safeguard = safeguard)
                   if (private$name != "deterministic") self$update_gamma(safeguard = safeguard)
                 })

SBMModelStoc$set("public", "VEStep",
                 function(threshold = 1e-6, fixPointIter = 100, safeguard = 1e-6){
                   condition = TRUE
                   it        = 0
                   tau_old   = private$var_param$tau
                   sigma_old = private$var_param$tau$I
                   while(condition){
                     ## sigma
                     sigma =
                       (private$M$I * private$X$I) %*%
                       tcrossprod(sigma_old, log(private$param$alpha)) +
                       (private$M$I * (1 - private$X$I)) %*%
                       tcrossprod(sigma_old, log(1 - private$param$alpha)) +
                       private$A %*% tcrossprod(tau_old, log(private$param$gamma))
                     if(private$Q == 1){
                       sigma = as.matrix(exp( apply(sigma, 1, x <- function(x) x - max(x))), ncol = 1 )
                     }else{
                       sigma = exp( t(apply(sigma, 1, x <- function(x) x - max(x))) )
                     }
                     sigma[sigma < safeguard] = safeguard
                     sigma = sigma/rowSums(sigma)
                     ## tau
                     tau =
                       matrix(log(private$param$rho), private$m, private$S, byrow = TRUE) +
                       (private$M$O * private$X$O) %*%
                       tcrossprod(tau_old, log(private$param$alpha$O)) +
                       (private$M$O * (1 - private$X$O)) %*%
                       tcrossprod(tau_old,  log(1 - private$param$alpha$O)) +
                       crossprod(private$A, sigma_old) %*% log(private$param$gamma)
                     if(private$S == 1){
                       tau = as.matrix(exp(apply(tau, 1, x <- function(x) x - max(x))), ncol = 1 )
                     }else{
                       tau = exp(t(apply(tau, 1, x <- function(x) x - max(x))) )
                     }
                     tau[tau < safeguard] = safeguard
                     tau = tau/rowSums(tau)
                     it = it + 1
                     # condition = (max(abs(tau - tau_old)) > threshold &&
                     #                 it <= fixPointIter)
                     condition = (max(dist_param(tau, tau_old),
                                      dist_param(sigma, sigma_old)) > threshold &&
                                    it <= fixPointIter)
                     tau_old   <- tau
                     sigma_old <- sigma
                   }
                   private$var_param$tau = tau
                   private$var_param$tau$I = sigma
                 })

SBMModelStoc$set("public", "doVEM",
                 function(init = "hierarchical", threshold = 1e-6, maxIter = 1000,
                          fixPointIter = 100, safeguard = 1e-6, Z = NULL,
                          bound = NA){
                   self$initClustering(method = init, safeguard = safeguard, Z = Z)
                   self$MStep(safeguard = safeguard)
                   self$vbound <-  c(self$vbound, self$bound)
                   condition   <-  TRUE
                   it          <-  0
                   if (private$Q != 1 | private$S!= 1) {
                     while (condition) {
                       param_old <- private$param
                       var_param_old   <- private$var_param
                       bound_old <- self$vbound[length(self$vbound)]
                       self$VEStep(safeguard = safeguard)
                       self$MStep(safeguard = safeguard)
                       if (bound_old > self$bound) {
                         private$var_param <- var_param_old
                         private$param     <- param_old
                         condition         <- FALSE
                       } else {
                         it          <-  it + 1
                         self$vbound <-  c(self$vbound, self$bound)
                         cat(it, " : ", self$bound, "\r" )
                         # condition = (max(c(max(abs(private$param$alpha - alpha_old)),
                         #                     max(abs(private$param$alpha$O  - beta_old)),
                         #                     max(abs(private$param$gamma  - gamma_old)))) > threshold &&
                         #                 it <= maxIter)
                         condition <- (max(c(dist_param(private$param$alpha, param_old$alpha),
                                             dist_param(private$param$alpha$O, param_old$alpha$O),
                                             dist_param(private$param$gamma, param_old$gamma)))
                                       > threshold &&
                                         it <= maxIter)
                       }
                     }
                     self$permute_empty_class()
                   }
                 })

SBMModelStoc$set("public", "permute_empty_class",
                 function(){
                   if(length(unique(self$Z$R)) < private$Q){
                     perm = c(unique(self$Z$R),
                              setdiff( 1:private$Q, self$Z$R))
                     private$var_param$tau$I     = private$var_param$tau$I[, perm]
                     private$param$alpha         = private$param$alpha[perm, perm]
                     private$param$gamma         = matrix(private$param$gamma[perm,], private$Q, private$S)
                   }
                   if(length(unique(self$Z$L)) < private$S){
                     perm = c(unique(self$Z$L),
                              setdiff( 1:private$S, self$Z$L))
                     private$var_param$tau  = private$var_param$tau[, perm]
                     private$param$beta     = private$param$beta[perm, perm]
                     private$param$rho      = private$param$rho[perm]
                     private$param$gamma    = matrix(private$param$gamma[,perm], private$Q, private$S)
                   }
                 }
)
