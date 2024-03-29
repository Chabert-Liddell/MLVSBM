#-------------------------------------------------------------------------------
#  SBM inference for upper or lower level
#-------------------------------------------------------------------------------
MLVSBM$set(
  "public",
  "estimate_level",
#' Fit a SBM on a given level of a MLVSBM object
#'
#' @param level "lower" or "upper"
#' @param Q_min Lower bound of exploration space
#' @param Q_max Upper bound of exploration space
#' @param Z Initial clustering if any
#' @param init Initialization method one of "hierarchical",
#' "spectral" or "merged_split"
#' @param depth Exploration depth (not used)
#' @param nb_cores Number of cores for parallel computing
#'
#' @return A list of FitSBM objects
  function(level = "lower",
           Q_min = 1,
           Q_max = 10,
           Z     = NULL,
           init = "hierarchical",
           depth = 1,
           nb_cores = NULL) {
    if (is.null(Z)) {
      model_list <-  vector("list", Q_max)
    } else {
      model_list <-  vector("list", max(max(Z), Q_max))
    }
    did_spectral <- rep(FALSE, Q_max)
    ICL <- rep(-Inf, length(model_list))
    bound <- rep(-Inf, length(model_list))
    print(paste0("Infering ", level, " level :"))
    # Ascending
    if (! is.null(Z)) {
      init <- "merge_split"
      best_fit <- self$estimate_sbm(level = level, Z = Z,
                                    Q = max(Z), init = init)
      model_list[[best_fit$nb_clusters]] <- best_fit
      ICL[best_fit$nb_clusters] <- best_fit$ICL
      bound[best_fit$nb_clusters] <- best_fit$bound
    } else {
      best_fit <- self$estimate_sbm(level = level, Z = Z,
                                      Q = Q_min, init = init)
      # if (! did_spectral[Q_min]) {
      #   spc_fit <- self$estimate_sbm(level = level, Z = Z,
      #                                Q = Q_min, init = "")
      #   if (spc_fit$ICL > best_fit$ICL) {
      #
      #   }
      # }
      model_list[[Q_min]] <- best_fit
      ICL[Q_min] <- best_fit$ICL
      bound[Q_min] <- best_fit$bound
    }
    condition <- TRUE
    print(paste0("# blocks: ", best_fit$nb_clusters,
                 ", ICL = ", best_fit$ICL, " !" ))
    while (condition) {
      fits <- self$estimate_sbm_neighbours(level = level,
                                           Q = best_fit$nb_clusters,
                                           Q_min = Q_min,
                                           Q_max = Q_max,
                                           fit = best_fit,
                                           nb_cores = nb_cores,
                                           init = init)
      new_fit  <-  fits[[which.max(sapply(seq_along(fits),
                                          function(x) fits[[x]]$ICL))]]
      # if (new_fit$ICL < spc_fit$ICL) new_fit <- spc_fit
      if (new_fit$ICL > best_fit$ICL) {
        best_fit  <-  new_fit
        model_list[[best_fit$nb_clusters]] <- best_fit
        ICL[best_fit$nb_clusters] <- best_fit$ICL
        bound[best_fit$nb_clusters] <- best_fit$bound
        print(paste0("# blocks: ", best_fit$nb_clusters,
                     ", ICL = ", best_fit$ICL, " !" ))
      } else {
        condition <- FALSE
      }
    }
    # Descending
    if (is.null(Z)) {
      best_fit <- self$estimate_sbm(level = level,
                                    Q = floor(Q_max/2),
                                    init = init)
      if (is.null(model_list[[best_fit$nb_clusters]]) |
          ICL[best_fit$nb_clusters] < best_fit$ICL) {
        model_list[[floor(Q_max/2)]] <- best_fit
        ICL[best_fit$nb_clusters] <- best_fit$ICL
        bound[best_fit$nb_clusters] <- best_fit$bound
      }
      condition <- TRUE
      print(paste0("# blocks: ", best_fit$nb_clusters,
                   ", ICL = ", best_fit$ICL, " !" ))
      while (condition) {
        fits <- self$estimate_sbm_neighbours(level = level,
                                             Q = best_fit$nb_clusters,
                                             Q_min = Q_min,
                                             Q_max = Q_max,
                                             fit = best_fit,
                                             nb_cores = nb_cores,
                                             init = init)
        new_fit  <-  fits[[
          which.max(sapply(seq_along(fits), function(x) fits[[x]]$ICL))]]
        # # spc_fit  <- self$estimate_sbm(level = level, init = init,
        # #                             Q = max(best_fit$nb_clusters -1, Q_max))
        # if (new_fit$ICL < spc_fit$ICL) new_fit <- spc_fit
        if (new_fit$ICL > best_fit$ICL) {
          best_fit  <-  new_fit
          if (is.null(model_list[[best_fit$nb_clusters]]) |
              ICL[best_fit$nb_clusters] < best_fit$ICL) {
            model_list[[best_fit$nb_clusters]] <- best_fit
            ICL[best_fit$nb_clusters] <- best_fit$ICL
            bound[best_fit$nb_clusters] <- best_fit$bound
          }
          print(paste0("# blocks: ", best_fit$nb_clusters,
                       ", ICL = ", best_fit$ICL, " !" ))
          # } else {
          #   print(paste0("Switching to neighbours mode!"))
          #   new_fit <- self$estimate_sbm_from_neighbours(
          #     Q = best_fit$nb_clusters,
          #     fits = fits)
          #   if ((!is.null(new_fit)) & (new_fit$bound > best_fit$bound)) {
          #     best_fit  <-  new_fit
          #     model_list <- c(model_list, best_fit)
          #     ICL <- c(ICL, best_fit$ICL)
          #     bound <- c(bound, best_fit$bound)
          #     print(paste0("# cluster : ", best_fit$nb_clusters,
          #                  ", ICL = ", best_fit$ICL, " !" ))
        } else {
          condition <- FALSE
          # }
        }
      }
    }
#     models <- lapply(seq(Q_min, Q_max), function(i) NULL)
# #    ICL <- rep(-Inf, Q_max - Q_min +1)
#     Q_list <- sapply(X = 1:length(model_list),
#                      FUN = function(i) model_list[[i]]$nb_clusters)
#     for (i in seq(length(models))) {
#       if(length(bound[Q_list == i]) >0) {
#       models[[i]] <- model_list[[which(bound == max(bound[Q_list == i]))[1]]]
# #        ICL[i] <- models[[i]]$ICL
#       }
#     }
    private$fitted_sbm[[level]] <- model_list
    private$ICLtab_sbm[[level]]    <- ICL
    return(list("models" = model_list, "ICL" = ICL))
  }
)

#-------------------------------------------------------------------------------
# Estimation neighbours for SBM
#-------------------------------------------------------------------------------
MLVSBM$set(
  "public",
  "estimate_sbm_neighbours",
  #' Fit models with size adjacent of a given model
  #'
  #' @param level "lower" or "upper"
  #' @param Q Size of the initial model
  #' @param Q_min Lower bound of exploration space
  #' @param Q_max Upper bound of exploration space
  #' @param fit A FitSBM object from where to explore
  #' @param init Initialization method for additional fits,
  #' one of "hierarchical", "spectral" or "merged_split"
  #' @param nb_cores Number of cores for parallel computing
  #'
  #' @return A list of FitSBM objects
  function(level = "lower",
           Q = NULL, Q_min = 1,
           Q_max = 10,
           fit = NULL,
           nb_cores = NULL,
           init = NULL) {
    os <- Sys.info()["sysname"]
    if (is.null(nb_cores)) {
      if (os != 'Windows') {
        nb_cores <-
          max(parallel::detectCores(all.tests = FALSE, logical = TRUE) %/% 2, 1)
        if (is.na(nb_cores)) nb_cores <- 1
      } else {
        nb_cores <-  1
      }
    }
    models <- list()
    if (Q > Q_min) {
      Z <- merge_clust(fit$Z, fit$nb_clusters)
      fits <- parallel::mclapply(
        X = seq(length(Z)),
        FUN = function(i) {
          self$estimate_sbm(level = level,
                            Q = Q -1,
                            Z = Z[[i]],
                            init = "merge_split")
        }, mc.cores = nb_cores)
      fits  <- c(fits, self$estimate_sbm(level = level,
                                Q = Q-1,
                                init = "spectral"))
      fits <-
        fits[which(sapply( fits, function(x) ! is.null(x)))]
      fits <-
        fits[which(sapply( fits, function(x) ! is.null(x$bound)))]
      fits <-
        fits[[which.max(sapply(seq_along(fits), function(x) fits[[x]]$bound))]]
      models <- c(models, fits)
    }
    if (Q < Q_max) {
      Z <- split_clust(X = fit$adjacency,
                       Z = fit$Z,
                       Q = Q + 1)
      fits <- parallel::mclapply(
        X = seq(length(Z)),
        FUN = function(i) {
          self$estimate_sbm(level = level,
                            Q = fit$nb_clusters + 1,
                            Z = Z[[i]],
                            init = "merge_split")
        }, mc.cores = nb_cores)
      fits  <- c(fits, self$estimate_sbm(level = level,
                                         Q = fit$nb_clusters+1,
                                         init = "spectral"))
      fits <-
        fits[which(sapply( fits, function(x) ! is.null(x)))]
      fits <-
        fits[which(sapply( fits, function(x) ! is.null(x$bound)))]
      fits <-
        fits[[which.max(sapply(seq_along(fits), function(x) fits[[x]]$bound))]]
      models <- c(models, fits)
    }
    return(models)
  })

#-------------------------------------------------------------------------------
# Estimation from neighbours for SBM
#-------------------------------------------------------------------------------
MLVSBM$set(
  "public",
  "estimate_sbm_from_neighbours",
  #' Fit a model of a given size by initiating from its neighbours
  #'
  #' @param level "lower" or "upper"
  #' @param Q Size of the model
  #' @param fits A list of FitSBM object from where to initialize the new model
  #' @param nb_cores Number of cores for parallel computing
  #'
  #' @return A  FitSBM object
  function(level = "lower",
           Q = NULL, fits = NULL,
           nb_cores = NULL) {
    os <- Sys.info()["sysname"]
    if (is.null(nb_cores)) {
      if (os != 'Windows') {
        nb_cores <-
          max(parallel::detectCores(all.tests = FALSE, logical = TRUE) %/% 2, 1)
        if (is.na(nb_cores)) nb_cores <- 1
      } else {
        nb_cores <-  1
      }
    }
    if (is.null(fits)) return(NULL)
    new_fits <- list()
    for (fit in fits) {
      if (Q == fit$nb_clusters - 1) {
        Z <- merge_clust(Z = fit$Z, Q = Q)
        fit_tmp <- parallel::mclapply(
          X = seq(length(Z)),
          FUN = function(i) {
            self$estimate_sbm(level = level,
                              Q = Q,
                              Z = Z[[i]],
                              init = "merge_split")
          }, mc.cores = nb_cores)
        new_fits <- c(new_fits, fit_tmp)
      }
      if (Q == fit$nb_clusters + 1) {
        Z <- split_clust(X = fit$adjacency, Z = fit$Z, Q = Q)
        fit_tmp <- parallel::mclapply(
          X = seq(length(Z)),
          FUN = function(i) {
            self$estimate_sbm(level = level,
                              Q = Q,
                              Z = Z[[i]],
                              init = "merge_split")
          }, mc.cores = nb_cores)
        new_fits <- c(new_fits, fit_tmp)
      }
    }
    new_fits  <- c(new_fits, self$estimate_sbm(level = level,
                                       Q = Q,
                                       init = "spectral"))
    new_fits <-
      new_fits[which(sapply( new_fits, function(x) ! is.null(x)))]
    new_fits <-
      new_fits[which(sapply( new_fits, function(x) ! is.null(x$bound)))]
    new_fit <-
      new_fits[[which.max(sapply(seq_along(fits), function(x) fits[[x]]$bound))]]
    return(new_fit)
  }
  )

#-------------------------------------------------------------------------------
# Estimation for one SBM
#-------------------------------------------------------------------------------
MLVSBM$set(
  "public",
  "estimate_sbm",
  function(level = "lower",
           Q = Q,
           Z = NULL,
           init = "hierarchical") {
    switch(
      EXPR = level,
      lower = {
        fit <- FitSBM$new(Q = Q,
                       X = private$X$I,
                       M = private$M$I,
                       directed = private$directed_$I,
                       distribution = private$distribution_$I)
      },
      upper = {
        fit <- FitSBM$new(Q = Q,
                       X = private$X$O,
                       M = private$M$O,
                       directed = private$directed_$O,
                       distribution = private$distribution_$O)
      },
      stop("Unknown level!!!")
    )
    if(is.null(Z)) {
      fit$do_vem(init = init)
    } else {
      fit$do_vem(init = "merge_split", Z = Z)
    }
    return(fit)
  }
  )
#-------------------------------------------------------------------------------
# Estimation for one MLVSBM
#-------------------------------------------------------------------------------
MLVSBM$set(
  "public",
  "mcestimate",
  function(Q, Z = NULL, init = "hierarchical", independent = FALSE){
  fit <- FitMLVSBM$new(Q = Q,
                    X = private$X,
                    A = private$A,
                    directed = private$directed_,
                    distribution = private$distribution_)
    if (is.null(Z)) {
      fit$do_vem(init = init)
      }  else {
        fit$do_vem(Z = Z, init = "merge_split")
    }
    return(fit)
  })
#-------------------------------------------------------------------------------
# MLVSBM estimate from neighbours
#-------------------------------------------------------------------------------
MLVSBM$set(
  "public",
  "estimate_from_neighbours",
  function(Q,
           models = NULL,
           independent = FALSE,
           nb_cores = nb_cores) {
    if (is.null(models)) return(NULL)
    os <- Sys.info()["sysname"]
    if (is.null(nb_cores)) {
      if (os != 'Windows') {
        nb_cores <-
          max(parallel::detectCores(all.tests = FALSE, logical = TRUE) %/% 2, 1)
        if (is.na(nb_cores)) nb_cores <- 1
      } else {
        nb_cores <-  1
      }
    }
    new_models <- list()
    for (fit in models) {
      if (Q$I == fit$nb_clusters$I + 1) {
        if (Q$O == fit$nb_clusters$O) {
          new_model <- self$mc_ms_estimate(
            Z = list("I" = split_clust(fit$adjacency_matrix$I,
                                       fit$Z$I, fit$nb_clusters$I),
                     "O" = list(fit$Z$O)),
            independent  = independent,
            nb_cores = nb_cores
            )
          new_models <- c(new_models, list(new_model))
        }
      }
      if (Q$I == fit$nb_clusters$I - 1) {
        if (Q$O == fit$nb_clusters$O) {
          new_model <- self$mc_ms_estimate(
            Z = list("I" = merge_clust(Z = fit$Z$I, Q = fit$nb_clusters$O),
                     "O" = list(fit$Z$O)),
            independent = independent,
            nb_cores = nb_cores)
          new_models <- c(new_models, list(new_model))
        }
      }
      if (Q$I == fit$nb_clusters$I) {
        if (Q$O == fit$nb_clusters$O + 1) {
          new_model <- self$mc_ms_estimate(
            Z = list("I" = list(fit$Z$I),
                     "O" = split_clust(fit$adjacency_matrix$O,
                                       fit$Z$O, fit$nb_clusters$O)),
            independent = independent,
            nb_cores = nb_cores)
          new_models <- c(new_models, list(new_model))
        }
      }
      if (Q$I == fit$nb_clusters$I) {
        if (Q$O == fit$nb_clusters$O - 1) {
          new_model <- self$mc_ms_estimate(
            Z = list("I" = list(fit$Z$I),
                     "O" = merge_clust(Z = fit$Z$O, Q = fit$nb_clusters$O)),
            independent = independent,
            nb_cores = nb_cores)
          new_models <- c(new_models, list(new_model))
        }
      }
    }
    new_models <- new_models[which(sapply( new_models,
                                           function(x) ! is.null(x)))]
    new_models <- new_models[which(sapply( new_models,
                                           function(x) ! is.null(x$bound)))]
    new_ICL  <-  which.max(
      sapply(seq_along(new_models), function(x) new_models[[x]]$bound ))
    new_model <- new_models[[new_ICL]]
    return(new_model)
  }
)
#-------------------------------------------------------------------------------
# MLVSBM estimate  neighbours
#-------------------------------------------------------------------------------
MLVSBM$set(
  "public",
  "estimate_neighbours",
  function (Q,
           fit = NULL,
           independent = independent,
           nb_cores = NULL) {
    if (is.null(fit)) return(NULL)
    os <- Sys.info()["sysname"]
    if (is.null(nb_cores)) {
      if (os != 'Windows') {
        nb_cores <-
          max(parallel::detectCores(all.tests = FALSE, logical = TRUE) %/% 2, 1)
        if (is.na(nb_cores)) nb_cores <- 1
      } else {
        nb_cores <-  1
      }
    }
    Z_tmp <-  self$merge_split_membership(fit)
    models <-  list(fit)
    if (! is.null(Z_tmp$I$split[[1]])) {
      if (! is.null(Z_tmp$O$same[[1]])) {
        fitted <- self$mc_ms_estimate(
          Z = list("I" = Z_tmp$I$split,
                   "O" = Z_tmp$O$same),
          independent = independent,
          nb_cores = nb_cores)
        models <-  c(models, list(fitted))
      }
    }
    models <- c(models,
              self$mcestimate(Q = list(I = Q$I +1, O = Q$O), init = "spectral"))
    if (! is.null(Z_tmp$O$split[[1]])) {
      if (! is.null(Z_tmp$I$same[[1]])) {
        fitted <- self$mc_ms_estimate(
          Z = list("I" = Z_tmp$I$same,
                   "O" = Z_tmp$O$split),
          independent = independent,
          nb_cores = nb_cores)
        models <- c(models, list(fitted))
      }
    }
    models <- c(models,
               self$mcestimate(Q = list(I = Q$I, O = Q$O+1), init = "spectral"))
    if (! is.null(Z_tmp$I$merge[[1]])) {
      if (! is.null(Z_tmp$O$same[[1]])) {
        fitted <- self$mc_ms_estimate(
          Z = list("I" = Z_tmp$I$merge,
                   "O" = Z_tmp$O$same),
          independent = independent,
          nb_cores = nb_cores)
        models <- c(models, list(fitted))
      }
    }
    if(Q$I != 1) {
      models <- c(models,
              self$mcestimate(Q = list(I = Q$I-1, O = Q$O), init = "spectral"))
    }
    if (! is.null(Z_tmp$O$merge[[1]])) {
      if (! is.null(Z_tmp$I$same[[1]])) {
        fitted <- self$mc_ms_estimate(
          Z = list("I" = Z_tmp$I$same,
                   "O" = Z_tmp$O$merge),
          independent = independent,
          nb_cores = nb_cores)
        models <- c(models, list(fitted))
      }
    }
    if(Q$O != 1) {
      models <- c(models,
              self$mcestimate(Q = list(I = Q$I, O = Q$O-1), init = "spectral"))
    }
    models <- models[which(sapply( models,function(x) ! is.null(x)))]
    models <- models[which(sapply( models, function(x) ! is.null(x$ICL)))]
    return(models)
  }
)

#-------------------------------------------------------------------------------
# Merge and split the Z for MLVSBM
#-------------------------------------------------------------------------------
MLVSBM$set(
  "public",
  "merge_split_membership",
  function (fitted = private$fitted[[length(private$fitted)]]) {
    estim_Q  <-  list("I" = length(unique(fitted$Z$I)),
                   "O" = length(unique(fitted$Z$O)))
    Z <-  list("I" = list(), "O" = list())
    Z$I$same  <-  list(fitted$Z$I)
    Z$O$same  <-  list(fitted$Z$O)
    if (estim_Q$I <= private$max_Q$I) {
      Z$I$split <- split_clust(fitted$adjacency_matrix$I, fitted$Z$I, estim_Q$I)
      } else {
        Z$I$split <- list(NULL)
        }
    if (estim_Q$O <= private$max_Q$O) {
      Z$O$split <- split_clust(fitted$adjacency_matrix$O, fitted$Z$O, estim_Q$O)
      } else {
        Z$O$split <- list(NULL)
        }
    if (estim_Q$I >= max(2, private$min_Q$I)) {
      Z$I$merge <- merge_clust(fitted$Z$I, estim_Q$I)
      } else {
        Z$I$merge <- list(NULL)
        }
    if (estim_Q$O >= max(2, private$min_Q$O)) {
      Z$O$merge <- merge_clust(fitted$Z$O, estim_Q$O)
      } else {
        Z$O$merge <- list(NULL)
        }
    return(Z)
    }
  )
#-------------------------------------------------------------------------------
# Dispatcher function for MLVSBM
#-------------------------------------------------------------------------------
MLVSBM$set(
  "public",
  "mc_ms_estimate",
  function (Z = NA,
           independent = FALSE,
           nb_cores = NULL) {
    if (is.null(Z$I[[1]]) | is.null(Z$O[[1]])) return(NULL)
    os <- Sys.info()["sysname"]
    if (is.null(nb_cores)) {
      if (os != 'Windows') {
        nb_cores <-
          max(parallel::detectCores(all.tests = FALSE, logical = TRUE) %/% 2, 1)
        if (is.na(nb_cores)) nb_cores <- 1
      } else {
        nb_cores <-  1
      }
    }
    nb_models <- list("I" = length(Z$I), "O" = length(Z$O))
    models  <-  parallel::mclapply(
      seq(nb_models$O * nb_models$I),
      function(x) {
        self$mcestimate(
          Q = list(
            I = max(Z$I[[
              ifelse(! x %% nb_models$I == 0,
                             x %% nb_models$I,
                             nb_models$I)]]) ,
            O = max(Z$O[[(x + nb_models$I - 1) %/% nb_models$I]])),
          Z       = list(
            I = Z$I[[ifelse(! x %% nb_models$I == 0,
                                    x %% nb_models$I,
                                    nb_models$I)]],
            O = Z$O[[(x + nb_models$I - 1) %/% nb_models$I]]),
          init    ="merge_split",
          independent = independent)
      }, mc.cores = nb_cores)
    models <- models[which(sapply( models, function(x) ! is.null(x)))]
    models <- models[which(sapply( models, function(x) ! is.null(x$bound)))]
    private$tmp_fitted <- c(private$tmp_fitted, models)
    best_model <- models[[which.max(
      sapply(seq_along(models), function(x) {models[[x]]$bound}))]]
    return(best_model)
    }
)
#-------------------------------------------------------------------------------
# Estimate with known numbers of clusters
#-------------------------------------------------------------------------------
MLVSBM$set(
  "public",
  "estimate_one",
#' Infer a MLVSBM for a given model size
#'
#' @param Q A list of integers, the model size
#' @param Z A list of vectors, the initial clustering (default to NULL)
#' @param independent A boolean, shall the levels be inferred independently
#'
#' @return A FitSBM object of the given model size
#' @export
  function(Q,
           Z = NULL,
           independent = FALSE,
           init = "hierarchical",
           nb_cores = NULL) {
    os <- Sys.info()["sysname"]
    if (is.null(nb_cores)) {
      if (os != 'Windows') {
        nb_cores <-
          max(parallel::detectCores(all.tests = FALSE, logical = TRUE) %/% 2, 1)
        if (is.na(nb_cores)) nb_cores <- 1
      } else {
        nb_cores <-  1
      }
    }
    if (is.null(Z)) {
      fit_tmp <- self$mcestimate(Q = Q,
                                 independent = independent)
    } else {
      fit_tmp <- self$mcestimate(Q = Q,
                                  Z = Z,
                                  init = init,
                                  independent = independent)
    }
    Z_tmp   <-  self$merge_split_membership(fit_tmp)
    models  <-  list(fit_tmp)
    if (! is.null(Z_tmp$I$split[[1]])) {
      if (! is.null(Z_tmp$O$same[[1]])) {
        fitted <-  self$mc_ms_estimate(
          Z = list("I" = Z_tmp$I$split,
                   "O" = Z_tmp$O$same),
          independent = independent,
          nb_cores = nb_cores)
        best   <-  self$mc_ms_estimate(
          Z = list("I" = merge_clust(fitted$Z$I, fitted$nb_clusters$I),
                   "O" = list(fitted$Z$O)),
          independent = independent,
          nb_cores = nb_cores)
        print(paste0("====== Searching neighbours...",
                     " # Individual clusters : ", best$nb_clusters$I,
                     " , # Organisation clusters ", best$nb_clusters$O,
                     ",  ICL : ", best$ICL, "========"))
        models <-  c(models, list(best))
      }
    }
    if (! is.null(Z_tmp$O$split[[1]])) {
      if (! is.null(Z_tmp$I$same[[1]])) {
        fitted <-  self$mc_ms_estimate(
          Z = list("I" = Z_tmp$I$same,
                   "O" = Z_tmp$O$split),
          independent = independent,
          nb_cores = nb_cores)
        best   <-  self$mc_ms_estimate(
          Z = list("I" = list(fitted$Z$I),
                   "O" = merge_clust(fitted$Z$O, fitted$nb_clusters$O)),
          independent = independent,
          nb_cores = nb_cores)
        models <-  c(models, list(best))
        print(paste0("====== Searching neighbours...",
                     " # Individual clusters : ", best$nb_clusters$I,
                     " , # Organisation clusters ", best$nb_clusters$O,
                     ",  ICL : ", best$ICL, "========"))
      }
    }
    if (! is.null(Z_tmp$I$merge[[1]])) {
      if (! is.null(Z_tmp$O$same[[1]])) {
        fitted <- self$mc_ms_estimate(
          Z = list("I" = Z_tmp$I$merge,
                   "O" = Z_tmp$O$same),
          independent = independent,
          nb_cores = nb_cores)
        best <- self$mc_ms_estimate(
          Z = list("I" = split_clust(fitted$adjacency_matrix$I,
                                     fitted$Z$I, fitted$nb_clusters$I),
                   "O"= list(fitted$Z$O)),
          independent = independent,
          nb_cores = nb_cores)
        models <- c(models, list(best))
        print(paste0("====== Searching neighbours...",
                     " # Individual clusters : ", best$nb_clusters$I,
                     " , # Organisation clusters ", best$nb_clusters$O,
                     ",  ICL : ", best$ICL, "========"))
      }
    }
    if (! is.null(Z_tmp$O$merge[[1]])) {
      if (! is.null(Z_tmp$I$same[[1]])) {
        fitted <- self$mc_ms_estimate(
          Z = list("I" = Z_tmp$I$same,
                   "O" = Z_tmp$O$merge),
          independent = independent,
          nb_cores = nb_cores)
        best  <- self$mc_ms_estimate(
          Z = list("I" = list(fitted$Z$I),
                   "O" = split_clust(fitted$adjacency_matrix$O,
                                     fitted$Z$O, fitted$nb_clusters$O)),
          independent = independent,
          nb_cores = nb_cores)
        models <- c(models, list(best))
        print(paste0("====== Searching neighbours...",
                     " # Individual clusters : ", best$nb_clusters$I,
                     " , # Organisation clusters ", best$nb_clusters$O,
                     ",  ICL : ", best$ICL, "========"))
      }
    }
    models <- models[which(sapply( models, function(x) ! is.null(x)))]
    models <- models[which(sapply( models, function(x) ! is.null(x$bound)))]
    best_ICL <- which.max(
      sapply(seq_along(models), function(x) {models[[x]]$bound}))
    best_model <- models[[best_ICL]]
    print(paste0("======= # Individual clusters : ", best_model$nb_clusters$I,
                 " , # Organisation clusters ", best_model$nb_clusters$O,
                 ",  ICL : ", best_model$ICL, "========"))
    self$addmodel(best_model)
    return(best_model)
  }
)

#-------------------------------------------------------------------------------
# estimate with initialization
#-------------------------------------------------------------------------------
MLVSBM$set(
  "public",
  "estimate_all_bm",
#' Infer a MLVSBM with a greedy algorithm to navigate between different size
#' of models
#'
#' @param Q A list of integers, the initial model size, default to NULL
#' @param Z A list of vectors, the initial clustering, default to NULL
#' @param independent A boolean, are the
#' @param clear A boolean, should all previous models list be cleared from the
#' object
#'
#' @return The FitSBM object with the best ICL
#' @export
  function (Q = NULL, Z = NULL, independent = FALSE,
            clear = TRUE, nb_cores = NULL) {
    os <- Sys.info()["sysname"]
    if (is.null(nb_cores)) {
      if (os != 'Windows') {
        nb_cores <-
          max(parallel::detectCores(all.tests = FALSE, logical = TRUE) %/% 2, 1)
        if (is.na(nb_cores)) nb_cores <- 1
      } else {
        nb_cores <-  1
      }
    }
    if (clear) self$clearmodels()
    best_model  <-  self$mcestimate(Q = Q, Z = Z, independent = independent)
    print(paste0("======= # Individual clusters : ", best_model$nb_clusters$I,
                 " , # Organisation clusters ", best_model$nb_clusters$O,
                 ",  ICL : ", best_model$ICL, "========"))
    self$addmodel(best_model)
    condition <- TRUE
    while (condition) {
      models <- self$estimate_neighbours(
        fit = best_model,
        Q = list("I" = best_model$nb_clusters$I,
                 "O" = best_model$nb_clusters$O),
        nb_cores = nb_cores
        )
# Z   = self$merge_split_membership(best_model)
# models = parallel::mcmapply(rep(1:3, times = 3), rep(1:3, each = 3),
#                             FUN = function(i, j) {
#                               return(self$mc_ms_estimate(Z$I[[i]], Z$O[[j]], model = model))
#                             },
#                             mc.cores = nb_cores)
# models <-  models[which(sapply( models, x <- function(x) ! is.null(x)))]
# models <-  models[which(sapply( models, x <- function(x) ! is.null(x$ICL)))]
#      private$tmp_fitted = c(private$tmp_fitted, models)
      new_best_ICL <-  which.max(
        sapply(seq_along(models), function(x) models[[x]]$ICL))
      new_best_model <- models[[new_best_ICL]]
      # new_best_model <-  models[[which.max(sapply(1:length(models), x <- function(x) {models[[x]]$ICL}))]]
      if (new_best_model$ICL > best_model$ICL) {
        best_model  <-  new_best_model
        self$addmodel(new_best_model)
      # } else {
      #   print(paste0("Switching to neighbours mode!"))
      #   new_model <- self$estimate_from_neighbours(
      #     Q = best_model$nRClusters,
      #     S = best_model$nLClusters,
      #     models = models,
      #     type = type)
      #   if ((!is.null(new_model)) & (new_model$bound > best_model$bound)) {
      #     best_model  <-  new_model
      #     self$addmodel(new_model)
                } else{
                  condition <- FALSE
                  # }
                }
      print(paste0("======= # Individual blocks : ", best_model$nb_clusters$I,
                   " , # Organizational blocks ", best_model$nb_clusters$O,
                   ",  ICL : ", best_model$ICL, "========"))
    }
  return(best_model)
  }
)
