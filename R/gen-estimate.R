#-------------------------------------------------------------------------------
#  SBM inference for upper or lower level
#-------------------------------------------------------------------------------
GenMLVSBM$set(
  "public",
  "estimate_level",
  #' Fit a SBM on a given level of a GenMLVSBM object
  #'
  #' @param level an Integer from 1 (upper level) to L (lower level)
  #' @param Q_min Lower bound of exploration space
  #' @param Q_max Upper bound of exploration space
  #' @param Z Initial clustering if any
  #' @param init Initialization method one of "hierarchical",
  #' "spectral" or "merged_split"
  #' @param depth Exploration depth (not used)
  #' @param nb_cores Number of cores for parallel computing
  #'
  #' @return A list of FitSBM objects
  function(level = NULL,
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
    print(paste0("Infering level ", level, ":"))
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
        } else {
          condition <- FALSE
          # }
        }
      }
    }
    private$fitted_sbm[[level]] <- model_list
    private$ICLtab_sbm[[level]]    <- ICL
    return(list("models" = model_list, "ICL" = ICL))
  }
)

#-------------------------------------------------------------------------------
# Estimation neighbours for SBM
#-------------------------------------------------------------------------------
GenMLVSBM$set(
  "public",
  "estimate_sbm_neighbours",
  #' Fit models with size adjacent of a given model
  #'
  #' @param level an Integer from 1 (upper level) to L (lower level)
  #' @param Q Size of the initial model
  #' @param Q_min Lower bound of exploration space
  #' @param Q_max Upper bound of exploration space
  #' @param fit A FitSBM object from where to explore
  #' @param init Initialization method for additional fits,
  #' one of "hierarchical", "spectral" or "merged_split"
  #' @param nb_cores Number of cores for parallel computing
  #'
  #' @return A list of FitSBM objects
  function(level = NULL,
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
GenMLVSBM$set(
  "public",
  "estimate_sbm_from_neighbours",
  #' Fit a model of a given size by initiating from its neighbours
  #'
  #' @param level an Integer from 1 (upper level) to L (lower level)
  #' @param Q Size of the model
  #' @param fits A list of FitSBM object from where to initialize the new model
  #' @param nb_cores Number of cores for parallel computing
  #'
  #' @return A  FitSBM object
  function(level = NULL,
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
        Z <- split_clust(X = fit$adjacency_matrix, Z = fit$Z, Q = Q)
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
GenMLVSBM$set(
  "public",
  "estimate_sbm",
  function(level = NULL,
           Q = Q,
           Z = NULL,
           init = "hierarchical") {
    fit <- FitSBM$new(Q = Q,
                      X = private$X[[level]],
                      M = private$M[[level]],
                      directed = private$directed_[level],
                      distribution = private$distribution_[level])
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
GenMLVSBM$set(
  "public",
  "mcestimate",
  function(Q, Z = NULL, init = "hierarchical", independent = FALSE){
    fit <- FitGenMLVSBM$new(Q = Q,
                         X = private$X,
                         A = private$A,
                         M = private$M,
                         directed = private$directed_,
                         distribution = private$distribution_)
    fit$fit_options <- self$fit_options
    if (is.null(Z)) {
      fit$do_vem(init = init)
    }  else {
      fit$do_vem(Z = Z, init = "merge_split")
    }
    return(fit)
  })

#-------------------------------------------------------------------------------
# MLVSBM estimate  neighbours
#-------------------------------------------------------------------------------
GenMLVSBM$set(
  "public",
  "estimate_neighbours",
  function (level,
            fit = NULL,
            Q,
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
    models <-  list(fit)
    if (Q[level] <= private$max_Q[level]){
      Z_split <- split_clust(fit$adjacency_matrix[[level]], fit$Z[[level]],
                             Q[level])
      Q_plus <- Q
      Q_plus[level] <- Q[level] +1
      models  <-  c(
        models,
        parallel::mclapply(
          Z_split,
          function(z) {
            Z <- fit$Z
            Z[[level]] <- z
            self$mcestimate(
              Q = Q_plus,
              Z = Z,
              init    ="merge_split",
              independent = independent)
          }, mc.cores = nb_cores)
      )
    }
    if (Q[level] >= max(2, private$min_Q[level])){
      Q_minus <- Q
      Z_merge <- merge_clust(fit$Z[[level]], Q[level])
      Q_minus[level] <- Q[level] - 1
      models  <-  c(
        models,
        parallel::mclapply(
          Z_merge,
          function(z) {
            Z <- fit$Z
            Z[[level]] <- z
            self$mcestimate(
              Q = Q_minus,
              Z = Z,
              init    ="merge_split",
              independent = independent)
          }, mc.cores = nb_cores)
      )
    }
    models <- models[which(sapply( models,function(x) ! is.null(x)))]
    models <- models[which(sapply( models, function(x) ! is.null(x$ICL)))]
    return(models)
  }
)
#-------------------------------------------------------------------------------
# Estimate with known numbers of clusters
#-------------------------------------------------------------------------------
GenMLVSBM$set(
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
    best_model <- fit_tmp
    print(paste0("====== Searching neighbours...========"))
    for (m in seq(private$L)) {
      models <- self$estimate_neighbours(
        level = m,
        fit = best_model,
        Q = best_model$nb_clusters,
        nb_cores = nb_cores
      )
      new_best_ICL <-  which.max(
        sapply(seq_along(models), function(x) models[[x]]$ICL))
      new_best_model <- models[[new_best_ICL]]
      if (new_best_model$ICL > best_model$ICL) {
        best_model  <-  new_best_model
        print(paste0("======= # Blocks : ", toString(best_model$nb_clusters),
                     ",  ICL : ", best_model$ICL, "========"))
      }
    }
    if (any(best_model$Q != Q)) {
      for (m in seq(private$L-1, 1)) {
        if(best_model$Q[m] != Q[m]) {
          models <- self$estimate_neighbours(
            level = m,
            fit = best_model,
            Q = best_model$nb_clusters,
            nb_cores = nb_cores
          )
          new_best_ICL <-  which.max(
            sapply(seq_along(models), function(x) models[[x]]$ICL))
          new_best_model <- models[[new_best_ICL]]
          if (new_best_model$ICL > best_model$ICL) {
            best_model  <-  new_best_model
            print(paste0("======= # Blocks : ", toString(best_model),
                         ",  ICL : ", best_model$ICL, "========"))
          }
        }
      }
    }
    print(paste0("====== Back to desired model size...======="))
    for (m in seq(private$L)) {
      if (best_model$nb_clusters[m] == Q[m] + 1) {
        Q_minus <- best_model$nb_clusters
        Z_merge <- merge_clust(best_model$Z[[m]], Q[m] + 1)
        Q_minus[m] <- Q_minus[m] - 1
        models  <-
          parallel::mclapply(
            Z_merge,
            function(z) {
              Z <- fit$Z
              Z[[level]] <- z
              self$mcestimate(
                Q = Q_minus,
                Z = Z,
                init    ="merge_split",
                independent = independent)
            }, mc.cores = nb_cores)
        new_best_ICL <-  which.max(
          sapply(seq_along(models), function(x) models[[x]]$ICL))
        best_model <- models[[new_best_ICL]]
      }
      if (best_model$nb_clusters[m] == Q[m] - 1) {
        Z_split <- split_clust(best_model$adjacency_matrix[[m]],
                               best_model$Z[[m]],
                               Q[m])
        Q_plus <- best_model$nb_clusters
        Q_plus[m] <- Q_plus[m] +1
        models  <-
          parallel::mclapply(
            Z_split,
            function(z) {
              Z <- fit$Z
              Z[[m]] <- z
              self$mcestimate(
                Q = Q_plus,
                Z = Z,
                init    ="merge_split",
                independent = independent)
            }, mc.cores = nb_cores)
        new_best_ICL <-  which.max(
          sapply(seq_along(models), function(x) models[[x]]$ICL))
        best_model <- models[[new_best_ICL]]
      }
    }
    models <- list(fit_tmp, best_model)
    models <- models[which(sapply( models, function(x) ! is.null(x)))]
    models <- models[which(sapply( models, function(x) ! is.null(x$bound)))]
    best_ICL <- which.max(
      sapply(seq_along(models), function(x) {models[[x]]$bound}))
    best_model <- models[[best_ICL]]
    print(paste0("======= # Blocks : ", toString(best_model$nb_clusters),
                 ",  ICL : ", best_model$ICL, "========"))
    self$addmodel(best_model)
    return(best_model)
  }
)

#-------------------------------------------------------------------------------
# estimate with initialization
#-------------------------------------------------------------------------------
GenMLVSBM$set(
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
    #browser()
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
    print(paste0("======= # Blocks : ", toString(best_model$nb_clusters),
                 ",  ICL : ", best_model$ICL, "========"))
    self$addmodel(best_model)
    condition <- TRUE
    while (condition) {
      improved <- FALSE
      for (m in seq(private$L)) {
        models <- self$estimate_neighbours(
          level = m,
          fit = best_model,
          Q = best_model$nb_clusters,
          nb_cores = nb_cores
        )
        new_best_ICL <-  which.max(
          sapply(seq_along(models), function(x) models[[x]]$ICL))
        new_best_model <- models[[new_best_ICL]]
        if (new_best_model$ICL > best_model$ICL) {
          best_model  <-  new_best_model
          self$addmodel(new_best_model)
          print(paste0("======= # Blocks : ", toString(best_model$nb_clusters),
                       ",  ICL : ", best_model$ICL, "========"))
          improved <- TRUE
        }
      }
      condition <- improved
    }
    return(best_model)
  }
)
