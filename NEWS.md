# MLVSBM (development version)

* Changes to `mlvsbm_estimate_generalized_network()`:
  + Fixed bug when setting the number of clusters.
  + Added init_points and Qmax to the fit_options argument.
  + The default algorithm is now slower but handles better underfitting 
(fit_option$init_points == "all"). Users may fix it to any other value to revert 
to the old behavior.
  + When fixing the number of clusters, the algorithm will better explore 
  neighbors models.
  + Returned blocks are reordered so that the affiliation parameters (gamma) is 
  as "diagonal" as possible.

# MLVSBM 0.3.0
* Added support for generalized multilevel networks:
* This includes three new user functions: `mlvsbm_simulate_generalized_network()`,
`mlvsbm_create_generalized_network()` and `mlvsbm_estimate_generalized_network()`.
* This function added the following features: "poisson" distribution of the 
interactions to deal with count data, support for network of more than two levels 
(temporal networks e.g.), support for networks with nodes with no affiliation.

# MLVSBM 0.2.4 (latest cran release)
* Minor bug fix

# MLVSBM 0.2.3
* Changed `plot` function and method of FitMLVSBM object to graphon view type.

# MLVSBM 0.2.2
* Improved initialization for networks with missing data.
* Fixed a bug for directed networks.

# MLVSBM 0.2.1

# MLVSBM 0.2.0
* Added S3 generic function `print`, `plot`, `predict` and 
`coef` for FitMLVSBM object.
* R6 class methods are better documented.

# MLVSBM 0.1.4
* Reverse the `sbm` dependency to `blockmodels` for future `sbm` compatibility.

# MLVSBM 0.1.3
* Fixed a bug with `NA` interactions on user functions.
* added `sbm` package dependency for better initialization.

# MLVSBM 0.1.2
* Added options for `mlvl_estimate_network()` function.
* It is now possible to simulate networks with no isolated nodes.

# MLVSBM 0.1.1
* Fixed bug with `parallel` for Windows.

# MLVSBM 0.1
* Finished user functions.
* Tutorial.

# MLVSBM 0.0.0.9000

* Added a `NEWS.md` file to track changes to the package.
