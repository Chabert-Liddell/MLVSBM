
<!-- README.md is generated from README.Rmd. Please edit that file -->
MLVSBM
======

<!-- badges: start -->
[![Travis build status](https://travis-ci.org/Chabert-Liddell/MLVSBM.svg?branch=master)](https://travis-ci.org/Chabert-Liddell/MLVSBM) [![Codecov test coverage](https://codecov.io/gh/Chabert-Liddell/MLVSBM/branch/master/graph/badge.svg)](https://codecov.io/gh/Chabert-Liddell/MLVSBM?branch=master) <!-- badges: end -->

The goal of MLVSBM is to ...

Installation
------------

You can install the released version of MLVSBM from [github](https://github.com) with:

``` r
devtools::install_github("Chabert-Liddell/MLVSBM")
```

Example
-------

This is a how to simulate a multilevel SBM network:

``` r
set.seed(1)
my_mlvsbm <- MLVSBM::mlvsbm_simulate_network(
  n = list(100, 50),
  Q = list(2, 2),
  pi = c(.5, .5),
  gamma = matrix(c(.8, .2, .1, .9), ncol = 2),
  alpha = list(matrix(c(.25, .1, .1, .25), nrow = 2, ncol = 2),
               matrix(c(.7, .4, .4, .1), nrow = 2, ncol = 2)),
  directed = list(FALSE, FALSE),
  affiliation = "preferential",
  no_empty_org = FALSE)
```

This is how to create a network from data:

``` r
lower_level <- my_mlvsbm$adjacency_matrix$I
upper_level <- my_mlvsbm$adjacency_matrix$O
affiliation <- my_mlvsbm$affiliation_matrix
my_mlvsbm2 <- MLVSBM::mlvsbm_create_network(X = list(lower_level, upper_level),
                                            A = affiliation)
```

And this is how to infer it:

``` r
fit <- MLVSBM:::mlvsbm_estimate_network(my_mlvsbm)
#> [1] "Infering lower level :"
#> [1] "# cluster : 1, ICL = -2372.16770052719 !"
#> 1  :  -2372.118 
1  :  -2333.463 
[1] "# cluster : 5, ICL = -2490.06594742428 !"
#> [1] "# cluster : 4, ICL = -2401.94555812121 !"
#> [1] "# cluster : 3, ICL = -2373.60113971979 !"
#> [1] "# cluster : 2, ICL = -2355.91131355913 !"
#> [1] "Infering upper level :"
#> [1] "# cluster : 1, ICL = -840.403791989084 !"
#> 1  :  -773.1229 
[1] "# cluster : 2, ICL = -796.208742941495 !"
#> 1  :  -760.2382 
[1] "# cluster : 5, ICL = -847.872907916854 !"
#> [1] "# cluster : 4, ICL = -811.946225345139 !"
#> [1] "# cluster : 3, ICL = -791.749205113861 !"
#> [1] "# cluster : 2, ICL = -775.438351381624 !"
#> 1  :  -3070.956 
2  :  -3070.93 
3  :  -3070.929 
4  :  -3070.929 
5  :  -3070.929 
6  :  -3070.929 
7  :  -3070.929 
8  :  -3070.929 
[1] "======= # Individual clusters : 2 , # Organisation clusters 2,  ICL : -3100.91668983441========"
#> 1  :  -3127.319 
2  :  -3127.302 
3  :  -3127.3 
4  :  -3127.3 
5  :  -3127.3 
6  :  -3127.3 
7  :  -3127.3 
8  :  -3127.3 
9  :  -3127.3 
10  :  -3127.3 
1  :  -3169.149 
2  :  -3169.077 
3  :  -3169.066 
4  :  -3169.064 
5  :  -3169.064 
6  :  -3169.064 
7  :  -3169.064 
8  :  -3169.064 
9  :  -3169.064 
10  :  -3169.064 
11  :  -3169.064 
12  :  -3169.064 
[1] "======= # R clusters : 2 , # L clusters 2,  ICL : -3100.91668983441========"
#> [1] "ICL for independent levels : -3131.34966494076"
#> [1] "ICL for interdependent levels : -3100.91668983441"
#> [1] "=====Interdepence is detected between the two level====="
```

<img src="man/figures/README-pressure-1.png" width="100%" />
