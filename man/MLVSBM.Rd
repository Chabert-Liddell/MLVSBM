% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MLVSBM-class.R
\name{MLVSBM}
\alias{MLVSBM}
\title{R6Class for multilevel object}
\description{
Store all simulation parameters and list of fittedmodels.
Methods for global inference and model selection are included.
}
\section{Active bindings}{
\if{html}{\out{<div class="r6-active-bindings">}}
\describe{
\item{\code{nb_nodes}}{List of the umber of nodes for each levels}

\item{\code{simulation_parameters}}{List of parameters of the MLVSBM}

\item{\code{affiliation_matrix}}{Access the affiliation matrix}

\item{\code{adjacency_matrix}}{Access the list of adjacency_matrix}

\item{\code{memberships}}{Access the list of the clusterings}

\item{\code{fittedmodels}}{Get the list of selected fitted FitMLVSBM objects}

\item{\code{ICL}}{A summary table of selected fitted models and ICL model selection criterion}

\item{\code{ICL_sbm}}{Summary table of ICL by levels}

\item{\code{tmp_fittedmodels}}{A list of all fitted FitMLVSBM objects}

\item{\code{fittedmodels_sbm}}{A list of selected fitted FitSBM objects of each levels}

\item{\code{max_clusters}}{Access the list of maximum model size}

\item{\code{min_clusters}}{Access the list of minimum model size}

\item{\code{directed}}{Access the list of boolean for levels  direction}

\item{\code{directed}}{Access the list of the distribution used for each levels}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-MLVSBM-estimate_level}{\code{MLVSBM$estimate_level()}}
\item \href{#method-MLVSBM-estimate_sbm_neighbours}{\code{MLVSBM$estimate_sbm_neighbours()}}
\item \href{#method-MLVSBM-estimate_sbm_from_neighbours}{\code{MLVSBM$estimate_sbm_from_neighbours()}}
\item \href{#method-MLVSBM-estimate_sbm}{\code{MLVSBM$estimate_sbm()}}
\item \href{#method-MLVSBM-mcestimate}{\code{MLVSBM$mcestimate()}}
\item \href{#method-MLVSBM-estimate_from_neighbours}{\code{MLVSBM$estimate_from_neighbours()}}
\item \href{#method-MLVSBM-estimate_neighbours}{\code{MLVSBM$estimate_neighbours()}}
\item \href{#method-MLVSBM-merge_split_membership}{\code{MLVSBM$merge_split_membership()}}
\item \href{#method-MLVSBM-mc_ms_estimate}{\code{MLVSBM$mc_ms_estimate()}}
\item \href{#method-MLVSBM-estimate_one}{\code{MLVSBM$estimate_one()}}
\item \href{#method-MLVSBM-estimate_all_bm}{\code{MLVSBM$estimate_all_bm()}}
\item \href{#method-MLVSBM-new}{\code{MLVSBM$new()}}
\item \href{#method-MLVSBM-findmodel}{\code{MLVSBM$findmodel()}}
\item \href{#method-MLVSBM-clearmodels}{\code{MLVSBM$clearmodels()}}
\item \href{#method-MLVSBM-addmodel}{\code{MLVSBM$addmodel()}}
\item \href{#method-MLVSBM-simulate}{\code{MLVSBM$simulate()}}
\item \href{#method-MLVSBM-clone}{\code{MLVSBM$clone()}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-MLVSBM-estimate_level"></a>}}
\if{latex}{\out{\hypertarget{method-MLVSBM-estimate_level}{}}}
\subsection{Method \code{estimate_level()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{MLVSBM$estimate_level(
  level = "lower",
  Q_min = 1,
  Q_max = 10,
  Z = NULL,
  init = "hierarchical",
  depth = 1,
  nb_cores = NULL
)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-MLVSBM-estimate_sbm_neighbours"></a>}}
\if{latex}{\out{\hypertarget{method-MLVSBM-estimate_sbm_neighbours}{}}}
\subsection{Method \code{estimate_sbm_neighbours()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{MLVSBM$estimate_sbm_neighbours(
  level = "lower",
  Q = NULL,
  Q_min = 1,
  Q_max = 10,
  fit = NULL,
  nb_cores = NULL,
  init = NULL
)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-MLVSBM-estimate_sbm_from_neighbours"></a>}}
\if{latex}{\out{\hypertarget{method-MLVSBM-estimate_sbm_from_neighbours}{}}}
\subsection{Method \code{estimate_sbm_from_neighbours()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{MLVSBM$estimate_sbm_from_neighbours(
  level = "lower",
  Q = NULL,
  fits = NULL,
  nb_cores = NULL
)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-MLVSBM-estimate_sbm"></a>}}
\if{latex}{\out{\hypertarget{method-MLVSBM-estimate_sbm}{}}}
\subsection{Method \code{estimate_sbm()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{MLVSBM$estimate_sbm(level = "lower", Q = Q, Z = NULL, init = "hierarchical")}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-MLVSBM-mcestimate"></a>}}
\if{latex}{\out{\hypertarget{method-MLVSBM-mcestimate}{}}}
\subsection{Method \code{mcestimate()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{MLVSBM$mcestimate(Q, Z = NULL, init = "hierarchical", independent = FALSE)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-MLVSBM-estimate_from_neighbours"></a>}}
\if{latex}{\out{\hypertarget{method-MLVSBM-estimate_from_neighbours}{}}}
\subsection{Method \code{estimate_from_neighbours()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{MLVSBM$estimate_from_neighbours(
  Q,
  models = NULL,
  independent = FALSE,
  nb_cores = nb_cores
)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-MLVSBM-estimate_neighbours"></a>}}
\if{latex}{\out{\hypertarget{method-MLVSBM-estimate_neighbours}{}}}
\subsection{Method \code{estimate_neighbours()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{MLVSBM$estimate_neighbours(
  Q,
  fit = NULL,
  independent = independent,
  nb_cores = NULL
)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-MLVSBM-merge_split_membership"></a>}}
\if{latex}{\out{\hypertarget{method-MLVSBM-merge_split_membership}{}}}
\subsection{Method \code{merge_split_membership()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{MLVSBM$merge_split_membership(
  fitted = private$fitted[[length(private$fitted)]]
)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-MLVSBM-mc_ms_estimate"></a>}}
\if{latex}{\out{\hypertarget{method-MLVSBM-mc_ms_estimate}{}}}
\subsection{Method \code{mc_ms_estimate()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{MLVSBM$mc_ms_estimate(Z = NA, independent = FALSE, nb_cores = NULL)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-MLVSBM-estimate_one"></a>}}
\if{latex}{\out{\hypertarget{method-MLVSBM-estimate_one}{}}}
\subsection{Method \code{estimate_one()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{MLVSBM$estimate_one(
  Q,
  Z = NULL,
  independent = FALSE,
  init = "hierarchical",
  nb_cores = NULL
)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-MLVSBM-estimate_all_bm"></a>}}
\if{latex}{\out{\hypertarget{method-MLVSBM-estimate_all_bm}{}}}
\subsection{Method \code{estimate_all_bm()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{MLVSBM$estimate_all_bm(
  Q = NULL,
  Z = NULL,
  independent = FALSE,
  clear = TRUE,
  nb_cores = NULL
)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-MLVSBM-new"></a>}}
\if{latex}{\out{\hypertarget{method-MLVSBM-new}{}}}
\subsection{Method \code{new()}}{
Constructor for R6 class MLVSBM
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{MLVSBM$new(
  n = NULL,
  X = NULL,
  A = NULL,
  Z = NULL,
  directed = NULL,
  sim_param = NULL,
  distribution = list("bernoulli", "bernoulli")
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{n}}{A list of size 2, the number of nodes}

\item{\code{X}}{A list of 2 adjacency matrices}

\item{\code{A}}{The affiliation matrix}

\item{\code{Z}}{A list of 2 vectors, the blocks membership}

\item{\code{directed}}{A list of 2 booleans}

\item{\code{sim_param}}{A list of MLVSBM parameters for simulating networks}

\item{\code{distribution}}{The distributions of the interactions ("bernoulli")}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
A MLVSBM object
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-MLVSBM-findmodel"></a>}}
\if{latex}{\out{\hypertarget{method-MLVSBM-findmodel}{}}}
\subsection{Method \code{findmodel()}}{
Find a fitted model of a given size
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{MLVSBM$findmodel(nb_clusters = NA, fit = NA)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{nb_clusters}}{A list of the size of the model}

\item{\code{fit}}{if fit = "best" return the best model according to the ICL}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
A FitMLVSBM object
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-MLVSBM-clearmodels"></a>}}
\if{latex}{\out{\hypertarget{method-MLVSBM-clearmodels}{}}}
\subsection{Method \code{clearmodels()}}{
delete all fitted models
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{MLVSBM$clearmodels()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-MLVSBM-addmodel"></a>}}
\if{latex}{\out{\hypertarget{method-MLVSBM-addmodel}{}}}
\subsection{Method \code{addmodel()}}{
Added a FitMLVSBM object to the list of fitted model
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{MLVSBM$addmodel(fit)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{fit}}{The FitMLVSBM object to be added}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-MLVSBM-simulate"></a>}}
\if{latex}{\out{\hypertarget{method-MLVSBM-simulate}{}}}
\subsection{Method \code{simulate()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{MLVSBM$simulate()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-MLVSBM-clone"></a>}}
\if{latex}{\out{\hypertarget{method-MLVSBM-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{MLVSBM$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}
