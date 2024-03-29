% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FitGenMLVSBM-class.R
\name{FitGenMLVSBM}
\alias{FitGenMLVSBM}
\title{An R6 Class object, a fitted generalized multilevel network once $dovem() is done}
\description{
An R6 Class object, a fitted generalized multilevel network once $dovem() is done

An R6 Class object, a fitted generalized multilevel network once $dovem() is done
}
\section{Public fields}{
\if{html}{\out{<div class="r6-fields">}}
\describe{
\item{\code{vbound}}{The vector of variational bound for monitoring convergence}
}
\if{html}{\out{</div>}}
}
\section{Active bindings}{
\if{html}{\out{<div class="r6-active-bindings">}}
\describe{
\item{\code{affiliation_matrix}}{Get the affiliation matrix}

\item{\code{adjacency_matrix}}{Get the list of adjacency matrices}

\item{\code{nb_nodes}}{Get the list of the number of nodes}

\item{\code{nb_clusters}}{Get the list of the number of blocks}

\item{\code{nb_levels}}{Get the number of levels}

\item{\code{block_proportions}}{Get the block proportions of each level}

\item{\code{parameters}}{Get the list of the model parameters}

\item{\code{membership}}{Get the list of the variational parameters}

\item{\code{independent}}{Are the levels independent?}

\item{\code{distribution}}{Emission distribution of each level}

\item{\code{directed}}{Are the levels directed?}

\item{\code{entropy}}{Get the entropy of the model}

\item{\code{bound}}{Get the variational bound of the model}

\item{\code{df_mixture}}{Get the degrees of freedom of the
mixture parameters}

\item{\code{df_connect}}{Get the degrees of freedom of the
connection parameters}

\item{\code{connect}}{Get the number of possible observed connections}

\item{\code{ICL}}{Get the ICL model selection criterion of the model}

\item{\code{full_penalty}}{Get the penalty used to compute the ICL}

\item{\code{Z}}{Get the list of block memberships (vector form)}

\item{\code{X_hat}}{Get the list of the matrices of probability connection
predictions}

\item{\code{map}}{Get the list of block memberships (matrix form)}

\item{\code{penalty}}{Get the ICL penalty}

\item{\code{likelihood}}{Compute the likelihood of both levels}

\item{\code{complete_likelihood}}{Get the complete likelihood of the model}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-FitGenMLVSBM-new}{\code{FitGenMLVSBM$new()}}
\item \href{#method-FitGenMLVSBM-update_alpha}{\code{FitGenMLVSBM$update_alpha()}}
\item \href{#method-FitGenMLVSBM-update_pi}{\code{FitGenMLVSBM$update_pi()}}
\item \href{#method-FitGenMLVSBM-update_gamma}{\code{FitGenMLVSBM$update_gamma()}}
\item \href{#method-FitGenMLVSBM-init_clustering}{\code{FitGenMLVSBM$init_clustering()}}
\item \href{#method-FitGenMLVSBM-clear}{\code{FitGenMLVSBM$clear()}}
\item \href{#method-FitGenMLVSBM-m_step}{\code{FitGenMLVSBM$m_step()}}
\item \href{#method-FitGenMLVSBM-ve_step}{\code{FitGenMLVSBM$ve_step()}}
\item \href{#method-FitGenMLVSBM-ve_step2}{\code{FitGenMLVSBM$ve_step2()}}
\item \href{#method-FitGenMLVSBM-update_mqr}{\code{FitGenMLVSBM$update_mqr()}}
\item \href{#method-FitGenMLVSBM-do_vem}{\code{FitGenMLVSBM$do_vem()}}
\item \href{#method-FitGenMLVSBM-permute_empty_class}{\code{FitGenMLVSBM$permute_empty_class()}}
\item \href{#method-FitGenMLVSBM-xz_loglikelihood}{\code{FitGenMLVSBM$xz_loglikelihood()}}
\item \href{#method-FitGenMLVSBM-za_loglikelihood}{\code{FitGenMLVSBM$za_loglikelihood()}}
\item \href{#method-FitGenMLVSBM-reorder}{\code{FitGenMLVSBM$reorder()}}
\item \href{#method-FitGenMLVSBM-plot}{\code{FitGenMLVSBM$plot()}}
\item \href{#method-FitGenMLVSBM-show}{\code{FitGenMLVSBM$show()}}
\item \href{#method-FitGenMLVSBM-print}{\code{FitGenMLVSBM$print()}}
\item \href{#method-FitGenMLVSBM-clone}{\code{FitGenMLVSBM$clone()}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-FitGenMLVSBM-new"></a>}}
\if{latex}{\out{\hypertarget{method-FitGenMLVSBM-new}{}}}
\subsection{Method \code{new()}}{
Constructor for the FitMLVSBM class
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{FitGenMLVSBM$new(
  Q = NULL,
  A = NULL,
  X = NULL,
  M = NULL,
  directed = NULL,
  distribution = NULL,
  independent = FALSE,
  no_affiliation = NULL
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{Q}}{Vector with the number of blocks}

\item{\code{A}}{List of affiliation matrice}

\item{\code{X}}{List of adjacency matrices}

\item{\code{M}}{List of Mask matrices}

\item{\code{directed}}{Vector of boolean}

\item{\code{distribution}}{Vector of string}

\item{\code{independent}}{Boolean}

\item{\code{no_affiliation}}{A vector of boolean.
For each level, are there any nodes with no affiliations?}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
A FitGenMLVSBM object
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-FitGenMLVSBM-update_alpha"></a>}}
\if{latex}{\out{\hypertarget{method-FitGenMLVSBM-update_alpha}{}}}
\subsection{Method \code{update_alpha()}}{
Update the connection parameters for the M step
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{FitGenMLVSBM$update_alpha(m, safeguard = 2 * .Machine$double.eps)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{safeguard}}{Parameter live in a compact [safeguard, 1-safeguard]}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-FitGenMLVSBM-update_pi"></a>}}
\if{latex}{\out{\hypertarget{method-FitGenMLVSBM-update_pi}{}}}
\subsection{Method \code{update_pi()}}{
Update the  mixture parameter for the M step of level m
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{FitGenMLVSBM$update_pi(m, safeguard = 0.001)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{safeguard}}{Parameter live in a compact [safeguard, 1-safeguard]}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-FitGenMLVSBM-update_gamma"></a>}}
\if{latex}{\out{\hypertarget{method-FitGenMLVSBM-update_gamma}{}}}
\subsection{Method \code{update_gamma()}}{
Update the hierarchical mixture parameter for the M step
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{FitGenMLVSBM$update_gamma(m, safeguard = 1e-06)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{safeguard}}{Parameter live in a compact [safeguard, 1-safeguard]}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-FitGenMLVSBM-init_clustering"></a>}}
\if{latex}{\out{\hypertarget{method-FitGenMLVSBM-init_clustering}{}}}
\subsection{Method \code{init_clustering()}}{
init_clustering Initial clustering for VEM algorithm
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{FitGenMLVSBM$init_clustering(
  safeguard = 2 * .Machine$double.eps,
  method = "hierarchical",
  Z = NULL
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{safeguard}}{Parameter live in a compact [safeguard, 1-safeguard]}

\item{\code{method}}{Algorithm used to initiate the clustering, either
"spectral", "hierarchical" or "merge_split" (if \code{Z} is provided)}

\item{\code{Z}}{Initial clustering if provided}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-FitGenMLVSBM-clear"></a>}}
\if{latex}{\out{\hypertarget{method-FitGenMLVSBM-clear}{}}}
\subsection{Method \code{clear()}}{
Reset all parameters
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{FitGenMLVSBM$clear()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-FitGenMLVSBM-m_step"></a>}}
\if{latex}{\out{\hypertarget{method-FitGenMLVSBM-m_step}{}}}
\subsection{Method \code{m_step()}}{
m_step Compute the M step of the VEM algorithm
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{FitGenMLVSBM$m_step(m, safeguard = 1e-06)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{safeguard}}{Parameter live in a compact [safeguard, 1-safeguard]}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-FitGenMLVSBM-ve_step"></a>}}
\if{latex}{\out{\hypertarget{method-FitGenMLVSBM-ve_step}{}}}
\subsection{Method \code{ve_step()}}{
Compute the VE step of the VEM algorithm
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{FitGenMLVSBM$ve_step(m, threshold = 1e-06, fixPointIter = 3, safeguard = 1e-06)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{m}}{The level to be updated}

\item{\code{threshold}}{The convergence threshold}

\item{\code{fixPointIter}}{The maximum number of fixed point iterations}

\item{\code{safeguard}}{Parameter live in a compact [safeguard, 1-safeguard]}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-FitGenMLVSBM-ve_step2"></a>}}
\if{latex}{\out{\hypertarget{method-FitGenMLVSBM-ve_step2}{}}}
\subsection{Method \code{ve_step2()}}{
Compute the VE step of the VEM algorithm
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{FitGenMLVSBM$ve_step2(threshold = 1e-06, fixPointIter = 5, safeguard = 1e-06)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{threshold}}{The convergence threshold}

\item{\code{fixPointIter}}{The maximum number of fixed point iterations}

\item{\code{safeguard}}{Parameter live in a compact [safeguard, 1-safeguard]}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-FitGenMLVSBM-update_mqr"></a>}}
\if{latex}{\out{\hypertarget{method-FitGenMLVSBM-update_mqr}{}}}
\subsection{Method \code{update_mqr()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{FitGenMLVSBM$update_mqr(m)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-FitGenMLVSBM-do_vem"></a>}}
\if{latex}{\out{\hypertarget{method-FitGenMLVSBM-do_vem}{}}}
\subsection{Method \code{do_vem()}}{
Launch a Variational EM algorithm
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{FitGenMLVSBM$do_vem(
  init = "hierarchical",
  threshold = 1e-06,
  maxIter = 1000,
  fixPointIter = 10,
  safeguard = 1e-06,
  Z = NULL
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{init}}{The method for \code{self$init_clustering}}

\item{\code{threshold}}{The convergence threshold}

\item{\code{maxIter}}{The max number of VEM iterations}

\item{\code{fixPointIter}}{The max number of fixed point iterations for VE step}

\item{\code{safeguard}}{Parameter live in a compact [safeguard, 1-safeguard]}

\item{\code{Z}}{Initial clustering if provided}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-FitGenMLVSBM-permute_empty_class"></a>}}
\if{latex}{\out{\hypertarget{method-FitGenMLVSBM-permute_empty_class}{}}}
\subsection{Method \code{permute_empty_class()}}{
permute_empty_class Put empty blocks numbers at the end
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{FitGenMLVSBM$permute_empty_class(m)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-FitGenMLVSBM-xz_loglikelihood"></a>}}
\if{latex}{\out{\hypertarget{method-FitGenMLVSBM-xz_loglikelihood}{}}}
\subsection{Method \code{xz_loglikelihood()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{FitGenMLVSBM$xz_loglikelihood(m)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-FitGenMLVSBM-za_loglikelihood"></a>}}
\if{latex}{\out{\hypertarget{method-FitGenMLVSBM-za_loglikelihood}{}}}
\subsection{Method \code{za_loglikelihood()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{FitGenMLVSBM$za_loglikelihood(m)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-FitGenMLVSBM-reorder"></a>}}
\if{latex}{\out{\hypertarget{method-FitGenMLVSBM-reorder}{}}}
\subsection{Method \code{reorder()}}{
Reorder the block memberships and parameters of the networks
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{FitGenMLVSBM$reorder(order = "affiliation")}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{order}}{One of c("affiliation", "degree")}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-FitGenMLVSBM-plot"></a>}}
\if{latex}{\out{\hypertarget{method-FitGenMLVSBM-plot}{}}}
\subsection{Method \code{plot()}}{
Plot of FitMLVSBM objects
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{FitGenMLVSBM$plot(type = c("matrix"), ...)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{type}}{A string for the type of plot, just "matrix" for now}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
a ggplot2 object
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-FitGenMLVSBM-show"></a>}}
\if{latex}{\out{\hypertarget{method-FitGenMLVSBM-show}{}}}
\subsection{Method \code{show()}}{
print method
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{FitGenMLVSBM$show(type = "Multilevel Stochastic Block Model")}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{type}}{character to tune the displayed name}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-FitGenMLVSBM-print"></a>}}
\if{latex}{\out{\hypertarget{method-FitGenMLVSBM-print}{}}}
\subsection{Method \code{print()}}{
print method
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{FitGenMLVSBM$print()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-FitGenMLVSBM-clone"></a>}}
\if{latex}{\out{\hypertarget{method-FitGenMLVSBM-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{FitGenMLVSBM$clone(deep = FALSE)}\if{html}{\out{</div>}}
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
