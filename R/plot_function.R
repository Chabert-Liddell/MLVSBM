#' Plot a matrix representation of the network.
#'
#' @description This function does not plot the affiliation relationship, only
#' the levels. The modes are ordered by blocks and a color gradient represent
#' the probability of each entries.
#' @param net A fitted network. Either a FitSBM or a FitMLVSBM object
#' @param level In \code{c("lower", "upper")}. Only for FitMLVSBM object.
#' If \code{NULL} then both level will be plotted.
#'
#' @return A \code{ggplot2} plot. The matrix representation of one or each level.
#'
#'
#' @import cowplot
#' @import ggraph
#' @importFrom  dplyr `%>%`
#' @export
#'
#' @examples
plot_matrix <- function(net, level = NULL) {
  if (!requireNamespace("ggraph", quietly = TRUE)) {
    stop("Package \"ggraph\" needed for this function to work. Please install it.",
         call. = FALSE)
  }

  ggdraw()+
    draw_plot(
      t(fit$X_hat$I * lower_level) %>%
        as_tbl_graph() %>%
        mutate(group = fit$Z$I) %>%
        ggraph('matrix', sort.by = group)+
        geom_edge_point(aes(color = weight))+
        scale_edge_color_gradient(name = 'Prob', low = "#deebf7", high = "#08519c") +
        geom_hline(yintercept = cumsum(table(fit$Z$I))[-fit$nb_clusters$I]+.5) +
        geom_vline(xintercept = cumsum(table(fit$Z$I))[-fit$nb_clusters$I]+.5) +
        scale_y_reverse() +
        coord_fixed() +
        theme_graph() , 0, 0, .5, 1
    ) +
    draw_plot(
      t(fit$X_hat$O * upper_level) %>%
        as_tbl_graph() %>%
        ggraph('matrix', sort.by = fit$Z$O) +
        geom_edge_point(aes(color = weight)) +
        geom_hline(yintercept = cumsum(table(fit$Z$O))[-fit$nb_clusters$O]+.5) +
        geom_vline(xintercept = cumsum(table(fit$Z$O))[-fit$nb_clusters$O]+.5) +
        scale_edge_colour_gradient(name = 'Prob',low = "#fcbba1", high = "#67000d") +
        scale_y_reverse() +
        coord_fixed() +
        theme_graph(), x = 0.5, y = 0, width = .5, height = 1
    ) +
    draw_plot_label(label = c("Lower level", "Upper level"), x = c(0, .5), y = c(.9, .9))
}
# my_fun <- function(a, b) {

#   }
# }
