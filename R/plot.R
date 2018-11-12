#' Functions for visualization
#'
#' @description
#' `fortify_cellpop` prepares plottable data.frame.
#' @param model output of `scale_branches()`
#' @param data tbl_tree
#' @param ... unused
#' @rdname plot
#' @export
fortify_cellpop = function(model, data, ...) {
    if (missing(data)) data = model
    mutant = dplyr::top_n(model, 1L, dplyr::desc(.data$p_driver))
    meta_info = group_clade(model, mutant$node) %>%
      dplyr::select(.data$node, .data$mutations, .data$exp_sibs, .data$p_driver, .data$group)
    data = if (is.data.frame(data)) {
      as_phylo(data)
    } else if (is.list(data)) {
      as_multiphylo(data)
    } else {
      stop("Unknown class(data): ", class(data))
    }
    data %>%
      ggtree::fortify() %>%
      dplyr::left_join(meta_info, by = "node")
}

#' @description
#' `plot_driver` draws tbl_tree.
#' @param verbose add verbose information
#' @rdname plot
#' @export
plot_driver = function(data, verbose = FALSE) {
  p = ggplot2::ggplot(data, ggplot2::aes_(~x, ~y)) +
    ggtree::geom_tree(ggplot2::aes_(colour = ~group)) +
    ggplot2::theme_void() +
    ggplot2::theme(legend.position = "none")
  if (verbose) {
    p = p +
      ggplot2::geom_point(ggplot2::aes_(alpha = ~-log10(p_driver), colour = ~group), size = 2) +
      ggtree::geom_text2(ggplot2::aes_(label = ~paste0(node, ":", label)), hjust=-.1, colour = '#666666') +
      ggtree::geom_text2(ggplot2::aes_(label = ~sprintf('%.02g', p_driver)), hjust=1.1) +
      ggplot2::scale_x_continuous(expand = ggplot2::expand_scale(mult = c(0.06, 0.16))) +
      ggplot2::scale_alpha(range = c(0, 1), guide = FALSE)
  }
  p
}
