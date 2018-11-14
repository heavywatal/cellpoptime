#' Functions for visualization
#'
#' @description
#' `fortify_cellpop` prepares plottable data.frame.
#' @param model output of `scale_branches()`
#' @param data tbl_tree
#' @rdname plot
#' @export
fortify_cellpop = function(model, data) {
    if (missing(data)) data = model
    mutant = filter_origins(model)$node
    meta_info = group_clade(model, mutant) %>%
      dplyr::select(.data$node, .data$mutations, .data$exp_sibs, .data$p_driver, .data$group)
    ggtree_fortify(data) %>%
      dplyr::left_join(meta_info, by = "node")
}

#' @description
#' `ggtree_fortify` prepares plottable data.frame.
ggtree_fortify = function(data) {
    data = if (is.data.frame(data)) {
      as_phylo(data)
    } else if (is.list(data)) {
      as_multiphylo(data)
    } else {
      stop("Unknown class(data): ", class(data))
    }
    ggtree::fortify(data)
}

#' @description
#' `plot_tree` draws tbl_tree.
#' @param ... passed to `ggtree::geom_tree()` aes mapping
#' @rdname plot
#' @export
plot_tree = function(data, ...) {
  ggplot2::ggplot(data, ggplot2::aes_(~x, ~y)) +
    ggtree::geom_tree(mapping = ggplot2::aes_(colour = ~try_null(group), ...))
}

#' @description
#' `geom_nodename` annotates tree nodes.
#' @rdname plot
#' @export
geom_nodename = function() {
  list(
    ggtree::geom_text2(ggplot2::aes_(label = ~paste0(node, ":", dplyr::coalesce(label, ""))), hjust=-.1, colour = '#666666'),
    ggplot2::scale_x_continuous(expand = ggplot2::expand_scale(mult = c(0.06, 0.16)))
  )
}

#' @description
#' `geom_driver` annotates driver mutations.
#' @rdname plot
#' @export
geom_driver = function() {
  list(
    ggplot2::geom_point(ggplot2::aes_(alpha = ~-log10(p_driver)), size = 2),
    ggtree::geom_text2(ggplot2::aes_(label = ~sprintf('%.02g', p_driver)), hjust=1.1),
    ggplot2::scale_alpha(range = c(0.05, 1), guide = FALSE, na.value = 0)
  )
}

plot_tree_dev = function(data) {
  plot_tree(data) +
  geom_nodename() +
  geom_driver() +
  ggplot2::theme(legend.position = "none")
}

try_null = function(...) tryCatch(..., error = function(e) NULL)
