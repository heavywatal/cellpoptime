#' Functions for visualization
#'
#' @details
#' `plot_tree` draws tbl_tree.
#' @param data tbl_tree
#' @param ... passed to `ggtree::geom_tree()` aes mapping
#' @rdname plot
#' @export
plot_tree = function(data, ...) {
  ggplot2::ggplot(data, ggplot2::aes_(~x, ~y)) +
    ggtree::geom_tree(mapping = ggplot2::aes_(colour = ~ try_null(group), ...))
}

#' @details
#' `geom_nodename` annotates tree nodes.
#' @rdname plot
#' @export
geom_nodename = function() {
  list(
    ggtree::geom_text2(ggplot2::aes_(label = ~ paste0(node, ":", dplyr::coalesce(label, ""))), hjust = -.1, colour = "#666666"),
    ggplot2::scale_x_continuous(expand = ggplot2::expand_scale(mult = c(0.06, 0.16)))
  )
}

#' @details
#' `geom_driver` annotates driver mutations.
#' @rdname plot
#' @export
geom_driver = function() {
  list(
    ggplot2::geom_point(ggplot2::aes_(alpha = ~ pmin(-log10(p_driver), 6)), size = 2),
    ggtree::geom_text2(ggplot2::aes_(label = ~ sprintf("%.02g", p_driver)), hjust = 1.1),
    ggplot2::scale_alpha(range = c(0.05, 1), limits = c(1, 6), guide = FALSE, na.value = 0)
  )
}

plot_tree_dev = function(data) {
  plot_tree(data) +
    geom_nodename() +
    geom_driver() +
    ggplot2::theme(legend.position = "none")
}

try_null = function(...) tryCatch(..., error = function(e) NULL)

annotate_node = function(p, node, shape = 15, size = 4, colour = "orange", after = 0L) {
  f = function(x) dplyr::filter(x, .data$node == !!node)
  insert_layer(p, ggplot2::geom_point(data = f, shape = shape, size = size, colour = colour), after = after)
}

insert_layer = function(p, ..., after = 0L) {
  if (after < 0L) {
    after = length(p$layers) + after
  }
  p$layers = append(p$layers, list(...), after = after)
  p
}
