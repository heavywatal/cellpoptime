#' Filter tbl_tree by various criteria
#'
#' @param .data tbl_tree
#' @param node integer ID
#' @rdname filter
#' @export
filter_root = function(.data) {
  .data[.data$parent == .data$node, ]
}

#' @rdname filter
#' @export
filter_children = function(.data, node) {
  .data[.data$parent == node & .data$node != node, ]
}

#' @details
#' `filter_offspring` and `filter_clade` are similar functions;
#' the former does not include the specified node itself, but the latter does.
#' @rdname filter
#' @export
filter_offspring = function(.data, node) {
  .data[is_offstring(.data, node), ]
}

#' @rdname filter
#' @export
filter_clade = function(.data, node) {
  .data[is_clade(.data, node), ]
}

is_offstring = function(.data, node) {
  # temporary vectors for speed
  data_parent = .data$parent
  data_node = .data$node
  included = added = (data_parent == node) & (data_node != node)
  while (any(added)) {
    added = data_parent %in% data_node[added]
    included = included | added
  }
  included
}

is_clade = function(.data, node) {
  is_offstring(.data, node) | (.data$node == node)
}

filter_root_children = function(.data) {
  filter_children(.data, filter_root(.data)$node)
}
