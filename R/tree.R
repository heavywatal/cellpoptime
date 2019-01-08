#' Functions to process tidytree objects
#'
#' @details
#' `as_tbl_tree` is a workaround until dplyr 0.8.0 is released.
#' @param x data.frame
#' @rdname tree
#' @export
as_tbl_tree = function(x) {
  class(x) = c("tbl_tree", "tbl_df", "tbl", "data.frame")
  x
}

#' @details
#' `as_phylo` is a shortcut to convert data.frame to ape::phylo.
#' @rdname tree
#' @export
as_phylo = function(x) {
  tidytree::as.phylo(as_tbl_tree(x))
}

#' @details
#' `as_multiphylo` adds class name.
#' @param tbls list of tbl_tree
#' @rdname tree
#' @export
as_multiphylo = function(tbls) {
  tbls = purrr::map(tbls, as_phylo)
  class(tbls) = c("multiPhylo", class(tbls))
  tbls
}

#' @details
#' `group_clade` is a shortcut to apply `tidytree::groupClade` to data.frame.
#' @param node ancestral node
#' @rdname tree
#' @export
group_clade = function(x, node) {
  tidytree::groupClade(as_tbl_tree(x), node)
}

#' @details
#' `upstream_corner` finds the joint between a node and its parent
#' @rdname tree
#' @export
upstream_corner = function(x, node) {
  row_n = dplyr::filter(x, .data$node == !!node)
  row_p = dplyr::filter(x, .data$node == row_n$parent)
  row_n %>%
    dplyr::transmute(.data$parent, .data$node, .data$y, .data$isTip) %>%
    dplyr::mutate(x = row_p$x)
}
