#' Take normal data.frame as tbl_tree for simplicity
#' @param x tbl_tree
#' @rdname tree
#' @export
as_phylo = function(x) {
  class(x) = c("tbl_tree", class(x))
  tidytree::as.phylo(x)
}

#' Add class name
#' @param tbls list of tbl_tree
#' @rdname tree
#' @export
as_multiphylo = function(tbls) {
  tbls = purrr::map(tbls, as_phylo)
  class(tbls) = c("multiPhylo", class(tbls))
  tbls
}

#' @param node ancestral node
#' @rdname tree
#' @export
group_clade = function(x, node) {
  class(x) = c("tbl_tree", class(x))
  tidytree::groupClade(x, node)
}

#' @rdname tree
#' @export
upstream_corner = function(x, node) {
  row_n = dplyr::filter(x, .data$node == !!node)
  row_p = dplyr::filter(x, .data$node == row_n$parent)
  row_n %>%
    dplyr::transmute(.data$parent, .data$node, .data$y, .data$isTip) %>%
    dplyr::mutate(x = row_p$x)
}

#' Remove dummy outgroup
#' @rdname tree
#' @export
remove_outgroup = function(x) {
  if ("0" %in% x$label) {
    dplyr::filter(x, .data$label != "0" | is.na(.data$label)) %>%
      dplyr::mutate(parent = .data$parent - 1L, node = .data$node - 1L)
  } else {
    x
  }
}
