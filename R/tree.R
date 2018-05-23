#' Infer rooted tree of sampled genotypes
#' @param genotypes string vector
#' @return tibble
#' @rdname tree
#' @export
infer_rooted_tree = function(genotypes) {
  genotypes %>%
    as_int_matrix() %>%
    add_outgroup() %>%
    stats::dist(method = "manhattan") %>%
    ape::fastme.bal() %>%
    ape::root("0") %>%
    tidytree::as_data_frame()
}

#' Take normal data.frame as tbl_tree for simplicity
#' @param x tbl
#' @return phylo
#' @rdname tree
#' @export
as_phylo = function(x) {
  class(x) = c("tbl_tree", class(x))
  tidytree::as.phylo(x)
}

#' Add class name
#' @param tbls list of tbl
#' @return multiPhylo
#' @rdname tree
#' @export
as_multiphylo = function(tbls) {
  tbls = purrr::map(tbls, as_phylo)
  class(tbls) = c("multiPhylo", class(tbls))
  tbls
}

#' Remove dummy outgroup
#' @return tibble
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

#' Rescale chilren branches recursively
#' @param scale numeric
#' @return tibble
#' @rdname tree
#' @export
rescale_children = function(x, scale) {
  if (is.null(x)) {
    x
  } else {
    dplyr::mutate(
      x,
      branch.length = .data$branch.length * scale,
      children = purrr::map(.data$children, rescale_children, scale)
    )
  }
}
