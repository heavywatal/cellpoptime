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
#' @importFrom ape as.phylo
#' @method as.phylo data.frame
#' @param x tbl
#' @return phylo
#' @rdname tree
#' @export
as.phylo.data.frame = function(x) {
  class(x) = c("tbl_tree", class(x))
  tidytree::as.phylo(x)
}
