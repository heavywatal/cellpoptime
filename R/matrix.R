#' Functions for genotype matrix
#'
#' @description
#' `infer_tree` infers a phylogenetic tree from sampled genotypes.
#' @param mtrx integer matrix; rows are samples, columns are sites.
#' @rdname matrix
#' @export
infer_tree = function(mtrx) {
  mtrx %>%
    add_outgroup() %>%
    infer_rooted_tree() %>%
    remove_outgroup()
}

infer_rooted_tree = function(mtrx) {
  mtrx %>%
    stats::dist(method = "manhattan") %>%
    ape::fastme.bal() %>%
    ape::root("0") %>%
    tidytree::as_tibble()
}

add_outgroup = function(mtrx) {
  rbind(`0` = 0L, mtrx)
}

remove_outgroup = function(x) {
  if ("0" %in% x$label) {
    dplyr::filter(x, .data$label != "0" | is.na(.data$label)) %>%
      dplyr::mutate(parent = .data$parent - 1L, node = .data$node - 1L)
  } else {
    x
  }
}
