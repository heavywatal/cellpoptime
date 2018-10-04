#' Functions for genotype matrix
#'
#' @description
#' `infer_rooted_tree` infers rooted tree of sampled genotypes.
#' @param mtrx integer matrix; rows are samples, columns are sites.
#' @rdname matrix
#' @export
infer_rooted_tree = function(mtrx) {
  mtrx %>%
    add_outgroup() %>%
    stats::dist(method = "manhattan") %>%
    ape::fastme.bal() %>%
    ape::root("0") %>%
    tidytree::as_data_frame()
}

#' @description
#' `add_outgroup` adds wild type row with genotype 000 and name 0.
#' @rdname matrix
#' @export
add_outgroup = function(mtrx) {
  mtrx = rbind(0L, mtrx)
  rownames(mtrx) = seq_len(nrow(mtrx)) - 1L
  mtrx
}
