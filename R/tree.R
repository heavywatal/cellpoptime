#' Infer rooted tree of sampled genotypes
#' @param genotypes string vector
#' @return tibble
#' @rdname tree
#' @export
infer_rooted_tree = function(genotypes) {
    genotypes %>%
    as_int_matrix() %>%
    add_outgroup() %>%
    stats::dist(method='manhattan') %>%
    ape::fastme.bal() %>%
    ape::root('0') %>%
    tidytree::as_data_frame()
}
