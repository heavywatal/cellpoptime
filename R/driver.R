#' Functions for detecting driver mutations
#
#' `detect_driver` returns a vector of p-values.
#' @param v a pair of branch lengths (number of mutations)
#' @rdname driver
#' @export
detect_driver = function(v) {
  stats::ppois(v, min(v), lower.tail=FALSE)
}

#' `filter_origins` filters candidate origins by adjusted p-values.
#' @param x scaled tbl_tree
#' @param method passed to `p.adjust()`
#' @param q threshold
#' @rdname driver
#' @export
filter_origins = function(x, method = "fdr", q = 0.05) {
  dplyr::filter(x, stats::p.adjust(.data$p_driver, method = method) < q)
}

#' @description
#' `mutant_labels` extracts mutated nodes in tbl_tree
#' @rdname driver
#' @export
mutant_labels = function(x) {
  origins = filter_origins(x)$node
  group_clade(x, origins) %>%
    dplyr::filter(.data$is_tip, .data$group != 0L) %>%
    dplyr::pull("label") %>%
    as.integer()
}
