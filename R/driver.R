#' Functions for detecting driver mutations
#
#' `detect_driver` returns a vector of p-values.
#' @param x a pair of branch lengths (number of mutations)
#' @rdname driver
#' @export
detect_driver = function(x) {
  ppois(x, min(x), lower.tail=FALSE)
}
