#' Functions for detecting driver mutations
#
#' `detect_driver` returns a vector of p-values.
#' @param exp_desc expected number of descendants
#' @rdname driver
#' @export
detect_driver = function(exp_desc) {
  p = exp_desc * 0 + 1
  p[which.max(exp_desc)] = p_binom(exp_desc - 1)
  p
}

#' @description
#' `p_binom` is a shortcut.
#' @param events integer
#' @rdname driver
#' @export
p_binom = function(events) {
  stats::pbinom(min(events), sum(events), 0.5)
}
