#' Detect driver mutations
#' @param exp_desc expected number of descendants
#' @rdname driver
#' @export
detect_driver = function(exp_desc) {
  p = exp_desc * 0 + 1
  p[which.max(exp_desc)] = p_binom(exp_desc - 1)
  p
}

#' @param events integer
#' @rdname driver
#' @export
p_binom = function(events) {
  stats::pbinom(min(events), sum(events), 0.5)
}
