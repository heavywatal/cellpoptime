#' Functions to judge false positive and false negative.
#'
#' @param predicted,actual vectors to compare
#' @rdname confusion
#' @export
false_negative = function(predicted, actual) {
  anyNA(match(actual, predicted))
}

#' @rdname confusion
#' @export
false_positive = function(predicted, actual) {
  anyNA(match(predicted, actual))
}
