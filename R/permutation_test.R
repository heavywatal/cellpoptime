#' Permutation Test
#'
#' @inheritParams stats::t.test
#' @param n performs bootstrapping if n > 0
#' @seealso `stats::t.test`, `stats::wilcox.test`
#' @rdname permutation_test
#' @export
permutation_test = function(...) UseMethod("permutation_test")

#' @rdname permutation_test
permutation_test.default = function(x, y, n = 1000L, alternative = c("two.sided", "greater", "less")) {
  alternative = match.arg(alternative)
  data.name = paste(deparse(substitute(x)), "and", deparse(substitute(y)))
  statistic = mean(x) - mean(y)
  min_len = min(length(x), length(y))
  v = c(x, y)
  bootstrapping = n > 0L
  if (bootstrapping) {
    diffs = .permutate_bootstrap(v, min_len, n = n)
  } else {
    diffs = .permutate_all(v, min_len)
    n = length(diffs)
  }
  p.value = switch(alternative,
    two.sided = 2 * min(sum(diffs < statistic), sum(statistic < diffs)) / n,
    greater = sum(statistic < diffs) / n,
    less = sum(diffs < statistic) / n
  )
  result = list(
    statistic = stats::setNames(statistic, "Difference in means"),
    parameter = stats::setNames(n, "n"),
    p.value = p.value,
    method = "Permutation Test",
    data.name = data.name
  )
  if (bootstrapping) {
    result$method = paste(result$method, " with Bootstrapping")
  }
  structure(result, class = "htest")
}

#' @rdname permutation_test
permutation_test.formula = function(formula, data, ...) {
  mf = stats::model.frame(formula, data = data)
  response = attr(attr(mf, "terms"), "response")
  group = factor(mf[[-response]])
  stopifnot(nlevels(group) == 2L)
  data = stats::setNames(split(mf[[response]], group), c("x", "y"))
  result = do.call("permutation_test", c(data, list(...)))
  result$data.name = paste(names(mf), collapse = " by ")
  result
}

.permutate_bootstrap = function(v, size, n) {
  len = length(v)
  diffs = numeric(n)
  for (i in seq_len(n)) {
    idx = sample.int(len, size, replace = FALSE)
    diffs[[i]] = mean(v[idx]) - mean(v[-idx])
  }
  diffs
}

.permutate_all = function(v, size) {
  len = length(v)
  combin = combinations(len, size)
  n = nrow(combin)
  diffs = numeric(n)
  for (i in seq_len(n)) {
    idx = combin[i, ]
    diffs[[i]] = mean(v[idx]) - mean(v[-idx])
  }
  diffs
}

# simplified version of gtools::combinations
combinations = function(n, r) {
  v = seq_len(n)
  subfun = function(n, r, v) {
    if (r == 0L)
      integer(0L)
    else if (r == 1L)
      matrix(v, n, 1L)
    else if (r == n)
      matrix(v, 1L, n)
    else
      rbind(
        cbind(v[1L], Recall(n - 1L, r - 1L, v[-1L])),
        Recall(n - 1L, r, v[-1L])
      )
  }
  subfun(n, r, v)
}
