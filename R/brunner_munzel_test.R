#' Brunner-Munzel Test
#'
#' The original code was copied from lawstat v3.2 on 2019-01-23.
#' @source <https://github.com/vlyubchich/lawstat>
#' @inheritParams stats::t.test
#' @seealso `stats::t.test`, `stats::wilcox.test`
#' @rdname brunner_munzel_test
#' @export
brunner_munzel_test = function(...) UseMethod("brunner_munzel_test")

#' @rdname brunner_munzel_test
#' @export
brunner_munzel_test.default = function(x, y,
                                       alternative = c("two.sided", "greater", "less"),
                                       conf.level = 0.95, ...) {
  alternative = match.arg(alternative)
  data.name = paste(deparse(substitute(x)), "and", deparse(substitute(y)))
  x = stats::na.omit(x)
  y = stats::na.omit(y)
  n1 = length(x)
  n2 = length(y)
  r1 = rank(x)
  r2 = rank(y)
  r = rank(c(x, y))
  m1 = mean(r[seq_len(n1)])
  m2 = mean(r[-seq_len(n1)])
  v1 = sum((r[seq_len(n1)] - r1 - m1 + (n1 + 1) / 2)^2) / (n1 - 1)
  v2 = sum((r[-seq_len(n1)] - r2 - m2 + (n2 + 1) / 2)^2) / (n2 - 1)
  statistic = n1 * n2 * (m2 - m1) / (n1 + n2) / sqrt(n1 * v1 + n2 * v2)
  df = ((n1 * v1 + n2 * v2)^2) / (((n1 * v1)^2) / (n1 - 1) + ((n2 * v2)^2) / (n2 - 1))
  p.value = switch(alternative,
    two.sided = 2 * stats::pt(abs(statistic), df, lower.tail = FALSE),
    greater = stats::pt(statistic, df, lower.tail = TRUE),
    less = stats::pt(statistic, df, lower.tail = FALSE)
  )
  estimate = (m2 - (n2 + 1) / 2) / n1
  conf_half = stats::qt(0.5 + 0.5 * conf.level, df) * sqrt(v1 / (n1 * n2^2) + v2 / (n2 * n1^2))
  conf.int = estimate + c(-conf_half, conf_half)
  attr(conf.int, "conf.level") = conf.level
  structure(list(
    statistic = stats::setNames(statistic, "Brunner-Munzel Test Statistic"),
    parameter = stats::setNames(df, "df"),
    p.value = p.value,
    conf.int = stats::setNames(conf.int, c("lower", "upper")),
    estimate = stats::setNames(estimate, "P(X<Y)+.5*P(X=Y)"),
    method = "Brunner-Munzel Test",
    data.name = data.name
  ), class = "htest")
}

#' @rdname brunner_munzel_test
#' @export
brunner_munzel_test.formula = function(formula, data, ...) {
  mf = stats::model.frame(formula, data = data)
  response = attr(attr(mf, "terms"), "response")
  group = factor(mf[[-response]])
  stopifnot(nlevels(group) == 2L)
  data = stats::setNames(split(mf[[response]], group), c("x", "y"))
  result = do.call("brunner_munzel_test", c(data, list(...)))
  result$data.name = paste(names(mf), collapse = " by ")
  result
}
