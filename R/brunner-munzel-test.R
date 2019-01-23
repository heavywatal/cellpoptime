# Copied from lawstat package on 2019-01-23
# https://github.com/vlyubchich/lawstat/blob/master/lawstat/R/brunner.munzel.test.R
# Cleaned by @heavywatal
brunner.munzel.test = function(x, y, alternative = c("two.sided", "greater", "less"), alpha = 0.05) {
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
  conf_half = stats::qt(1 - alpha / 2, df) * sqrt(v1 / (n1 * n2^2) + v2 / (n2 * n1^2))
  conf.int = estimate + c(-conf_half, conf_half)
  attr(conf.int, "conf.level") = (1 - alpha)
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
