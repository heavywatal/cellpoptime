// [[Rcpp::plugins(cpp14)]]
#include <Rcpp.h>

//' Vectorized version of `is.null`
//' @param x list
//' @rdname rcpp
//' @export
// [[Rcpp::export]]
std::vector<bool> apply_is_null(const Rcpp::List& x) {
  size_t n = x.size();
  std::vector<bool> res(n);
  for (size_t i = 0u; i < n; ++i) {
    res[i] = (x[i] == R_NilValue);
  }
  return res;
}
