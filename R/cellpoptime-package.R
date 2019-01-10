#' cellpoptime: analyze and visualize phylogenetic tree of cell population
#' @useDynLib cellpoptime, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @keywords internal
"_PACKAGE"

.onUnload = function(libpath) {
  library.dynam.unload("cellpoptime", libpath)
}
