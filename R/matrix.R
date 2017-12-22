#' Convert 01 string vector to integer matrix
#' @param x string vector like '101'
#' @return integer matrix
#' @rdname matrix
#' @export
as_int_matrix = function(x) {
    x = stringr::str_split(x, '') %>% purrr::map(as.integer)
    matrix(unlist(x), nrow=length(x), byrow=TRUE)
}

#' Add wild type row with genotype 000 and name 0
#' @param mtrx integer matrix
#' @return integer matrix
#' @rdname matrix
#' @export
add_outgroup = function(mtrx) {
    mtrx = rbind(0L, mtrx)
    rownames(mtrx) = seq_len(nrow(mtrx)) - 1L
    mtrx
}
