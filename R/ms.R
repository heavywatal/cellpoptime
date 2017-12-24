#' Execute ms as an external program
#' @param nsam number of samples
#' @param nreps number of repeats
#' @param theta population mutation rate
#' @return string vector
#' @rdname ms
#' @export
ms = function(nsam=4L, nreps=2L, theta=5.0) {
  cmd = paste("ms", nsam, nreps, "-t", theta)
  system(cmd, intern = TRUE, ignore.stderr = TRUE) %>%
    parse_ms()
}

#' Read ms output and return list of '01' string vector
#' @param lines output of system('ms')
#' @return list of string vector
#' @rdname ms
parse_ms = function(lines) {
  mobj = stringr::str_match(lines[1L], "^ms\\s+(\\d+)\\s+(\\d+)")[1L, -1L]
  nsam = as.integer(mobj[1L])
  nreps = as.integer(mobj[2L])
  utils::tail(lines, -2L) %>%
    split(rep(seq_len(nreps), each = nsam + 4L)) %>%
    purrr::map(utils::tail, -4L)
}
