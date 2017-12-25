library(tidyverse)
library(wtl)
library(ggtree)
loadNamespace("cowplot")
refresh('cellpoptime')

samples = ms(6L, 4L, theta = 100) %>% print()

trees = samples %>% purrr::map(infer_rooted_tree) %>% print()

.plot = function(.tbl) {
  as.phylo(.tbl) %>%
    ggplot(aes(x, y)) +
    geom_tree() +
    geom_tiplab() +
    theme_bw()
}

.plts = purrr::map(trees, .plot)

cowplot::plot_grid(plotlist = .plts, ncol = 2L)


remove_dummy_root = function(x) {
  if ("0" %in% x$label) {
    dplyr::filter(x, label != "0" | is.na(label)) %>%
      dplyr::mutate(parent = parent - 1L, node = node - 1L)
  } else {
    x
  }
}

.base = trees[[1]] %>% remove_dummy_root() %>% print()


nest_tippairs = function(x) {
  paired = dplyr::filter(x, !is.na(branch.length)) %>%
    dplyr::group_by(parent) %>%
    dplyr::filter(all(isTip)) %>% print() %>%
    dplyr::mutate(
      branch.length = ifelse(branch.length > 0, min(branch.length + coalesce(term_length, 0)) * branch.length / (branch.length + coalesce(term_length, 0)), 0),
      term_length = min(branch.length + coalesce(term_length, 0))) %>%
    dplyr::ungroup() %>%
    print()
  nested = paired %>% tidyr::nest(-parent, -term_length, .key="children") %>% print()
  x %>%
    dplyr::filter(!node %in% paired$node) %>%
    dplyr::left_join(nested, by=c(node="parent"), suffix=c("", ".y")) %>%
    dplyr::mutate(isTip = isTip | node %in% nested$parent,
      term_length = coalesce(term_length, term_length.y),
      children = ifelse(purrr::map_lgl(children, is.null), children.y, children),
      term_length.y = NULL, children.y = NULL)
  #TODO: shrink children as well
}

.nested = .base %>%
  dplyr::mutate(mutations = branch.length, term_length=NA_real_, children=list(NULL)) %>%
  nest_tippairs() %>% print() %>%
  nest_tippairs() %>% print() %>%
  nest_tippairs() %>% print()

unnest_pairs = function(x) {
  .nested = x %>% dplyr::transmute(parent=node, children) %>% dplyr::filter(purrr::map_lgl(children, ~!is.null(.x)))
  bind_rows(
    x %>% dplyr::mutate(children=list(NULL)),
    .nested %>% tidyr::unnest()
  )
}

.shrunk = .nested %>%
  unnest_pairs() %>% print() %>%
  unnest_pairs() %>% print() %>%
  unnest_pairs() %>% print() %>%
  dplyr::select(-term_length, -children) %>%
  print()

cowplot::plot_grid(.base %>% .plot(), .shrunk %>% .plot(), ncol = 2L)
