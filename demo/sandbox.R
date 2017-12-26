library(tidyverse)
library(wtl)
library(ggtree)
refresh('cellpoptime')

samples = ms(6L, 4L, theta = 100) %>% print()

trees = samples %>% purrr::map(infer_rooted_tree) %>% print()

.ggtree = function(.tbl) {
    ggtree(.tbl) +
    geom_tree() +
    geom_tiplab() +
    theme_bw()
}

trees %>%
  as_multiphylo() %>%
  .ggtree() +
  facet_wrap(~.id)

filter_scale_tips = function(x) {
  x %>%
    dplyr::filter(!is.na(branch.length)) %>%
    dplyr::group_by(parent) %>%
    dplyr::filter(all(isTip)) %>%
    print() %>%
    dplyr::mutate(total_length = min(branch.length + term_length)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      scale = total_length / (branch.length + term_length),
      branch.length = branch.length * scale,
      term_length = term_length * scale,
      children = purrr::map2(children, scale, rescale_children)
    ) %>%
    print() %>%
    dplyr::mutate(term_length = total_length) %>%
    dplyr::select(-scale, -total_length)
}

nest_tippairs = function(x) {
  nested = filter_scale_tips(x) %>%
    print() %>%
    tidyr::nest(-parent, -term_length, .key="children") %>%
    print()
  x %>%
    dplyr::filter(is.na(branch.length) | !parent %in% nested$parent) %>%
    dplyr::left_join(nested, by=c(node="parent"), suffix=c("", ".y")) %>%
    dplyr::mutate(
      isTip = isTip | node %in% nested$parent,
      term_length = pmax(term_length, term_length.y, na.rm=TRUE),
      children = ifelse(purrr::map_lgl(children, is.null), children.y, children),
      term_length.y = NULL, children.y = NULL)
}

unnest_children = function(x) {
  .outer = x %>%
    dplyr::transmute(parent=node, children) %>%
    dplyr::filter(!purrr::map_lgl(children, is.null)) %>%
    tidyr::unnest()
  .inner = x %>% dplyr::mutate(children=list(NULL))
  bind_rows(.outer, .inner) %>%
    dplyr::mutate(isTip = !is.na(label))
}

rescale_branches = function(x) {
  x = x %>%
    dplyr::arrange(parent) %>%
    dplyr::mutate(
      mutations = branch.length,
      branch.length = pmax(branch.length + 0.01),
      term_length = 0,
      children = list(NULL)
    )
  num_edges = nrow(x)
  while (nrow(x) > 1L) {
    x = nest_tippairs(x)
  }
  while (nrow(x) < num_edges) {
    x = unnest_children(x)
  }
  x %>%
    dplyr::select(-term_length, -children) %>%
    dplyr::arrange(node)
}

.base = trees[[4]] %>% remove_outgroup() %>% print()
.shrunk = .base %>% rescale_branches() %>% print()

list(.base, .shrunk) %>%
  as_multiphylo() %>%
  .ggtree() +
  facet_wrap(~.id, ncol=1)
