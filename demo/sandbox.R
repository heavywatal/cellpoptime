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

scale_children = function(x, scale) {
  if (is.null(x)) {
    x
  } else {
    dplyr::mutate(x,
      branch.length = branch.length * scale,
      children = purrr::map(children, scale_children, scale)
    )
  }
}

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
      children = purrr::map2(children, scale, scale_children)
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

.base = trees[[1]] %>% remove_dummy_root() %>% print()

.nested = .base %>%
  dplyr::arrange(parent) %>%
  dplyr::mutate(mutations = branch.length, term_length=0, children=list(NULL)) %>%
  nest_tippairs() %>% print() %>%
  nest_tippairs() %>% print() %>%
  nest_tippairs() %>% print() %>%
  nest_tippairs() %>% print()

unnest_pairs = function(x) {
  .outer = x %>%
    dplyr::transmute(parent=node, children) %>%
    dplyr::filter(!purrr::map_lgl(children, is.null)) %>%
    tidyr::unnest()
  .inner = x %>% dplyr::mutate(children=list(NULL))
  bind_rows(.outer, .inner) %>%
    dplyr::mutate(isTip = !is.na(label))
}

.shrunk = .nested %>%
  unnest_pairs() %>% print() %>%
  unnest_pairs() %>% print() %>%
  unnest_pairs() %>% print() %>%
  unnest_pairs() %>% print() %>%
  dplyr::select(-term_length, -children) %>%
  dplyr::arrange(node) %>%
  print()

cowplot::plot_grid(
  .base %>% .plot() + coord_cartesian(xlim = c(0, 400)),
  .shrunk %>% .plot() + coord_cartesian(xlim = c(0, 400)),
  nrow = 2L)
