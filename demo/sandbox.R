library(tidyverse)
library(wtl)
library(tumopp)
library(ggtree)
refresh('cellpoptime')

.ggtree = function(.tbl) {
  ggtree(.tbl) +
  geom_tree() +
  geom_tiplab() +
  wtl::theme_wtl()
}

.param = c('-D2', '-k100', '-N16', '-U8', '--mb=99', '--ms1mut')
samples = purrr::rerun(4L, {
  (tumopp::mslike(16L, .param) %>% wtl::parse_ms())[[1L]]
}) %>% print()

# samples = wtl::ms(6L, 4L, theta = 100) %>% wtl::parse_ms() %>% print()
trees = samples %>% purrr::map(cellpoptime::infer_rooted_tree) %>% print()
trees %>%
  cellpoptime::as_multiphylo() %>%
  .ggtree() +
  facet_wrap(~.id)

.base = trees[[3]] %>% cellpoptime::remove_outgroup() %>% print()

add_extra_columns = function(x) {
  x %>%
    dplyr::arrange(parent) %>%
    dplyr::mutate(
      is_tip = !is.na(label),
      mutations = branch.length,
      branch.length = pmax(branch.length, 0.01),
      term_length = 0,
      children = list(NULL)
    )
}
.base %>% add_extra_columns()

filter_scale_tips = function(x) {
  x %>%
    dplyr::filter(is_tip) %>%
    dplyr::group_by(parent) %>%
    dplyr::filter(n() > 1L) %>%
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
.base %>% add_extra_columns() %>% filter_scale_tips()

nest_tippairs = function(x) {
  nested = filter_scale_tips(x) %>%
    tidyr::nest(-parent, -term_length, .key="children") %>%
    print()
  x %>%
    dplyr::filter(is.na(branch.length) | !parent %in% nested$parent) %>%
    print() %>%
    dplyr::left_join(nested, by=c(node="parent"), suffix=c("", ".y")) %>%
    dplyr::mutate(
      is_tip = is_tip | node %in% nested$parent,
      term_length = pmax(term_length, term_length.y, na.rm=TRUE),
      children = ifelse(purrr::map_lgl(children, is.null), children.y, children),
      term_length.y = NULL, children.y = NULL)
}
.base %>% add_extra_columns() %>% nest_tippairs() #%>% nest_tippairs() %>% nest_tippairs()

unnest_children = function(x) {
  .outer = x %>%
    dplyr::transmute(parent=node, children) %>%
    dplyr::filter(!purrr::map_lgl(children, is.null)) %>%
    tidyr::unnest()
  .inner = x %>% dplyr::mutate(children=list(NULL))
  dplyr::bind_rows(.outer, .inner) %>%
    dplyr::mutate(is_tip = !is.na(label))
}

rescale_branches = function(x) {
  num_edges = nrow(x)
  while (nrow(x) > 1L) {
    x = nest_tippairs(x)
  }
  while (nrow(x) < num_edges) {
    x = unnest_children(x)
  }
  x %>%
    # dplyr::select(parent, node, branch.length, label) %>%
    dplyr::arrange(node)
}

.shrunk = .base %>% add_extra_columns() %>% rescale_branches() %>% print()

.p = list(.base, .shrunk) %>%
  as_multiphylo() %>%
  .ggtree() +
  facet_wrap(~.id, ncol=1)
.p
ggsave('scaletree.pdf', .p, width=7, height=9.9)
