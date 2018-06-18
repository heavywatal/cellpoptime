library(tidyverse)
library(wtl)
library(tumopp)
library(ggtree)
refresh('rtumopp')
refresh('cellpoptime')

.ggtree = function(.tbl) {
  ggtree(.tbl) +
  geom_tree() +
  geom_tiplab() +
  wtl::theme_wtl()
}

N = 20L
.param = c('-D2', '-k100', paste0('-N', N), '-U10', '--mb=99', '--ms1mut', '-u0')
samples = purrr::rerun(4L, {
  (tumopp::mslike(N, .param) %>% wtl::parse_ms(byrow=TRUE))[[1L]]
}) %>% print()

# samples = wtl::ms(6L, 4L, theta = 100) %>% wtl::parse_ms() %>% print()
trees = samples %>% purrr::map(cellpoptime::infer_rooted_tree) %>% print()
trees %>%
  cellpoptime::as_multiphylo() %>%
  .ggtree() +
  facet_wrap(~.id)

.tree = trees[[1]]
.base = .tree %>% cellpoptime::remove_outgroup() %>% print()

.base %>% add_extra_columns()
.base %>% add_extra_columns() %>% filter_scale_tips()

p_binom(c(1, 8))
detect_driver(c(1, 8))

.base %>% add_extra_columns() %>% nest_tippairs() #%>% nest_tippairs() %>% nest_tippairs()

rescale_branches = function(x) {
  num_edges = nrow(x)
  while (nrow(x) > 1L) {
    x = nest_tippairs(x)
  }
  while (nrow(x) < num_edges) {
    x = unnest_children(x)
  }
  x %>% dplyr::arrange(node)
}

.shrunk = .base %>% add_extra_columns() %>% rescale_branches() %>% print()
.shrunk_cols = .shrunk %>% dplyr::select(-parent, -branch.length, -label, -is_tip)

.root = .base %>% dplyr::filter(is.na(.data$branch.length)) %>% dplyr::pull("parent") %>% print()
.graph = .base[,c(1L, 2L)] %>% as.matrix() %>% igraph::graph_from_edgelist()
.gdist = igraph::distances(.graph, as.character(.root), mode = "out") %>% as.vector()

.fortified = list(.base, .shrunk) %>%
  as_multiphylo() %>%
  ggtree::fortify() %>%
  dplyr::left_join(.shrunk_cols, by = "node") %>%
  dplyr::mutate(step = .gdist[.data$node]) %>%
  print()

.p_min = min(.fortified$p_driver, na.rm=TRUE)
filter_p_driver = function(x) {dplyr::filter(x, p_driver == .p_min)}

.p = ggplot(.fortified, aes(x, y)) +
  geom_tree() +
  geom_nodepoint(data = filter_p_driver, colour = "orangered", size = 4) +
  geom_tiplab() +
  facet_wrap(~.id, ncol=1) +
  wtl::theme_wtl() +
  wtl::erase(axis.title)
.p
ggsave('scaletree.pdf', .p, width=7, height=9.9)

.tree1 = .fortified %>% dplyr::filter(.id == "Tree #1") %>% print()
.tree2 = .fortified %>% dplyr::filter(.id == "Tree #2") %>% print()

.plot_sequencially = function(df) {
  .xlim = range(df$x)
  .ylim = range(df$y)
  purrr::map(seq_len(max(df$step)), ~{
    ggplot(df %>% dplyr::filter(step <= .x), aes(x, y)) +
      geom_tree() +
      geom_nodepoint(data = filter_p_driver, colour = "orangered", size = 4) +
      geom_tiplab() +
      coord_cartesian(xlim = .xlim, ylim = .ylim) +
      theme_void()
      # wtl::theme_wtl() + wtl::erase(axis.title)
  })
}
.plts = .plot_sequencially(.tree1)
.tmpl = "~/git/slides-draft/content/tumopp2018ncc/figure-ncc/tree1.%d.png"
purrr::iwalk(.plts, ~ggsave(sprintf(.tmpl, .y), .x, width=4, height=2))

system("magick -loop 1 -delay 50 tree*.png tree1-raw.gif")
system("magick -loop 1 tree1-raw.gif -layers Optimize tree1.gif")

.rec_rescale_branches = function(x) {
  num_edges = nrow(x)
  l = NULL
  while (nrow(x) > 1L) {
    x = nest_tippairs(x)
    l = c(l, list(x))
  }
  purrr::map(l, ~{
    while (nrow(.x) < num_edges) {
      .x = unnest_children(.x)
    }
    .x %>% dplyr::arrange(node)
  })
}
.shrunkl = .base %>% add_extra_columns() %>% .rec_rescale_branches() %>% print()

.fortified_l = .shrunkl %>% purrr::map(~{
  .x %>%
    as_phylo() %>%
    ggtree::fortify() %>%
    dplyr::left_join(.shrunk_cols, by = "node")
})

.plot_scale = function(.x, xlim = .xlim) {
  ggplot(.x, aes(x, y)) +
    geom_tree() +
    geom_nodepoint(data = filter_p_driver, colour = "orangered", size = 4) +
    geom_tiplab() +
    coord_cartesian(xlim = xlim) +
    theme_void()
    # wtl::theme_wtl() + wtl::erase(axis.title)
}

.xlim = range(.fortified$x) %>% print()
.plts_scale = .fortified_l %>% purrr::map(.plot_scale)
.tmpl = "~/git/slides-draft/content/tumopp2018ncc/figure-ncc/scale.%02d.png"
purrr::iwalk(.plts_scale, ~ggsave(sprintf(.tmpl, .y), .x, width=3, height=1, scale=2.4, dpi=240))
.final = .fortified_l %>% purrr::pluck(length(.)) %>% .plot_scale(NULL)
ggsave(sprintf(.tmpl, 99L), .final, width=3, height=1, scale=2.4, dpi=240)

system("magick -loop 1 -delay 50 scale*.png scale-raw.gif")
system("magick -loop 1 scale-raw.gif -layers Optimize scale.gif")
