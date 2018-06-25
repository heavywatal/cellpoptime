library(tidyverse)
library(wtl)
library(tumopp)
library(ggtree)
refresh('rtumopp')
refresh('cellpoptime')

.wd = "~/Dropbox/working/tumopp/cell-diversity-iwate/"
setwd(.wd)

.ggtree = function(.tbl) {
  ggtree(.tbl) +
  geom_tree() +
  geom_tiplab() +
  wtl::theme_wtl()
}

N = 20L
.param = c('-D2', '-k100', paste0('-N', N), paste0('-U', N / 2L), '--mb=99', '--ms1mut', '-u0')
samples = purrr::rerun(4L, {
  (tumopp::mslike(N, .param) %>% wtl::parse_ms(byrow=TRUE))[[1L]]
}) %>% print()

N = 127L
.param = c('-D2', '-k1000', paste0('-N', N), paste0('-U', N %/% 2L - 1L), '--mb=99', '--ms1mut', '-u0')
samples_full = purrr::rerun(4L, {
  (tumopp::mslike(N, .param) %>% wtl::parse_ms(byrow=TRUE))[[1L]]
})
wtl::save_as(samples_full)
# samples_full = readRDS("samples_full.rds")

samples_part = purrr::rerun(4L, {
  (tumopp::mslike(16L, .param) %>% wtl::parse_ms(byrow=TRUE))[[1L]]
})
wtl::save_as(samples_part)
# samples_part = readRDS("samples_part.rds")

samples = samples_full
samples = samples_part

# samples = wtl::ms(6L, 4L, theta = 100) %>% wtl::parse_ms() %>% print()
trees = samples %>% purrr::map(cellpoptime::infer_rooted_tree) %>% print()
trees %>%
  cellpoptime::as_multiphylo() %>%
  .ggtree() +
  facet_wrap(~.id)

.tree = trees[[1L]]
.base = .tree %>% cellpoptime::remove_outgroup() %>% print()

# refresh('cellpoptime')
.shrunk = .base %>% add_extra_columns() %>% scale_branches() %>% print()
.shrunk %>% dplyr::group_by(is_tip) %>% summarise(sum(exp_sibs))

.mutant = .shrunk %>% dplyr::top_n(1L, desc(p_driver)) %>% print()

.shrunk_cols = .shrunk %>%
  group_clade(.mutant$node) %>%
  dplyr::select(-parent, -branch.length, -label, -is_tip)

.root = .base %>% dplyr::filter(is.na(.data$branch.length)) %>% dplyr::pull("parent") %>% print()
.graph = .base[,c(1L, 2L)] %>% igraph::graph_from_data_frame()
igraph::degree(.graph, mode = "in")
.gdist = igraph::distances(.graph, as.character(.root), mode = "out") %>%
  as.vector() %>%
  setNames(names(igraph::V(.graph))) %>%
  print()

.fortified = list(.base, .shrunk) %>%
  as_multiphylo() %>%
  ggtree::fortify() %>%
  dplyr::left_join(.shrunk_cols, by = "node") %>%
  dplyr::mutate(step = .gdist[as.character(.data$node)]) %>%
  print()

.filter_driver_mutation = function(x) {
  x %>%
    tidyr::nest(-.id) %>%
    dplyr::transmute(.data$.id, data = purrr::map(data, ~{
      row = dplyr::top_n(.x, 1L, dplyr::desc(.data$p_driver))
      upstream_corner(.x, row$node)
    })) %>%
    tidyr::unnest()
}
.fortified %>% .filter_driver_mutation()

.fortified %>% .plot_snapshot(NULL) + facet_wrap(~.id, ncol = 1L)

.plot_snapshot = function(df, xlim, ylim = NULL, step = 65535) {
  df %>% dplyr::filter(.data$step <= !!step) %>%
    dplyr::mutate(isTip = ! node %in% parent) %>%
    ggplot(aes(x, y)) +
    geom_tree(aes(colour = group)) +
    ggtree::geom_tippoint(aes(colour = group), size = 2, alpha = 0.5) +
    scale_colour_manual(values = c(`0` = "#f768a1", `1` = "#7a0177"), guide = FALSE) +
    # geom_nodepoint(data = .driver_data, colour = "dodgerblue", size = 4) +
    # geom_tiplab() +
    # geom_text2(aes(subset=!isTip, label=node), hjust=-.3) +
    coord_cartesian(xlim = xlim, ylim = ylim) +
    theme_void()
    # wtl::theme_wtl() + wtl::erase(axis.title)
}

.ggsave = function(filename, plot) {
  ggsave(filename, plot, width=3, height=1, scale=2.4, dpi=200)
}

# #######1#########2#########3#########4

.tree1 = .fortified %>% dplyr::filter(.id == "Tree #1") %>% print()
.tree2 = .fortified %>% dplyr::filter(.id == "Tree #2") %>% print()

.tree1 %>% .plot_snapshot(NULL)
.tree1 %>% dplyr::filter(step <= 4) %>% .plot_snapshot(NULL)

# in cell-division scale
.plot_sequencially = function(df) {
  .xlim = range(df$x)
  .ylim = range(df$y)
  seq_len(max(df$step)) %>%
    purrr::map(.plot_snapshot, df = df, xlim = .xlim, ylim = .ylim)
}

.plts = .plot_sequencially(.tree1)
.tmpl_g = "scale-genetic-%02d.pdf"
# .tmpl_g = "part-scale-genetic-%02d.pdf"
purrr::iwalk(.plts, ~.ggsave(sprintf(.tmpl_g, .y), .x))

system("magick -loop 1 -delay 50 tree*.png tree1-raw.gif")
system("magick -loop 1 tree1-raw.gif -layers Optimize tree1.gif")

# #######1#########2#########3#########4

.shrunkl = .base %>% add_extra_columns() %>% scale_branches_record() %>% print()

.fortified_l = .shrunkl %>% purrr::map(~{
  .x %>%
    as_phylo() %>%
    ggtree::fortify() %>%
    dplyr::left_join(.shrunk_cols, by = "node") %>%
    dplyr::mutate(step = .gdist[as.character(.data$node)])
})

.xlim = range(.fortified$x) %>% print()
.plts_scale = .fortified_l %>% purrr::map(.plot_snapshot, xlim = .xlim)
.tmpl_t = "scale-time-%02d.pdf"
# .tmpl_t = "part-scale-time-%02d.pdf"
purrr::iwalk(.plts_scale, ~.ggsave(sprintf(.tmpl_t, .y), .x))
.p_final = .fortified_l %>% purrr::pluck(length(.)) %>% .plot_snapshot(xlim = NULL)
.ggsave(sprintf(.tmpl_t, 90L), .p_final)

system("magick -loop 1 -delay 50 scale*.png scale-raw.gif")
system("magick -loop 1 scale-raw.gif -layers Optimize scale.gif")

# #######1#########2#########3#########4#########5#########6#########7#########
# Triangle

refresh("cellpoptime")

.final_tbl = .fortified_l %>% purrr::pluck(length(.)) %>% print()

.fastnodes = tumopp:::paths_to_sink(.graph, as.character(.mutant$node))[[1L]] %>% as.integer() %>% print()
.fastvertices = .final_tbl %>%
  dplyr::filter(isTip, node %in% .fastnodes) %>%
  dplyr::filter(y %in% range(y)) %>%
  dplyr::pull(node) %>%
  print()

.joint = .final_tbl %>% upstream_corner(.mutant$node) %>% print()

.fast_tbl = .final_tbl %>%
  dplyr::filter(node %in% .fastvertices) %>%
  dplyr::mutate(y = y + ifelse(y < mean(y), -0.5, 0.5)) %>%
  dplyr::bind_rows(.joint) %>%
  print()

.slowvertices = .final_tbl %>%
  dplyr::filter(isTip, !node %in% .fastnodes) %>%
  dplyr::filter(y %in% range(y)) %>%
  dplyr::pull(node) %>%
  print()

.mrca = .final_tbl %>%
  dplyr::filter(parent == node) %>%
  # dplyr::mutate(x = -0.5 * max(.final_tbl$x)) %>%
  print()

.slow_tbl = .final_tbl %>%
  dplyr::filter(.data$node %in% .slowvertices) %>%
  dplyr::mutate(y = y + ifelse(y < mean(y), -0.5, 0.5)) %>%
  dplyr::arrange(desc(y)) %>%
  dplyr::bind_rows(.mrca, .joint) %>%
  print()

.p_triangle = .final_tbl %>%
  .plot_snapshot(xlim = NULL) +
  # geom_point(aes(size = exp_sibs), alpha = 0.4) +
  geom_polygon(data = .fast_tbl, alpha = 0.6, fill = "#7a0177")+
  geom_polygon(data = .slow_tbl, alpha = 0.6, fill = "#f768a1")
.p_triangle
.ggsave(sprintf(.tmpl_t, 99L), .p_triangle)
