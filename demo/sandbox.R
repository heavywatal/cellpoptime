library(tidyverse)
library(wtl)
library(ggtree)
wtl::refresh('rtumopp')
wtl::refresh('cellpoptime')

.wd = "~/Dropbox/working/tumopp/cell-diversity-iwate/"
setwd(.wd)

N = 20L
.param = c('-D2', '-k100', paste0('-N', N), paste0('-U', N / 2L), '--mb=99')
samples = purrr::rerun(4L, {
  tumopp::tumopp(.param)$graph[[1L]] %>% tumopp::make_sample(N, -1)
}) %>% print()

N = 127L
.param = c('-D2', '-k1000', paste0('-N', N), paste0('-U', N %/% 2L - 1L), '--mb=99')
samples_full = purrr::rerun(4L, {
  tumopp::tumopp(.param)$graph[[1L]] %>% tumopp::make_sample(N, -1)
})
wtl::save_as(samples_full)
# samples_full = readRDS("samples_full.rds")

samples_part = purrr::rerun(4L, {
  tumopp::tumopp(.param)$graph[[1L]] %>% tumopp::make_sample(16L, -1)
})
wtl::save_as(samples_part)
# samples_part = readRDS("samples_part.rds")

samples = samples_full
samples = samples_part

# samples = wtl::ms(6L, 4L, theta = 100) %>% wtl::parse_ms() %>% print()
trees = samples %>% purrr::map(cellpoptime::infer_rooted_tree) %>% print()
trees %>%
  cellpoptime::as_multiphylo() %>%
  ggtree() +
  geom_tree() +
  geom_tiplab() +
  wtl::theme_wtl() +
  facet_wrap(~.id)

.tree = trees[[1L]]
.base = .tree %>% cellpoptime::remove_outgroup() %>% print()

.scaled_rec = .base %>% scale_branches_record()
.shrunk = .scaled_rec %>% dplyr::last() %>% print()
# .shrunk = .base %>% scale_branches() %>% print()
.shrunk %>% dplyr::group_by(is_tip) %>% summarise(sum(exp_sibs))

.mutant = .shrunk %>% dplyr::top_n(1L, dplyr::desc(p_driver)) %>% print()
.meta_info = .shrunk %>%
  group_clade(.mutant$node) %>%
  dplyr::select(-parent, -branch.length, -label, -is_tip) %>%
  print()

.fortified = .scaled_rec %>%
  as_multiphylo() %>%
  ggtree::fortify() %>%
  dplyr::left_join(.meta_info, by = "node") %>%
  print()

.fortified %>% .plot_snapshot() + facet_wrap(~.id, ncol = 1L)

.plot_snapshot = function(df, xlim = NULL, ylim = NULL, xmax = 65535) {
  df %>% dplyr::filter(.data$x <= xmax) %>%
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
# Draw scaling steps

.xlim = range(.fortified$x) %>% print()
.nested = .fortified %>% tidyr::nest(-.id) %>% print()
.plts_scale = .nested$data %>% purrr::map(.plot_snapshot, xlim = .xlim)
.tmpl_t = "scale-time-%02d.pdf"
# .tmpl_t = "part-scale-time-%02d.pdf"
purrr::iwalk(.plts_scale, ~.ggsave(sprintf(.tmpl_t, .y), .x))
.p_final = .nested$data %>% purrr::pluck(length(.)) %>% .plot_snapshot()
.ggsave(sprintf(.tmpl_t, 90L), .p_final)

system("magick -loop 1 -delay 50 scale*.png scale-raw.gif")
system("magick -loop 1 scale-raw.gif -layers Optimize scale.gif")

# #######1#########2#########3#########4
# Draw partial tree from small x

.fortified1 = .nested$data[[1L]] %>% print()
.fortified1 %>% .plot_snapshot(NULL)
.fortified1 %>% dplyr::filter(x <= 3) %>% .plot_snapshot()

.plot_sequencially = function(df) {
  .xlim = range(df$x)
  .ylim = range(df$y)
  unique(df$x) %>% sort() %>% tail(-1L) %>%
    purrr::map(.plot_snapshot, df = df, xlim = .xlim, ylim = .ylim)
}

.plts = .plot_sequencially(.fortified1)
.tmpl_g = "scale-genetic-%02d.pdf"
# .tmpl_g = "part-scale-genetic-%02d.pdf"
purrr::iwalk(.plts, ~.ggsave(sprintf(.tmpl_g, .y), .x))

system("magick -loop 1 -delay 50 tree*.png tree1-raw.gif")
system("magick -loop 1 tree1-raw.gif -layers Optimize tree1.gif")

# #######1#########2#########3#########4#########5#########6#########7#########
# Triangle

refresh("cellpoptime")

.fortified_last = .nested$data %>% purrr::pluck(length(.)) %>% print()

.fastvertices = .fortified_last %>%
  dplyr::filter(isTip, group == 1) %>%
  dplyr::filter(y %in% range(y)) %>%
  dplyr::pull(node) %>%
  print()

.joint = .fortified_last %>% upstream_corner(.mutant$node) %>% print()

.fast_tbl = .fortified_last %>%
  dplyr::filter(node %in% .fastvertices) %>%
  dplyr::mutate(y = y + ifelse(y < mean(y), -0.5, 0.5)) %>%
  dplyr::bind_rows(.joint) %>%
  print()

.slowvertices = .fortified_last %>%
  dplyr::filter(isTip, group == 0) %>%
  dplyr::filter(y %in% range(y)) %>%
  dplyr::pull(node) %>%
  print()

.mrca = .fortified_last %>%
  dplyr::filter(parent == node) %>%
  # dplyr::mutate(x = -0.5 * max(.fortified_last$x)) %>%
  print()

.slow_tbl = .fortified_last %>%
  dplyr::filter(.data$node %in% .slowvertices) %>%
  dplyr::mutate(y = y + ifelse(y < mean(y), -0.5, 0.5)) %>%
  dplyr::arrange(desc(y)) %>%
  dplyr::bind_rows(.mrca, .joint) %>%
  print()

.p_triangle = .fortified_last %>%
  .plot_snapshot(xlim = NULL) +
  # geom_point(aes(size = exp_sibs), alpha = 0.4) +
  geom_polygon(data = .fast_tbl, alpha = 0.6, fill = "#7a0177")+
  geom_polygon(data = .slow_tbl, alpha = 0.6, fill = "#f768a1")
.p_triangle
.ggsave(sprintf(.tmpl_t, 99L), .p_triangle)
