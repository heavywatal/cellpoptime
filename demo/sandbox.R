library(tidyverse)
library(wtl)
library(ggtree)
loadNamespace('cowplot')

samples = ms(6L, 4L, theta=100) %>% print()

trees = samples %>% purrr::map(infer_rooted_tree) %>% print()

.plot = function(.tbl) {
    as.phylo(.tbl) %>%
    ggplot(aes(x, y))+
    geom_tree()+
    geom_tiplab()+
    theme_bw()
}

.plts = purrr::map(trees, .plot)

cowplot::plot_grid(plotlist=.plts, ncol=2L)


remove_dummy_root = function(x) {
    if ('0' %in% x$label) {
        dplyr::filter(x, label != '0' | is.na(label)) %>%
        dplyr::mutate(parent= parent - 1L, node= node - 1L)
    } else {x}
}

.base = trees[[1]] %>% remove_dummy_root() %>% print()

.plot(.base)

filter_tippairs = function(x) {
    dplyr::group_by(x, parent) %>%
    dplyr::filter(all(isTip)) %>%
    dplyr::mutate(branch.length= min(branch.length))
    #TODO: min is not enough
}

filter_remaining = function(x, new_tips) {
    dplyr::group_by(x, parent) %>%
    dplyr::filter(!all(isTip)) %>%
    dplyr::mutate(isTip= isTip | (node %in% new_tips))
}

shrink = function(x) {
    tipinfo = dplyr::select(x, node, isTip)
    done = dplyr::filter(x, is.na(branch.length))
    x = dplyr::filter(x, !is.na(branch.length))
    while (nrow(x) > 1L) {
        done = bind_rows(done, filter_tippairs(x)) %>% print()
        x = filter_remaining(x, done$parent) %>% print()
    }
    bind_rows(done, x) %>%
    dplyr::select(-isTip) %>%
    dplyr::left_join(tipinfo)
}

.mod = shrink(.base) %>% print()
.mod %>% .plot()
