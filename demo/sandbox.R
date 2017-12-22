library(tidyverse)
library(wtl)
library(ggtree)
loadNamespace('cowplot')

samples = ms(6L, 4L, theta=100) %>% print()

trees = samples %>% purrr::map(infer_rooted_tree) %>% print()

.plts = purrr::map(trees, ~{
    as.phylo(.x) %>%
    ggplot(aes(x, y))+
    geom_tree()+
    geom_tiplab()+
    theme_bw()
})

cowplot::plot_grid(plotlist=.plts, ncol=2L)
