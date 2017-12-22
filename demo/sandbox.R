library(tidyverse)
library(wtl)
library(ape)
library(ggtree)

results = ms(6L, 5L, theta=100) %>% print()

.m = results[[1]] %>% as_int_matrix() %>% add_outgroup()
.d = .m %>% dist(method='manhattan') %>% print()
.tree = .d %>% ape::fastme.bal()
.rooted = .tree %>% ape::root('0')

attributes(.rooted)
str(.rooted)
.rooted$edge
.rooted$edge.length
.rooted$tip.label
.rooted$Nnode

.rooted %>% plot()

ggplot(.rooted, aes(x, y))+
geom_tree()+
geom_tiplab()+
theme_tree()
