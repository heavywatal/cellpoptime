library(tidyverse)
library(wtl)

make_data = function(tmin, tmax, p=0.4, sd=0.3) {
    tmean = mean(c(tmax, tmin))
    tibble(x= seq(tmin, tmax, length.out=20),
           ymin= p - p * pnorm(x, tmean, sd),
           ymax= p + (1 - p) * pnorm(x, tmean, sd)) %>%
    add_row(x=10, ymin=0, ymax=1)
}
make_data(2, 4)

.tbl = tibble(tmin=c(0, 3.0), tmax=c(5.5, 5.9), p=c(0.5, 0.8), sd=c(0.55, 0.30)) %>%
    purrr::pmap_dfr(make_data, .id='id') %>%
    print()

.p = ggplot(.tbl, aes(x))+
    geom_ribbon(aes(ymin = ymin, ymax = ymax, fill = id), alpha = 1)+
    scale_fill_manual(values = c("#f768a1", "#7a0177"), guide = FALSE)+
    coord_cartesian(xlim = c(0.5, 4.5), ylim = c(0, 1), expand = FALSE)+
    # labs(x='Time')+
    wtl::erase(axis.title, axis.text, axis.ticks, panel.grid)
.p

# .outfile = '~/git/slides-draft/content/tumopp2018ncc/figure-ncc/muller.png'
.outfile = '~/Dropbox/working/tumopp/cell-diversity-iwate/muller.pdf'
ggsave(.outfile, .p, width=3, height=1, scale=1.2, dpi=250)

# #######1#########2#########3#########4#########5#########6#########7#########

library(ggmuller)

edges = tibble(Parent = c("A", "A"), Identity = c("B", "C")) %>% print()
pop = tibble(Identity = c("A", "B", "C"), lambda = c(1.05, 1.2, 1.5), start = c(0, 30, 50)) %>%
  mutate(data = purrr::map2(lambda, start, ~{
    tibble(Generation = seq(0, 100, 2), Population = 10 * .x ** (Generation - .y))
  })) %>%
  tidyr::unnest() %>%
  print()

Muller_df = get_Muller_df(edges, pop) %>% as_tibble() %>% print()
Muller_plot(Muller_df, add_legend = TRUE, palette = "Blues")
