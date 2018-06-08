library(tidyverse)
library(wtl)
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
