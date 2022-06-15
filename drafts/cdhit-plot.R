library(ggplot2)
library(dplyr)
library(tidyr)

png("cdhit-sizes.png", width = 680*4, height = 680*4, res = 190)
read.csv("data/CD-HIT_reduction_results.csv") %>% 
  pivot_longer(cols = -Dataset) %>% 
  mutate(name = paste0(as.numeric(gsub("X", "", name)) * 100, "%")) %>% 
  ggplot(aes(x = name, y  = value, label = value)) +
  geom_col() +
  geom_label(aes(y = value/2)) +
  facet_wrap(~ Dataset, scales = "free_y") +
  scale_x_discrete("CD-Hit reduction") +
  scale_y_continuous("Number of sequences") +
  theme_bw()
dev.off()
