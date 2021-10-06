library(dplyr)
library(ggplot2)
library(tidyr)

dat <- read.csv("./drafts/Architectures_mean_performance.csv")

arch_dat <- mutate(dat, architecture = sapply(strsplit(model, "-"), first),
                   filtering = grepl(pattern = "Filtering", x = model, ignore.case = FALSE)) %>% 
  mutate(architecture = sapply(strsplit(architecture, "_"), function(i) i[2])) %>% 
  select(architecture, filtering, model, mean_AU1U, mean_N_IM_sensitivity, mean_N_OM_sensitivity, mean_N_TL_SEC_sensitivity) %>% 
  pivot_longer(cols = !c(architecture, model, filtering)) 


filter(arch_dat, name != "mean_AU1U") %>% 
  ggplot(aes(x = filtering, y = value, color = name)) +
  geom_boxplot() +
  facet_wrap(~ architecture, nrow = 1)


filter(arch_dat, name != "mean_AU1U") %>% 
  ggplot(aes(x = architecture, y = value, color = name)) +
  geom_boxplot() +
  facet_wrap(~ filtering, nrow = 1)
