library(dplyr)
library(biogram)
library(readxl)
library(tidyr)
library(ggplot2)

plot_tax_rep <- function(dat, set) {
  p <- ggplot(dat, aes(x = Organism, y = count)) +
    geom_col(color = "black", size = 0.3, fill = "#7fba70") +
    geom_text(aes(x = Organism, y = count, label = count),
              position = position_dodge(width = 1), hjust = -0.25, size = 3) +
    coord_flip() +
    theme_bw() +
    ggtitle(set) 
  ggsave(filename = paste0("plots/Tax_rep_", set, ".png"), plot = p, width = 10, height = nrow(dat)/4)
}


ntm <- read_fasta("N_TM.fa")

dplg <- ntm[which(sapply(ntm, function(i) grepl("DPLG", paste(i, collapse = ""))))]

write_fasta(dplg, file = "dplg.fa")


df <- read_xlsx("/media/kasia/Data/Dropbox/Projekty/BioNgramProjects/PlastoGram/Annotations/All_proteins.xlsx")
ntl <- read_fasta("N_TL_TAT.fa")

tax_rep <- df %>% 
  filter(!(Final_dataset %in% c("-", "P_TL_SEC"))) %>% 
  group_by(Final_dataset, Organism) %>% 
  summarise(count = n()) 
tax_rep_df <- tax_rep %>% 
  pivot_wider(names_from = Final_dataset, values_from = count)

tax_rep %>% 
  group_by(Final_dataset) %>% 
  summarise(n_org = n())

tax_rep_to_plot <- df %>% 
  filter(!(Final_dataset %in% c("-", "P_TL_SEC"))) %>% 
  group_by(Final_dataset, Organism) %>% 
  summarise(count = n()) %>% 
  mutate(Organism = ifelse(count == 1, "Other (organisms with 1 protein)", Organism)) %>% 
  group_by(Final_dataset, Organism) %>% 
  summarise(count = ifelse(n() == 1, count, n()))
  

lapply(unique(filter(tax_rep_to_plot, Final_dataset %in% c("N_OM", "N_TM", "N_S", "N_TL_SEC", "N_TL_TAT"))[["Final_dataset"]]), function(ith_set) {
  filter(tax_rep_to_plot, Final_dataset == ith_set) %>% 
  plot_tax_rep(ith_set)
})

plot_tax_rep(filter(tax_rep, Final_dataset == "N_IM"), "N_IM")


tax_rep_to_plot_p <- df %>% 
  filter(!(Final_dataset %in% c("-", "P_TL_SEC"))) %>% 
  group_by(Final_dataset, Organism) %>% 
  summarise(count = n()) %>% 
  mutate(Organism = case_when(count == 1 ~ "Organisms with 1 protein", 
                              count == 2 ~ "Organisms with 2 proteins",
                              count > 2 ~ Organism)) %>% 
  group_by(Final_dataset, Organism) %>% 
  summarise(count = ifelse(n() == 1, count, n()))


lapply(unique(filter(tax_rep_to_plot, Final_dataset %in% c("P_TM", "P_S", "P_IM"))[["Final_dataset"]]), function(ith_set) {
  filter(tax_rep_to_plot_p, Final_dataset == ith_set) %>% 
    plot_tax_rep(ith_set)
})
