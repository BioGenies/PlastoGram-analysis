library(biogram)
library(dplyr)
library(readxl)

source("./functions/do_cdhit.R")

data_path <- "~/Dropbox/Projekty/BioNgramProjects/PlastoGram/"

annotation_df <- read_xlsx(paste0(data_path, "Annotations/All_proteins.xlsx"))

# Combine all sequences downloaded from UniProt
all_seqs <- unlist(lapply(list.files(paste0(data_path, "Annotations/UniProt_Raw_data"),
                                     pattern = ".fasta", full.names = TRUE),
                          function(ith_file) {
                            read_fasta(ith_file)
                          }), recursive = FALSE)
names(all_seqs) <- sapply(names(all_seqs), function(i) strsplit(i, "|", fixed = TRUE)[[1]][2])
write_fasta(all_seqs, paste0(data_path, "Sequences/All_sequences.fa"))


# Extract sequences corresponding to each dataset and count numbers of sequences
seq_numbers <- lapply(c("N_OM", "N_IM", "P_IM", "N_S", "P_S", "N_TM", "P_TM", "N_TL_SEC", "N_TL_TAT", "P_TL_SEC"), function(ith_set) {
  accessions <- filter(annotation_df, Final_dataset == ith_set)[["Entry"]]
  seqs <- all_seqs[which(names(all_seqs) %in% accessions)]
  write_fasta(seqs, paste0(data_path, "Sequences/", ith_set, ".fa"))
  data.frame(Dataset = ith_set,
             Sequences = length(seqs))
}) %>% bind_rows()


# Check numbers of sequences after CD-HIT reduction with different cut-offs
seq_numbers_cdhit <- lapply(c("N_OM", "N_IM", "P_IM", "N_S", "P_S", "N_TM", "P_TM", "N_TL_SEC", "N_TL_TAT", "P_TL_SEC"), function(ith_group) {
  seqs <- read_fasta(paste0(data_path, "Sequences/", ith_group, ".fa"))
  filtered_1 <- filter_with_cdhit(seqs, 1)
  filtered_0.95 <- filter_with_cdhit(seqs, 0.95)
  filtered_0.9 <- filter_with_cdhit(seqs, 0.9)
  data.frame(Dataset = ith_group,
             all = length(seqs),
             cdhit_1 = length(filtered_1),
             cdhit_0.95 = length(filtered_0.95),
             cdhit_0.9 = length(filtered_0.9))
}) %>% bind_rows()
