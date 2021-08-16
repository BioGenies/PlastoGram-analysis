create_datasets <- function(sequence_file, annotation_file, output_dir) {
  seqs <- read_fasta(sequence_file)
  dat <- read_xlsx(annotation_file)
  datasets <- c("N_OM", "N_IM", "P_IM", "N_S", "P_S", "N_TM", "P_TM", "N_TL_SEC", "N_TL_TAT")
  
  lapply(datasets, function(ith_dataset) {
    selected <- unique(filter(dat, Final_dataset == ith_dataset)[["Entry"]])
    selected_seqs <- seqs[which(names(seqs) %in% selected)]
    write_fasta(selected_seqs, paste0(output_dir, ith_dataset, ".fa"))
    selected_seqs
  }) %>% setNames(datasets)
}
