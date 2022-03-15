create_datasets <- function(sequence_file, annotation_file, output_dir, datasets) {
  seqs <- read_fasta(sequence_file)
  dat <- read_xlsx(annotation_file)
  
  lapply(datasets, function(ith_dataset) {
    selected <- unique(filter(dat, Final_dataset == ith_dataset)[["Entry"]])
    selected_seqs <- seqs[which(names(seqs) %in% selected)]
    write_fasta(selected_seqs, paste0(output_dir, ith_dataset, ".fa"))
    selected_seqs
  }) %>% setNames(datasets)
}


process_for_graphpart <- function(sequence_dataset_list) {
  labelled_seqs <- lapply(names(sequence_dataset_list), function(ith_set) {
    seqs <- sequence_dataset_list[[ith_set]]
    names(seqs) <- paste0(names(seqs), "|label=", ith_set)
    seqs
  }) %>% unlist(recursive = FALSE)
  labelled_seqs
}

run_graphpart <- function(graphpart_input, sequence_dir, n_threads = 24) {
  write_fasta(graphpart_input, paste0(sequence_dir, "Datasets_for_graph-part_cdhit_0.9.fa"))
  system(paste0("graphpart needle --fasta-file ", sequence_dir, "Datasets_for_graph-part_cdhit_0.9.fa --threshold 0.4 --out-file ", 
                sequence_dir, "Graph-part_results.csv --labels-name label --test-ratio 0 --val-ratio 0.15 --threads 24 --no-moving"))
  label_df <- data.frame(dataset = unique(sapply(names(graphpart_input), function(i) gsub("label=", "", strsplit(i, "|", fixed = TRUE)[[1]][2]))),
                         `label.val` = 0:8)
  res <- read.csv(paste0(sequence_dir, "Graph-part_results.csv")) 
  left_join(res, label_df)
}

