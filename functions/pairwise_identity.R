calculate_pairwise_identity <- function(traintest, benchmark, sequences_cv_graphpart, sequences_independent_graphpart,
                                        envelope_target_df, target_df_graphpart, data_path) {
  
  sequence_datasets <- list("Holdout" = list("traintest" = traintest, 
                                             "independent" = benchmark),
                            "Partitioning" = list("traintest" = sequences_cv_graphpart, 
                                                  "independent" = sequences_independent_graphpart))
  
  lapply(names(sequence_datasets), function(ith_set) {
    print(paste0(ith_set))
    lapply(unique(envelope_target_df[["dataset"]]), function(ith_class) {
      print(paste0(ith_class))
      if(ith_set == "Holdout") {
        dat <- envelope_target_df
        t_seq_names <- filter(dat, seq_name %in% names(sequence_datasets[[ith_set]][["traintest"]]), dataset == ith_class)[["seq_name"]]
        i_seq_names <- filter(dat, seq_name %in% names(sequence_datasets[[ith_set]][["independent"]]), dataset == ith_class)[["seq_name"]]
      } else {
        dat <- target_df_graphpart
        t_seq_names <- filter(dat, cluster != 1, dataset == ith_class)[["seq_name"]]
        i_seq_names <- filter(dat, cluster == 1, dataset == ith_class)[["seq_name"]]
      }
      t_seqs <- sequence_datasets[[ith_set]][["traintest"]][which(names(sequence_datasets[[ith_set]][["traintest"]]) %in% t_seq_names)]
      i_seqs <- sequence_datasets[[ith_set]][["independent"]][which(names(sequence_datasets[[ith_set]][["independent"]]) %in% i_seq_names)]
      write_fasta(t_seqs, "b.fa")
      needle_res <- lapply(1:length(i_seqs), function(i) {
        write_fasta(i_seqs[i], "a.fa")
        system(paste0("needle -asequence a.fa -bsequence b.fa -gapopen 10 -gapextend 0.5 -outfile res.needle -sprotein"))
        res <- read_needle_res("res.needle")
        file.remove(c("a.fa", "res.needle"))
        res
      }) %>% bind_rows()
      file.remove("b.fa")
      write.csv(needle_res, paste0(data_path, "Pairwise_identity/", ith_set, "_", ith_class, ".csv"), row.names = FALSE)
      mutate(needle_res, 
             Dataset = ith_class,
             Version = ith_set)
    }) %>% bind_rows()
  }) %>% bind_rows()
}

read_needle_res <- function(needle_res_file) {
  x <- readLines(needle_res_file)
  first <- sapply(x[grepl("# 1:", x)], function(i) gsub("# 1: ", "", i), USE.NAMES = FALSE)
  second <- sapply(x[grepl("# 2:", x)], function(i) gsub("# 2: ", "", i), USE.NAMES = FALSE)
  ident <- sapply(x[grepl("Identity:", x)], function(i) gsub("^ |%)", "", strsplit(i, "(", fixed = TRUE)[[1]][2]), USE.NAMES = FALSE)
  data.frame(independent = first,
             traintest = second,
             identity = ident) 
}

calculate_mean_pairwise_identity <- function(pairwise_identity_df, data_path) {
  res <- pairwise_identity_df %>% 
    group_by(Version, Dataset) %>% 
    summarise(mean = mean(as.numeric(identity)), 
              median = median(as.numeric(identity)),
              max = max(as.numeric(identity)),
              min = min(as.numeric(identity))) 
  write.csv(filter(res, Version == "Holdout"), paste0(data_path, "mean_pairwise_identity_holdout.csv"), row.names = FALSE)
  write.csv(filter(res, Version == "Partitioning"), paste0(data_path, "mean_pairwise_identity_partitioning.csv"), row.names = FALSE)
  res
}