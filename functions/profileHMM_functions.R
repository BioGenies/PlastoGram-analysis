train_profileHMM <- function(train_seqs, model_name, remove_files = TRUE) {
  write_fasta(train_seqs, "train_seqs.fa")
  system('mafft  --quiet  --localpair  --maxiterate 16 --reorder "train_seqs.fa" > "train_seqs_aln.fa"')
  system(paste0("hmmbuild ", model_name, ".hmm train_seqs_aln.fa"))
  if(remove_files == TRUE) {
    file.remove(c("train_seqs.fa", "train_seqs_aln.fa"))
  }
}

predict_profileHMM <- function(test_seqs, model_name, remove_files = TRUE) {
  invisible({
  write_fasta(test_seqs, "test_seqs.fa")
  system(paste0("hmmsearch --tblout hmmer_res ", model_name, ".hmm test_seqs.fa")) 
  })
  res <- bind_rows(read_tblout("hmmer_res")) 
  if(remove_files == TRUE) {
    file.remove(c("test_seqs.fa", "hmmer_res", paste0(model_name, ".hmm")))
  }
  res 
}

combine_tat_sec_results <- function(tat_res, sec_res, fold_df_tat, fold_df_sec, fold) {
 tat_res <- tat_res %>% 
    mutate(Tat = (2^sequence_score) / (1+2^sequence_score),
           seq_name = domain_name) %>% 
    select(c("seq_name", "Tat"))
  sec_res <- sec_res %>% 
    mutate(Sec = (2^sequence_score) / (1+2^sequence_score),
           seq_name = domain_name) %>% 
    select(c("seq_name", "Sec"))
  all_res <- filter(bind_rows(fold_df_sec, fold_df_tat), fold == fold) %>% 
    left_join(full_join(sec_res, tat_res, by = "seq_name"), by = "seq_name")
  all_res[is.na(all_res)] <- 0
  all_res
}

do_hmmer_cv <- function(tat_seqs, sec_seqs) {
  folded_tat <- cvFolds(length(tat_seqs), K = 5)
  fold_df_tat <- data.frame(seq_name = names(tat_seqs)[folded_tat[["subsets"]]], 
                            fold = folded_tat[["which"]],
                            target = "Tat")
  folded_sec <- cvFolds(length(sec_seqs), K = 5)
  fold_df_sec <- data.frame(seq_name = names(sec_seqs)[folded_sec[["subsets"]]], 
                            fold = folded_sec[["which"]],
                            target = "Sec")
  
  lapply(1:5, function(ith_fold) {
    tat_train <- tat_seqs[filter(fold_df_tat, fold != ith_fold)[["seq_name"]]]
    tat_test <- tat_seqs[filter(fold_df_tat, fold == ith_fold)[["seq_name"]]]
    sec_train <- sec_seqs[filter(fold_df_sec, fold != ith_fold)[["seq_name"]]]
    sec_test <- sec_seqs[filter(fold_df_sec, fold == ith_fold)[["seq_name"]]]
    
    train_profileHMM(tat_train, "tat_model")
    train_profileHMM(sec_train, "sec_model")
    
    tat_res <- predict_profileHMM(c(tat_test, sec_test), "tat_model")
    sec_res <- predict_profileHMM(c(sec_test, tat_test), "sec_model")
    
    combine_tat_sec_results(tat_res, sec_res, 
                            filter(fold_df_tat, fold == ith_fold), 
                            filter(fold_df_sec, fold == ith_fold), 
                            ith_fold)
  }) %>% bind_rows()
}

summarise_hmmer_results <- function(cv_res, type) {
  cv_res %>% 
    select(c("seq_name", "fold", "target", type)) %>% 
    setNames(c("seq_name", "fold", "target", "prob")) %>% 
    mutate(target = ifelse(target == type, TRUE, FALSE),
           pred = ifelse(prob > 0.5, TRUE, FALSE)) %>% 
    get_cv_res_summary("TRUE", ngrams = FALSE)
}
