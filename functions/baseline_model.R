do_baseline_cv <- function(ngram_matrix, data_df) {
  lapply(unique(data_df[["fold"]]), function(ith_fold) {
    dat <- ngram_matrix[data_df[["fold"]] != ith_fold, ]
    test_dat <- ngram_matrix[data_df[["fold"]] == ith_fold, ]
    
    trained_model <- ranger(x = dat, y = as.factor(data_df[["dataset"]][data_df[["fold"]] != ith_fold]),
                            write.forest = TRUE, probability = TRUE, num.trees = 500, 
                            verbose = FALSE, seed = 108567)
    
    data_df %>% 
      filter(fold == ith_fold) %>% 
      select(seq_name, fold, dataset) %>% 
      cbind(., predict(trained_model, test_dat)[["predictions"]]) %>% 
      mutate(pred = c("N_IM", "N_OM", "N_S", "N_TL_SEC", "N_TL_TAT", "N_TM", "P_IM", "P_S", "P_TM")
             [max.col(.[, c("N_IM", "N_OM", "N_S", "N_TL_SEC", "N_TL_TAT", "N_TM", "P_IM", "P_S", "P_TM")])])
    
  }) %>% bind_rows()
} 