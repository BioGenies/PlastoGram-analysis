do_baseline_cv <- function(ngram_matrix, data_dfs_cv) {
  lapply(1:length(data_dfs_cv), function(i) {
    data_df <- data_dfs_cv[[i]]
    lapply(unique(data_df[["fold"]]), function(ith_fold) {
      dat <- ngram_matrix[data_df[["fold"]] != ith_fold, ]
      test_dat <- ngram_matrix[data_df[["fold"]] == ith_fold, ]
      
      trained_model <- ranger(x = dat, y = as.factor(data_df[["dataset"]][data_df[["fold"]] != ith_fold]),
                              write.forest = TRUE, probability = TRUE, num.trees = 500, 
                              verbose = FALSE, seed = 427244)
      
      data_df %>% 
        filter(fold == ith_fold) %>% 
        select(seq_name, fold, dataset) %>% 
        cbind(., predict(trained_model, test_dat)[["predictions"]]) %>% 
        mutate(Localization = c("N_E", "N_S", "N_TL_SEC", "N_TL_TAT", "N_TM", "P_IM", "P_S", "P_TM")
               [max.col(.[, c("N_E", "N_S", "N_TL_SEC", "N_TL_TAT", "N_TM", "P_IM", "P_S", "P_TM")])],
               rep = i)
    }) %>% bind_rows()
  }) %>% bind_rows()
} 


train_baseline_model <- function(ngrams, data_df) {
  ranger(x = ngrams, y = as.factor(data_df[["dataset"]]),
         write.forest = TRUE, probability = TRUE, num.trees = 500, 
         verbose = FALSE, seed = 427244)
}


evaluate_baseline_model <- function(baseline_model, test_ngrams, test_dat, imp_baseline_ngrams) {
  
  test_dat %>% 
    add_missing_ngrams(baseline_model[["forest"]][["independent.variable.names"]]) %>% 
    select(seq_name, dataset) %>% 
    cbind(., predict(baseline_model, test_ngrams)[["predictions"]]) %>% 
    mutate(Localization = c("N_E", "N_S", "N_TL_SEC", "N_TL_TAT", "N_TM", "P_IM", "P_S", "P_TM")
           [max.col(.[, c("N_E", "N_S", "N_TL_SEC", "N_TL_TAT", "N_TM", "P_IM", "P_S", "P_TM")])])
}

