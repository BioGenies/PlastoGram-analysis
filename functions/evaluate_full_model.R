get_all_imp_ngrams <- function(ngram_models) {
  sapply(names(ngram_models), function(i) ngram_models[[i]][["forest"]][["independent.variable.names"]]) %>% 
    unname() %>% 
    unlist() %>% 
    unique()
}

add_missing_ngrams <- function(ngrams, imp_ngrams) {
  missing <- imp_ngrams[which(!(imp_ngrams %in% colnames(ngrams)))]
  if(length(missing) > 0) {
    cbind(ngrams, setNames(lapply(missing, function(x) x = 0), missing))
  } else {
    ngrams
  }
}



predict_with_PlastoGram <- function(ngram_models, hmm_models, higher_level_model, test_ngrams, test_sequences, test_df, informative_ngrams) {
  test_ngrams <- add_missing_ngrams(test_ngrams, informative_ngrams)
  ngram_preds_list <- lapply(names(ngram_models), function(ith_model) {
    df <- data.frame(seq_name = test_df[["seq_name"]],
                     pred = predict(ngram_models[[ith_model]], test_ngrams)[["predictions"]][, "TRUE"])
    colnames(df) <- c("seq_name", ith_model)
    df
  }) 
  ngram_models_res <- reduce(ngram_preds_list, full_join, by = "seq_name")
  
  tat_res <- predict_profileHMM(test_sequences, hmm_models[["Tat_model"]], remove_files = FALSE) 
  sec_res <- predict_profileHMM(test_sequences, hmm_models[["Sec_model"]], remove_files = FALSE) 
  
  tat_res <- tat_res %>% 
    mutate(Tat_model = (2^sequence_score) / (1+2^sequence_score),
           seq_name = domain_name) %>% 
    select(c("seq_name", "Tat_model"))
  sec_res <- sec_res %>% 
    mutate(Sec_model = (2^sequence_score) / (1+2^sequence_score),
           seq_name = domain_name) %>% 
    select(c("seq_name", "Sec_model"))
  hmm_models_res <- full_join(sec_res, tat_res, by = "seq_name")
  
  all_res <- left_join(ngram_models_res, hmm_models_res, by = "seq_name")
  all_res[is.na(all_res)] <- 0
  
  glm_preds <- data.frame(seq_name = test_df[["seq_name"]],
                          dataset = test_df[["dataset"]],
                          predict(higher_level_model, all_res, type = "probs"),
                          Localization = predict(higher_level_model, all_res))
  
  list("Lower-order_models_preds" = all_res,
       "Higher-order_model_preds" = glm_preds,
       "Final_results" = data.frame(glm_preds[, c("seq_name", "dataset", "Localization")],
                                    Probability = sapply(1:nrow(glm_preds[, 3:(ncol(glm_preds)-1)]), function(i) max(glm_preds[i, 3:(ncol(glm_preds)-1)]))))
}

predict_scaled_with_PlastoGram <- function(ngram_models, hmm_models, higher_level_model, test_ngrams, test_sequences, test_df, informative_ngrams) {
  test_ngrams <- add_missing_ngrams(test_ngrams, informative_ngrams)
  ngram_preds_list <- lapply(names(ngram_models), function(ith_model) {
    df <- data.frame(seq_name = test_df[["seq_name"]],
                     pred = predict(ngram_models[[ith_model]], test_ngrams)[["predictions"]][, "TRUE"])
    colnames(df) <- c("seq_name", ith_model)
    df
  }) 
  ngram_models_res <- reduce(ngram_preds_list, full_join, by = "seq_name")
  
  tat_res <- predict_profileHMM(test_sequences, hmm_models[["Tat_model"]], remove_files = FALSE) 
  sec_res <- predict_profileHMM(test_sequences, hmm_models[["Sec_model"]], remove_files = FALSE) 
  
  tat_res <- tat_res %>% 
    mutate(Tat_model = (2^sequence_score) / (1+2^sequence_score),
           seq_name = domain_name) %>% 
    select(c("seq_name", "Tat_model"))
  sec_res <- sec_res %>% 
    mutate(Sec_model = (2^sequence_score) / (1+2^sequence_score),
           seq_name = domain_name) %>% 
    select(c("seq_name", "Sec_model"))
  hmm_models_res <- full_join(sec_res, tat_res, by = "seq_name")
  
  all_res <- left_join(ngram_models_res, hmm_models_res, by = "seq_name")
  all_res[is.na(all_res)] <- 0
  
  scaled_all_res <- scale(select(all_res, -seq_name), center = higher_level_model[["center"]],
                          scale = higher_level_model[["scale"]])
  
  glm_preds <- data.frame(seq_name = test_df[["seq_name"]],
                          dataset = test_df[["dataset"]],
                          predict(higher_level_model[["model"]], scaled_all_res, type = "probs"),
                          Localization = predict(higher_level_model[["model"]], scaled_all_res))
  
  list("Lower-order_models_preds" = all_res,
       "Higher-order_model_preds" = glm_preds,
       "Final_results" = data.frame(glm_preds[, c("seq_name", "dataset", "Localization")],
                                    Probability = sapply(1:nrow(glm_preds[, 3:(ncol(glm_preds)-1)]), function(i) max(glm_preds[i, 3:(ncol(glm_preds)-1)]))))
}
