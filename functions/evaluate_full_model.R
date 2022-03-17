predict_with_PlastoGram <- function(ngram_models, hmm_models, higher_level_model, test_ngrams, test_sequences, test_df) {
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
    mutate(Tat = (2^sequence_score) / (1+2^sequence_score),
           seq_name = domain_name) %>% 
    select(c("seq_name", "Tat"))
  sec_res <- sec_res %>% 
    mutate(Sec = (2^sequence_score) / (1+2^sequence_score),
           seq_name = domain_name) %>% 
    select(c("seq_name", "Sec"))
  hmm_models_res <- full_join(sec_res, tat_res, by = "seq_name")
  
  all_res <- left_join(ngram_models_res, hmm_models_res, by = "seq_name")
  all_res[is.na(all_res)] <- 0
  
  glm_preds <- data.frame(seq_name = test_df[["seq_name"]],
                          predict(higher_level_model, all_res, type = "probs"),
                          Localization = predict(higher_level_model, all_res))
  
  list("Lower-order_models_preds" = all_res,
       "Higher-order_model_preds" = glm_preds,
       "Final_results" = data.frame(glm_preds[, c("seq_name", "Localization")],
                                    Probability = sapply(1:nrow(glm_preds[, 2:(ncol(glm_preds)-1)]), function(i) max(glm_preds[i, 2:(ncol(glm_preds)-1)]))))
}

