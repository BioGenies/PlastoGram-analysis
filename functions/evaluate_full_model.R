evaluate_full_model <- function(ngram_matrix, target_df, n_fold) {
  
  folds <- lapply(unique(target_df[["dataset"]]), function(ith_dataset) {
    selected <- filter(target_df, dataset == ith_dataset)
    folded <- cvFolds(nrow(selected), K = n_fold)
    fold_df <- data.frame(seq_name = selected[["seq_name"]][folded[["subsets"]]], 
                          fold = folded[["which"]],
                          stringsAsFactors = FALSE)
  }) %>% bind_rows()
  
  data_df <- left_join(target_df, folds, by = "seq_name")
  
  lapply(unique(data_df[["fold"]]), function(ith_fold) {
    
    # Train all models 
    train_dat <- ngram_matrix[data_df[["fold"]] != ith_fold, ]
    test_dat <- ngram_matrix[data_df[["fold"]] == ith_fold, ]
    
    NP_imp_ngrams <- calc_imp_ngrams(train_dat, data_df[["NP_target"]][data_df[["fold"]] != ith_fold])
    NP_model <- train_rf(train_dat, 
                         data_df[["NP_target"]][data_df[["fold"]] != ith_fold], 
                         NP_imp_ngrams)
    
    Membrane_imp_ngrams <- calc_imp_ngrams(train_dat, data_df[["membrane_all_target"]][data_df[["fold"]] != ith_fold])
    Membrane_model <- train_rf(train_dat, 
                               data_df[["membrane_all_target"]][data_df[["fold"]] != ith_fold], 
                               Membrane_imp_ngrams)
    
    Plastid_membrane_imp_ngrams <- calc_imp_ngrams(ngram_matrix[data_df[["fold"]] != ith_fold & data_df[["NP_target"]] == FALSE & data_df[["membrane_target"]] %in% c("IM", "TM"), ],
                                                   data_df[["membrane_IM_target"]][data_df[["fold"]] != ith_fold & data_df[["NP_target"]] == FALSE & data_df[["membrane_target"]] %in% c("IM", "TM")])
    Plastid_membrane_model <- train_rf(ngram_matrix[data_df[["fold"]] != ith_fold & data_df[["NP_target"]] == FALSE & data_df[["membrane_target"]] %in% c("IM", "TM"), ],
                                       data_df[["membrane_IM_target"]][data_df[["fold"]] != ith_fold & data_df[["NP_target"]] == FALSE & data_df[["membrane_target"]] %in% c("IM", "TM")], 
                                       Plastid_membrane_imp_ngrams)
    
    Nuclear_membrane_imp_ngrams <- get_imp_ngrams_mc(ngram_matrix[data_df[["fold"]] != ith_fold & data_df[["NP_target"]] == TRUE & data_df[["membrane_target"]] != "other", ],
                                                     data_df[data_df[["fold"]] != ith_fold & data_df[["NP_target"]] == TRUE & data_df[["membrane_target"]] != "other", ], 
                                                     "membrane_target")
    Nuclear_membrane_model <- train_rf(ngram_matrix[data_df[["fold"]] != ith_fold & data_df[["NP_target"]] == TRUE & data_df[["membrane_target"]] != "other", ],
                                       data_df[["membrane_target"]][data_df[["fold"]] != ith_fold & data_df[["NP_target"]] == TRUE & data_df[["membrane_target"]] != "other"],
                                       unique(unlist(unname(Nuclear_membrane_imp_ngrams))))
    
    train_tat <- N_TL_TAT_seqs[names(N_TL_TAT_seqs) %in% filter(data_df, fold != ith_fold & Tat_target == TRUE)[["seq_name"]]]
    train_sec <- N_TL_SEC_seqs[names(N_TL_SEC_seqs) %in% filter(data_df, fold != ith_fold & Sec_target == TRUE)[["seq_name"]]]
    Tat_model <- train_profileHMM(train_tat, "tat_model") 
    Sec_model <- train_profileHMM(train_sec, "sec_model") 
    
    # Predict test sequences with 1st tier models
    NP_res <- predict(NP_model, test_dat)[["predictions"]][, "TRUE"]
    Membrane_res <- predict(Membrane_model, test_dat)[["predictions"]][, "TRUE"]
    
    NP_Membrane_results <- filter(data_df, fold == ith_fold) %>% 
      mutate(Nuclear = NP_res,
             Membrane = Membrane_res) 
    
    # Predict test sequences with 2nd tier models
    N_notmembrane <- c(N_seqs, P_seqs)[filter(NP_Membrane_results, Nuclear > 0.5 & Membrane <= 0.5)[["seq_name"]]]
    Tat_res <- predict_profileHMM(N_notmembrane, "tat_model") %>% 
      mutate(Tat = (2^sequence_score) / (1+2^sequence_score)) %>% 
      select(c("domain_name", "Tat")) %>% 
      setNames(c("seq_name", "Tat"))
    Sec_res <- predict_profileHMM(N_notmembrane, "sec_model") %>% 
      mutate(Sec = (2^sequence_score) / (1+2^sequence_score)) %>% 
      select(c("domain_name", "Sec")) %>% 
      setNames(c("seq_name", "Sec"))
    hmm_res <- left_join(data.frame(seq_name = names(N_notmembrane)), 
                         full_join(Tat_res, Sec_res, by = "seq_name"),
                         by = "seq_name")
    hmm_res[is.na(hmm_res)] <- 0
    
    N_Membrane <- filter(NP_Membrane_results, Nuclear > 0.5 & Membrane > 0.5)[["seq_name"]]
    N_Membrane_ngrams <- ngram_matrix[which(data_df[["seq_name"]] %in% N_Membrane), ]
    N_Membrane_res <- cbind(data.frame(seq_name = data_df[["seq_name"]][which(data_df[["seq_name"]] %in% N_Membrane)]),
                            predict(Nuclear_membrane_model, N_Membrane_ngrams)[["predictions"]])
    
    P_Membrane <- filter(NP_Membrane_results, Nuclear <= 0.5 & Membrane > 0.5)[["seq_name"]]
    P_Membrane_ngrams <- ngram_matrix[which(data_df[["seq_name"]] %in% P_Membrane), ]
    P_Membrane_res <- cbind(data.frame(seq_name = data_df[["seq_name"]][which(data_df[["seq_name"]] %in% P_Membrane)]),
                            P_IM = predict(Plastid_membrane_model, P_Membrane_ngrams)[["predictions"]][, "TRUE"])
    all_Membrane_res <- full_join(N_Membrane_res, 
                                  P_Membrane_res,
                                  by = "seq_name")
    left_join(NP_Membrane_results, full_join(hmm_res, all_Membrane_res, by = "seq_name"), by = "seq_name") %>% 
      mutate(Decision = case_when(
        Nuclear > 0.5 & Membrane <= 0.5 & Tat > Sec ~ "N_TL_TAT",
        Nuclear > 0.5 & Membrane <= 0.5 & Sec > Tat ~ "N_TL_SEC",
        Nuclear > 0.5 & Membrane <= 0.5 & Sec == 0 & Tat == 0 ~ "N_S",
        Nuclear <= 0.5 & Membrane <= 0.5 ~ "P_S",
        Nuclear <= 0.5 & Membrane > 0.5 & P_IM > 0.5 ~ "P_IM",
        Nuclear <= 0.5 & Membrane > 0.5 & P_IM <= 0.5 ~ "P_TM",
        Nuclear > 0.5 & Membrane > 0.5 & OM > IM & OM > TM ~ "N_OM",
        Nuclear > 0.5 & Membrane > 0.5 & IM > OM & IM > TM ~ "N_IM",
        Nuclear > 0.5 & Membrane > 0.5 & TM > IM & TM > IM ~ "N_TM"
      ))
    
  }) %>% bind_rows()
  
}
