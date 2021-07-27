evaluate_plastogram <- function(ngram_matrix, target_df, n_fold, N_seqs, P_seqs, 
                                N_TL_TAT_seqs, N_TL_SEC_seqs, with_case_weights = FALSE) {
  
  folds <- lapply(unique(target_df[["dataset"]]), function(ith_dataset) {
    selected <- filter(target_df, dataset == ith_dataset)
    folded <- cvFolds(nrow(selected), K = n_fold)
    data.frame(seq_name = selected[["seq_name"]][folded[["subsets"]]], 
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
                         NP_imp_ngrams,
                         with_case_weights = with_case_weights)
    
    Membrane_imp_ngrams <- calc_imp_ngrams(train_dat, data_df[["membrane_all_target"]][data_df[["fold"]] != ith_fold])
    Membrane_model <- train_rf(train_dat, 
                               data_df[["membrane_all_target"]][data_df[["fold"]] != ith_fold], 
                               Membrane_imp_ngrams,
                               with_case_weights = with_case_weights)
    
    Plastid_membrane_imp_ngrams <- calc_imp_ngrams(ngram_matrix[data_df[["fold"]] != ith_fold & data_df[["NP_target"]] == FALSE & data_df[["membrane_target"]] %in% c("IM", "TM"), ],
                                                   data_df[["membrane_IM_target"]][data_df[["fold"]] != ith_fold & data_df[["NP_target"]] == FALSE & data_df[["membrane_target"]] %in% c("IM", "TM")])
    Plastid_membrane_model <- train_rf(ngram_matrix[data_df[["fold"]] != ith_fold & data_df[["NP_target"]] == FALSE & data_df[["membrane_target"]] %in% c("IM", "TM"), ],
                                       data_df[["membrane_IM_target"]][data_df[["fold"]] != ith_fold & data_df[["NP_target"]] == FALSE & data_df[["membrane_target"]] %in% c("IM", "TM")], 
                                       Plastid_membrane_imp_ngrams,
                                       with_case_weights = with_case_weights)
    
    Nuclear_membrane_imp_ngrams <- get_imp_ngrams_mc(ngram_matrix[data_df[["fold"]] != ith_fold & data_df[["NP_target"]] == TRUE & data_df[["membrane_target"]] != "other", ],
                                                     data_df[data_df[["fold"]] != ith_fold & data_df[["NP_target"]] == TRUE & data_df[["membrane_target"]] != "other", ], 
                                                     "membrane_target")
    Nuclear_membrane_model <- train_rf(ngram_matrix[data_df[["fold"]] != ith_fold & data_df[["NP_target"]] == TRUE & data_df[["membrane_target"]] != "other", ],
                                       data_df[["membrane_target"]][data_df[["fold"]] != ith_fold & data_df[["NP_target"]] == TRUE & data_df[["membrane_target"]] != "other"],
                                       unique(unlist(unname(Nuclear_membrane_imp_ngrams))),
                                       with_case_weights = with_case_weights)
    
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
    left_join(NP_Membrane_results, full_join(hmm_res, all_Membrane_res, by = "seq_name"), by = "seq_name")  %>% 
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


evaluate_plastogram_with_additional_Nuclear_model <- function(ngram_matrix, target_df, n_fold, 
                                                              N_seqs, P_seqs, N_TL_TAT_seqs, N_TL_SEC_seqs,
                                                              with_case_weights = FALSE) {
  
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
                         NP_imp_ngrams,
                         with_case_weights = with_case_weights)
    
    Membrane_imp_ngrams <- calc_imp_ngrams(train_dat, data_df[["membrane_all_target"]][data_df[["fold"]] != ith_fold])
    Membrane_model <- train_rf(train_dat, 
                               data_df[["membrane_all_target"]][data_df[["fold"]] != ith_fold], 
                               Membrane_imp_ngrams,
                               with_case_weights = with_case_weights)
    
    Plastid_membrane_imp_ngrams <- calc_imp_ngrams(ngram_matrix[data_df[["fold"]] != ith_fold & data_df[["NP_target"]] == FALSE & data_df[["membrane_target"]] %in% c("IM", "TM"), ],
                                                   data_df[["membrane_IM_target"]][data_df[["fold"]] != ith_fold & data_df[["NP_target"]] == FALSE & data_df[["membrane_target"]] %in% c("IM", "TM")])
    Plastid_membrane_model <- train_rf(ngram_matrix[data_df[["fold"]] != ith_fold & data_df[["NP_target"]] == FALSE & data_df[["membrane_target"]] %in% c("IM", "TM"), ],
                                       data_df[["membrane_IM_target"]][data_df[["fold"]] != ith_fold & data_df[["NP_target"]] == FALSE & data_df[["membrane_target"]] %in% c("IM", "TM")], 
                                       Plastid_membrane_imp_ngrams,
                                       with_case_weights = with_case_weights)
    
    Nuclear_membrane_imp_ngrams <- get_imp_ngrams_mc(ngram_matrix[data_df[["fold"]] != ith_fold & data_df[["NP_target"]] == TRUE & data_df[["membrane_target"]] != "other", ],
                                                     data_df[data_df[["fold"]] != ith_fold & data_df[["NP_target"]] == TRUE & data_df[["membrane_target"]] != "other", ], 
                                                     "membrane_target")
    Nuclear_membrane_model <- train_rf(ngram_matrix[data_df[["fold"]] != ith_fold & data_df[["NP_target"]] == TRUE & data_df[["membrane_target"]] != "other", ],
                                       data_df[["membrane_target"]][data_df[["fold"]] != ith_fold & data_df[["NP_target"]] == TRUE & data_df[["membrane_target"]] != "other"],
                                       unique(unlist(unname(Nuclear_membrane_imp_ngrams))),
                                       with_case_weights = with_case_weights)
    
    Nuclear_membrane_decision_model_dat <- cbind(data.frame(dataset = as.factor(data_df[["membrane_target"]][data_df[["fold"]] != ith_fold & data_df[["NP_target"]] == TRUE & data_df[["membrane_target"]] != "other"]),
                                                            Nuclear = predict(NP_model, ngram_matrix[data_df[["fold"]] != ith_fold & data_df[["NP_target"]] == TRUE & data_df[["membrane_target"]] != "other", ])[["predictions"]][, "TRUE"],
                                                            Membrane = predict(Membrane_model, ngram_matrix[data_df[["fold"]] != ith_fold & data_df[["NP_target"]] == TRUE & data_df[["membrane_target"]] != "other", ])[["predictions"]][, "TRUE"]),
                                                 predict(Nuclear_membrane_model, ngram_matrix[data_df[["fold"]] != ith_fold & data_df[["NP_target"]] == TRUE & data_df[["membrane_target"]] != "other", ])[["predictions"]])
    Nuclear_membrane_decision_model <- ranger(data = Nuclear_membrane_decision_model_dat, dependent.variable.name = "dataset",
                                              write.forest = TRUE, probability = TRUE, verbose = FALSE, 
                                              class.weights = sapply(levels(Nuclear_membrane_decision_model_dat[["dataset"]]), function(i) 1/sum(Nuclear_membrane_decision_model_dat[["dataset"]] == i)/nrow(Nuclear_membrane_decision_model_dat), USE.NAMES = FALSE),
                                              case.weights = calc_case_weights(Nuclear_membrane_decision_model_dat[["dataset"]]))
    
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
                            predict(Nuclear_membrane_decision_model,
                                    cbind(select(filter(NP_Membrane_results, Nuclear > 0.5 & Membrane > 0.5), c("Nuclear", "Membrane")),
                                          predict(Nuclear_membrane_model, N_Membrane_ngrams)[["predictions"]]))[["predictions"]])
    
    P_Membrane <- filter(NP_Membrane_results, Nuclear <= 0.5 & Membrane > 0.5)[["seq_name"]]
    P_Membrane_ngrams <- ngram_matrix[which(data_df[["seq_name"]] %in% P_Membrane), ]
    P_Membrane_res <- cbind(data.frame(seq_name = data_df[["seq_name"]][which(data_df[["seq_name"]] %in% P_Membrane)]),
                            P_IM = predict(Plastid_membrane_model, P_Membrane_ngrams)[["predictions"]][, "TRUE"])
    all_Membrane_res <- full_join(N_Membrane_res, 
                                  P_Membrane_res,
                                  by = "seq_name")
    left_join(NP_Membrane_results, full_join(hmm_res, all_Membrane_res, by = "seq_name"), by = "seq_name")  %>% 
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


evaluate_plastogram_with_additional_model <- function(ngram_matrix, target_df, n_fold, 
                                                      N_seqs, P_seqs, N_TL_TAT_seqs, N_TL_SEC_seqs,
                                                      with_case_weights = FALSE) {
  
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
                         NP_imp_ngrams,
                         with_case_weights = with_case_weights)
    
    Membrane_imp_ngrams <- calc_imp_ngrams(train_dat, data_df[["membrane_all_target"]][data_df[["fold"]] != ith_fold])
    Membrane_model <- train_rf(train_dat, 
                               data_df[["membrane_all_target"]][data_df[["fold"]] != ith_fold], 
                               Membrane_imp_ngrams,
                               with_case_weights = with_case_weights)
    
    Plastid_membrane_imp_ngrams <- calc_imp_ngrams(ngram_matrix[data_df[["fold"]] != ith_fold & data_df[["NP_target"]] == FALSE & data_df[["membrane_target"]] %in% c("IM", "TM"), ],
                                                   data_df[["membrane_IM_target"]][data_df[["fold"]] != ith_fold & data_df[["NP_target"]] == FALSE & data_df[["membrane_target"]] %in% c("IM", "TM")])
    Plastid_membrane_model <- train_rf(ngram_matrix[data_df[["fold"]] != ith_fold & data_df[["NP_target"]] == FALSE & data_df[["membrane_target"]] %in% c("IM", "TM"), ],
                                       data_df[["membrane_IM_target"]][data_df[["fold"]] != ith_fold & data_df[["NP_target"]] == FALSE & data_df[["membrane_target"]] %in% c("IM", "TM")], 
                                       Plastid_membrane_imp_ngrams,
                                       with_case_weights = with_case_weights)
    
    Nuclear_membrane_imp_ngrams <- get_imp_ngrams_mc(ngram_matrix[data_df[["fold"]] != ith_fold & data_df[["NP_target"]] == TRUE & data_df[["membrane_target"]] != "other", ],
                                                     data_df[data_df[["fold"]] != ith_fold & data_df[["NP_target"]] == TRUE & data_df[["membrane_target"]] != "other", ], 
                                                     "membrane_target")
    Nuclear_membrane_model <- train_rf(ngram_matrix[data_df[["fold"]] != ith_fold & data_df[["NP_target"]] == TRUE & data_df[["membrane_target"]] != "other", ],
                                       data_df[["membrane_target"]][data_df[["fold"]] != ith_fold & data_df[["NP_target"]] == TRUE & data_df[["membrane_target"]] != "other"],
                                       unique(unlist(unname(Nuclear_membrane_imp_ngrams))),
                                       with_case_weights = with_case_weights)
    
    train_tat <- N_TL_TAT_seqs[names(N_TL_TAT_seqs) %in% filter(data_df, fold != ith_fold & Tat_target == TRUE)[["seq_name"]]]
    train_sec <- N_TL_SEC_seqs[names(N_TL_SEC_seqs) %in% filter(data_df, fold != ith_fold & Sec_target == TRUE)[["seq_name"]]]
    Tat_model <- train_profileHMM(train_tat, "tat_model") 
    Sec_model <- train_profileHMM(train_sec, "sec_model") 
    
    # Predict test sequences with 1st tier models
    NP_res <- predict(NP_model, ngram_matrix)[["predictions"]][, "TRUE"]
    Membrane_res <- predict(Membrane_model, ngram_matrix)[["predictions"]][, "TRUE"]
    
    NP_Membrane_results <- data_df %>% 
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
    all_results <- left_join(NP_Membrane_results, full_join(hmm_res, all_Membrane_res, by = "seq_name"), by = "seq_name") 
    all_results[is.na(all_results)] <- 0
    
    decision_model_train_dat <- all_results %>% 
      filter(fold != ith_fold) %>% 
      select(c(Nuclear, Membrane, Tat, Sec, IM, OM, TM, P_IM, dataset)) %>% 
      mutate(dataset = as.factor(dataset))
    decision_model <- ranger(data = decision_model_train_dat, dependent.variable.name = "dataset",
                             write.forest = TRUE, probability = TRUE, verbose = FALSE,
                             class.weights = sapply(levels(decision_model_train_dat[["dataset"]]), function(i) 1/sum(decision_model_train_dat[["dataset"]] == i)/nrow(decision_model_train_dat), USE.NAMES = FALSE))
    decision_preds <- cbind(all_results[all_results[["fold"]] == ith_fold, c("seq_name", "dataset")],
                            predict(decision_model, filter(all_results, fold == ith_fold))[["predictions"]])
    mutate(decision_preds, prediction = c("N_IM", "N_OM", "N_S", "N_TL_SEC", "N_TL_TAT", "N_TM", "P_IM", "P_S", "P_TM")
           [max.col(decision_preds[c("N_IM", "N_OM", "N_S", "N_TL_SEC", "N_TL_TAT", "N_TM", "P_IM", "P_S", "P_TM")])])
    
  }) %>% bind_rows()
  
}


evaluate_plastogram_with_additional_model_for_all <- function(ngram_matrix, target_df, n_fold, 
                                                              N_seqs, P_seqs, N_TL_TAT_seqs, N_TL_SEC_seqs,
                                                              with_case_weights = FALSE) {
  
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
                         NP_imp_ngrams,
                         with_case_weights = with_case_weights)
    
    Membrane_imp_ngrams <- calc_imp_ngrams(train_dat, data_df[["membrane_all_target"]][data_df[["fold"]] != ith_fold])
    Membrane_model <- train_rf(train_dat, 
                               data_df[["membrane_all_target"]][data_df[["fold"]] != ith_fold], 
                               Membrane_imp_ngrams,
                               with_case_weights = with_case_weights)
    
    Plastid_membrane_imp_ngrams <- calc_imp_ngrams(ngram_matrix[data_df[["fold"]] != ith_fold & data_df[["NP_target"]] == FALSE & data_df[["membrane_target"]] %in% c("IM", "TM"), ],
                                                   data_df[["membrane_IM_target"]][data_df[["fold"]] != ith_fold & data_df[["NP_target"]] == FALSE & data_df[["membrane_target"]] %in% c("IM", "TM")])
    Plastid_membrane_model <- train_rf(ngram_matrix[data_df[["fold"]] != ith_fold & data_df[["NP_target"]] == FALSE & data_df[["membrane_target"]] %in% c("IM", "TM"), ],
                                       data_df[["membrane_IM_target"]][data_df[["fold"]] != ith_fold & data_df[["NP_target"]] == FALSE & data_df[["membrane_target"]] %in% c("IM", "TM")], 
                                       Plastid_membrane_imp_ngrams,
                                       with_case_weights = with_case_weights)
    
    Nuclear_membrane_imp_ngrams <- get_imp_ngrams_mc(ngram_matrix[data_df[["fold"]] != ith_fold & data_df[["NP_target"]] == TRUE & data_df[["membrane_target"]] != "other", ],
                                                     data_df[data_df[["fold"]] != ith_fold & data_df[["NP_target"]] == TRUE & data_df[["membrane_target"]] != "other", ], 
                                                     "membrane_target")
    Nuclear_membrane_model <- train_rf(ngram_matrix[data_df[["fold"]] != ith_fold & data_df[["NP_target"]] == TRUE & data_df[["membrane_target"]] != "other", ],
                                       data_df[["membrane_target"]][data_df[["fold"]] != ith_fold & data_df[["NP_target"]] == TRUE & data_df[["membrane_target"]] != "other"],
                                       unique(unlist(unname(Nuclear_membrane_imp_ngrams))),
                                       with_case_weights = with_case_weights)
    
    train_tat <- N_TL_TAT_seqs[names(N_TL_TAT_seqs) %in% filter(data_df, fold != ith_fold & Tat_target == TRUE)[["seq_name"]]]
    train_sec <- N_TL_SEC_seqs[names(N_TL_SEC_seqs) %in% filter(data_df, fold != ith_fold & Sec_target == TRUE)[["seq_name"]]]
    Tat_model <- train_profileHMM(train_tat, "tat_model") 
    Sec_model <- train_profileHMM(train_sec, "sec_model") 
    
    # Predict test sequences with 1st tier models
    NP_res <- predict(NP_model, ngram_matrix)[["predictions"]][, "TRUE"]
    Membrane_res <- predict(Membrane_model, ngram_matrix)[["predictions"]][, "TRUE"]
    
    NP_Membrane_results <- data_df %>% 
      mutate(Nuclear = NP_res,
             Membrane = Membrane_res) 
    
    # Predict test sequences with 2nd tier models
    #N_notmembrane <- c(N_seqs, P_seqs)[filter(NP_Membrane_results, Nuclear > 0.5 & Membrane <= 0.5)[["seq_name"]]]
    Tat_res <- predict_profileHMM(c(N_seqs, P_seqs), "tat_model") %>% 
      mutate(Tat = (2^sequence_score) / (1+2^sequence_score)) %>% 
      select(c("domain_name", "Tat")) %>% 
      setNames(c("seq_name", "Tat"))
    Sec_res <- predict_profileHMM(c(N_seqs, P_seqs), "sec_model") %>% 
      mutate(Sec = (2^sequence_score) / (1+2^sequence_score)) %>% 
      select(c("domain_name", "Sec")) %>% 
      setNames(c("seq_name", "Sec"))
    hmm_res <- left_join(data.frame(seq_name = names(c(N_seqs, P_seqs))), 
                         full_join(Tat_res, Sec_res, by = "seq_name"),
                         by = "seq_name")
    hmm_res[is.na(hmm_res)] <- 0
    
    #N_Membrane <- filter(NP_Membrane_results, Nuclear > 0.3 & Membrane > 0.3)[["seq_name"]]
    #N_Membrane_ngrams <- ngram_matrix[which(data_df[["seq_name"]] %in% N_Membrane), ]
    N_Membrane_res <- cbind(data.frame(seq_name = data_df[["seq_name"]]),
                            predict(Nuclear_membrane_model, ngram_matrix)[["predictions"]])
    
    #P_Membrane <- filter(NP_Membrane_results, Nuclear <= 0.5 & Membrane > 0.5)[["seq_name"]]
    #P_Membrane_ngrams <- ngram_matrix[which(data_df[["seq_name"]] %in% P_Membrane), ]
    P_Membrane_res <- cbind(data.frame(seq_name = data_df[["seq_name"]]),
                            P_IM = predict(Plastid_membrane_model, ngram_matrix)[["predictions"]][, "TRUE"])
    all_Membrane_res <- full_join(N_Membrane_res, 
                                  P_Membrane_res,
                                  by = "seq_name")
    all_results <- left_join(NP_Membrane_results, full_join(hmm_res, all_Membrane_res, by = "seq_name"), by = "seq_name") 
    all_results[is.na(all_results)] <- 0
    
    decision_model_train_dat <- all_results %>% 
      filter(fold != ith_fold) %>% 
      select(c(Nuclear, Membrane, Tat, Sec, IM, OM, TM, P_IM, dataset)) %>% 
      mutate(dataset = as.factor(dataset))
    decision_model <- ranger(data = decision_model_train_dat, dependent.variable.name = "dataset",
                             write.forest = TRUE, probability = TRUE, verbose = FALSE,
                             class.weights = sapply(levels(decision_model_train_dat[["dataset"]]), function(i) 1/sum(decision_model_train_dat[["dataset"]] == i)/nrow(decision_model_train_dat), USE.NAMES = FALSE))
    decision_preds <- cbind(all_results[all_results[["fold"]] == ith_fold, c("seq_name", "dataset")],
                            predict(decision_model, filter(all_results, fold == ith_fold))[["predictions"]])
    mutate(decision_preds, prediction = c("N_IM", "N_OM", "N_S", "N_TL_SEC", "N_TL_TAT", "N_TM", "P_IM", "P_S", "P_TM")
           [max.col(decision_preds[c("N_IM", "N_OM", "N_S", "N_TL_SEC", "N_TL_TAT", "N_TM", "P_IM", "P_S", "P_TM")])])
    
    
  }) %>% bind_rows()
  
}


evaluate_plastogram_with_stroma_model_and_additional_model <- function(ngram_matrix, target_df, n_fold, 
                                                                       N_seqs, P_seqs, N_TL_TAT_seqs, N_TL_SEC_seqs,
                                                                       with_case_weights = FALSE) {
  
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
                         NP_imp_ngrams,
                         with_case_weights = with_case_weights)
    
    Membrane_imp_ngrams <- calc_imp_ngrams(train_dat, data_df[["membrane_all_target"]][data_df[["fold"]] != ith_fold])
    Membrane_model <- train_rf(train_dat, 
                               data_df[["membrane_all_target"]][data_df[["fold"]] != ith_fold], 
                               Membrane_imp_ngrams,
                               with_case_weights = with_case_weights)
    
    Stroma_imp_ngrams <- calc_imp_ngrams(ngram_matrix[data_df[["fold"]] != ith_fold, ],
                                         data_df[["S_target"]][data_df[["fold"]] != ith_fold])
    Stroma_model <- train_rf(ngram_matrix[data_df[["fold"]] != ith_fold, ],
                             data_df[["S_target"]][data_df[["fold"]] != ith_fold], 
                             Stroma_imp_ngrams,
                             with_case_weights = with_case_weights)
    
    Plastid_membrane_imp_ngrams <- calc_imp_ngrams(ngram_matrix[data_df[["fold"]] != ith_fold & data_df[["NP_target"]] == FALSE & data_df[["membrane_target"]] %in% c("IM", "TM"), ],
                                                   data_df[["membrane_IM_target"]][data_df[["fold"]] != ith_fold & data_df[["NP_target"]] == FALSE & data_df[["membrane_target"]] %in% c("IM", "TM")])
    Plastid_membrane_model <- train_rf(ngram_matrix[data_df[["fold"]] != ith_fold & data_df[["NP_target"]] == FALSE & data_df[["membrane_target"]] %in% c("IM", "TM"), ],
                                       data_df[["membrane_IM_target"]][data_df[["fold"]] != ith_fold & data_df[["NP_target"]] == FALSE & data_df[["membrane_target"]] %in% c("IM", "TM")], 
                                       Plastid_membrane_imp_ngrams,
                                       with_case_weights = with_case_weights)
    
    Nuclear_membrane_imp_ngrams <- get_imp_ngrams_mc(ngram_matrix[data_df[["fold"]] != ith_fold & data_df[["NP_target"]] == TRUE & data_df[["membrane_target"]] != "other", ],
                                                     data_df[data_df[["fold"]] != ith_fold & data_df[["NP_target"]] == TRUE & data_df[["membrane_target"]] != "other", ], 
                                                     "membrane_target")
    Nuclear_membrane_model <- train_rf(ngram_matrix[data_df[["fold"]] != ith_fold & data_df[["NP_target"]] == TRUE & data_df[["membrane_target"]] != "other", ],
                                       data_df[["membrane_target"]][data_df[["fold"]] != ith_fold & data_df[["NP_target"]] == TRUE & data_df[["membrane_target"]] != "other"],
                                       unique(unlist(unname(Nuclear_membrane_imp_ngrams))),
                                       with_case_weights = with_case_weights)
    
    train_tat <- N_TL_TAT_seqs[names(N_TL_TAT_seqs) %in% filter(data_df, fold != ith_fold & Tat_target == TRUE)[["seq_name"]]]
    train_sec <- N_TL_SEC_seqs[names(N_TL_SEC_seqs) %in% filter(data_df, fold != ith_fold & Sec_target == TRUE)[["seq_name"]]]
    Tat_model <- train_profileHMM(train_tat, "tat_model") 
    Sec_model <- train_profileHMM(train_sec, "sec_model") 
    
    # Predict test sequences with 1st tier models
    NP_res <- predict(NP_model, ngram_matrix)[["predictions"]][, "TRUE"]
    Membrane_res <- predict(Membrane_model, ngram_matrix)[["predictions"]][, "TRUE"]
    Stroma_res <- predict(Stroma_model, ngram_matrix)[["predictions"]][, "TRUE"]
    
    NP_Membrane_results <- data_df %>% 
      mutate(Nuclear = NP_res,
             Membrane = Membrane_res,
             Stroma = Stroma_res) 
    
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
    all_results <- left_join(NP_Membrane_results, full_join(hmm_res, all_Membrane_res, by = "seq_name"), by = "seq_name") 
    all_results[is.na(all_results)] <- 0
    
    decision_model_train_dat <- all_results %>% 
      filter(fold != ith_fold) %>% 
      select(c(Nuclear, Membrane, Stroma, Tat, Sec, IM, OM, TM, P_IM, dataset)) %>% 
      mutate(dataset = as.factor(dataset))
    decision_model <- ranger(data = decision_model_train_dat, dependent.variable.name = "dataset",
                             write.forest = TRUE, probability = TRUE, verbose = FALSE,
                             class.weights = sapply(levels(decision_model_train_dat[["dataset"]]), function(i) 1/sum(decision_model_train_dat[["dataset"]] == i)/nrow(decision_model_train_dat), USE.NAMES = FALSE))
    decision_preds <- cbind(all_results[all_results[["fold"]] == ith_fold, c("seq_name", "dataset")],
                            predict(decision_model, filter(all_results, fold == ith_fold))[["predictions"]])
    mutate(decision_preds, prediction = c("N_IM", "N_OM", "N_S", "N_TL_SEC", "N_TL_TAT", "N_TM", "P_IM", "P_S", "P_TM")
           [max.col(decision_preds[c("N_IM", "N_OM", "N_S", "N_TL_SEC", "N_TL_TAT", "N_TM", "P_IM", "P_S", "P_TM")])])
    
  }) %>% bind_rows()
  
}

evaluate_plastogram_with_OM_stroma_model_and_additional_model <- function(ngram_matrix, target_df, n_fold, 
                                                                       N_seqs, P_seqs, N_TL_TAT_seqs, N_TL_SEC_seqs,
                                                                       with_case_weights = FALSE) {
  
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
                         NP_imp_ngrams,
                         with_case_weights = with_case_weights)
    
    Membrane_imp_ngrams <- calc_imp_ngrams(train_dat, data_df[["membrane_all_target"]][data_df[["fold"]] != ith_fold])
    Membrane_model <- train_rf(train_dat, 
                               data_df[["membrane_all_target"]][data_df[["fold"]] != ith_fold], 
                               Membrane_imp_ngrams,
                               with_case_weights = with_case_weights)
    
    OM_Stroma_imp_ngrams <- calc_imp_ngrams(ngram_matrix[data_df[["fold"]] != ith_fold & (data_df[["membrane_target"]] == "OM" | data_df[["S_target"]] == "TRUE"), ],
                                         data_df[["membrane_OM_target"]][data_df[["fold"]] != ith_fold & (data_df[["membrane_target"]] == "OM" | data_df[["S_target"]] == "TRUE")])
    OM_Stroma_model <- train_rf(ngram_matrix[data_df[["fold"]] != ith_fold & (data_df[["membrane_target"]] == "OM" | data_df[["S_target"]] == "TRUE"), ],
                             data_df[["membrane_OM_target"]][data_df[["fold"]] != ith_fold & (data_df[["membrane_target"]] == "OM" | data_df[["S_target"]] == "TRUE")], 
                             OM_Stroma_imp_ngrams,
                             with_case_weights = with_case_weights)
    
    Plastid_membrane_imp_ngrams <- calc_imp_ngrams(ngram_matrix[data_df[["fold"]] != ith_fold & data_df[["NP_target"]] == FALSE & data_df[["membrane_target"]] %in% c("IM", "TM"), ],
                                                   data_df[["membrane_IM_target"]][data_df[["fold"]] != ith_fold & data_df[["NP_target"]] == FALSE & data_df[["membrane_target"]] %in% c("IM", "TM")])
    Plastid_membrane_model <- train_rf(ngram_matrix[data_df[["fold"]] != ith_fold & data_df[["NP_target"]] == FALSE & data_df[["membrane_target"]] %in% c("IM", "TM"), ],
                                       data_df[["membrane_IM_target"]][data_df[["fold"]] != ith_fold & data_df[["NP_target"]] == FALSE & data_df[["membrane_target"]] %in% c("IM", "TM")], 
                                       Plastid_membrane_imp_ngrams,
                                       with_case_weights = with_case_weights)
    
    Nuclear_membrane_imp_ngrams <- get_imp_ngrams_mc(ngram_matrix[data_df[["fold"]] != ith_fold & data_df[["NP_target"]] == TRUE & data_df[["membrane_target"]] != "other", ],
                                                     data_df[data_df[["fold"]] != ith_fold & data_df[["NP_target"]] == TRUE & data_df[["membrane_target"]] != "other", ], 
                                                     "membrane_target")
    Nuclear_membrane_model <- train_rf(ngram_matrix[data_df[["fold"]] != ith_fold & data_df[["NP_target"]] == TRUE & data_df[["membrane_target"]] != "other", ],
                                       data_df[["membrane_target"]][data_df[["fold"]] != ith_fold & data_df[["NP_target"]] == TRUE & data_df[["membrane_target"]] != "other"],
                                       unique(unlist(unname(Nuclear_membrane_imp_ngrams))),
                                       with_case_weights = with_case_weights)
    
    train_tat <- N_TL_TAT_seqs[names(N_TL_TAT_seqs) %in% filter(data_df, fold != ith_fold & Tat_target == TRUE)[["seq_name"]]]
    train_sec <- N_TL_SEC_seqs[names(N_TL_SEC_seqs) %in% filter(data_df, fold != ith_fold & Sec_target == TRUE)[["seq_name"]]]
    Tat_model <- train_profileHMM(train_tat, "tat_model") 
    Sec_model <- train_profileHMM(train_sec, "sec_model") 
    
    # Predict test sequences with 1st tier models
    NP_res <- predict(NP_model, ngram_matrix)[["predictions"]][, "TRUE"]
    Membrane_res <- predict(Membrane_model, ngram_matrix)[["predictions"]][, "TRUE"]
    
    NP_Membrane_results <- data_df %>% 
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
    results <- left_join(NP_Membrane_results, full_join(hmm_res, all_Membrane_res, by = "seq_name"), by = "seq_name") 
    Stroma_NMembrane <- c(N_Membrane, filter(hmm_res, Sec == 0 & Tat == 0)[["seq_name"]])
    Stroma_NMembrane_ngrams <- ngram_matrix[which(data_df[["seq_name"]] %in% Stroma_NMembrane), ]
    Stroma_NMembrane_res <- data.frame(seq_name = Stroma_NMembrane,
                                       OM_Stroma = predict(OM_Stroma_model, Stroma_NMembrane_ngrams)[["predictions"]][, "TRUE"])
    all_results <- left_join(results, Stroma_NMembrane_res)
    all_results[is.na(all_results)] <- 0
    
    decision_model_train_dat <- all_results %>% 
      filter(fold != ith_fold) %>% 
      select(c(Nuclear, Membrane, OM_Stroma, Tat, Sec, IM, OM, TM, P_IM, dataset)) %>% 
      mutate(dataset = as.factor(dataset))
    decision_model <- ranger(data = decision_model_train_dat, dependent.variable.name = "dataset",
                             write.forest = TRUE, probability = TRUE, verbose = FALSE,
                             class.weights = sapply(levels(decision_model_train_dat[["dataset"]]), function(i) 1/sum(decision_model_train_dat[["dataset"]] == i)/nrow(decision_model_train_dat), USE.NAMES = FALSE))
    decision_preds <- cbind(all_results[all_results[["fold"]] == ith_fold, c("seq_name", "dataset")],
                            predict(decision_model, filter(all_results, fold == ith_fold))[["predictions"]])
    mutate(decision_preds, prediction = c("N_IM", "N_OM", "N_S", "N_TL_SEC", "N_TL_TAT", "N_TM", "P_IM", "P_S", "P_TM")
           [max.col(decision_preds[c("N_IM", "N_OM", "N_S", "N_TL_SEC", "N_TL_TAT", "N_TM", "P_IM", "P_S", "P_TM")])])
    
  }) %>% bind_rows()
  
}


evaluate_plastogram_with_stroma_model_and_additional_model_for_all <- function(ngram_matrix, target_df, n_fold, 
                                                                               N_seqs, P_seqs, N_TL_TAT_seqs, N_TL_SEC_seqs,
                                                                               with_case_weights = FALSE) {
  
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
                         NP_imp_ngrams,
                         with_case_weights = with_case_weights)
    
    Membrane_imp_ngrams <- calc_imp_ngrams(train_dat, data_df[["membrane_all_target"]][data_df[["fold"]] != ith_fold])
    Membrane_model <- train_rf(train_dat, 
                               data_df[["membrane_all_target"]][data_df[["fold"]] != ith_fold], 
                               Membrane_imp_ngrams,
                               with_case_weights = with_case_weights)
    
    Plastid_membrane_imp_ngrams <- calc_imp_ngrams(ngram_matrix[data_df[["fold"]] != ith_fold & data_df[["NP_target"]] == FALSE & data_df[["membrane_target"]] %in% c("IM", "TM"), ],
                                                   data_df[["membrane_IM_target"]][data_df[["fold"]] != ith_fold & data_df[["NP_target"]] == FALSE & data_df[["membrane_target"]] %in% c("IM", "TM")])
    Plastid_membrane_model <- train_rf(ngram_matrix[data_df[["fold"]] != ith_fold & data_df[["NP_target"]] == FALSE & data_df[["membrane_target"]] %in% c("IM", "TM"), ],
                                       data_df[["membrane_IM_target"]][data_df[["fold"]] != ith_fold & data_df[["NP_target"]] == FALSE & data_df[["membrane_target"]] %in% c("IM", "TM")], 
                                       Plastid_membrane_imp_ngrams,
                                       with_case_weights = with_case_weights)
    
    Nuclear_membrane_imp_ngrams <- get_imp_ngrams_mc(ngram_matrix[data_df[["fold"]] != ith_fold & data_df[["NP_target"]] == TRUE & data_df[["membrane_target"]] != "other", ],
                                                     data_df[data_df[["fold"]] != ith_fold & data_df[["NP_target"]] == TRUE & data_df[["membrane_target"]] != "other", ], 
                                                     "membrane_target")
    Nuclear_membrane_model <- train_rf(ngram_matrix[data_df[["fold"]] != ith_fold & data_df[["NP_target"]] == TRUE & data_df[["membrane_target"]] != "other", ],
                                       data_df[["membrane_target"]][data_df[["fold"]] != ith_fold & data_df[["NP_target"]] == TRUE & data_df[["membrane_target"]] != "other"],
                                       unique(unlist(unname(Nuclear_membrane_imp_ngrams))),
                                       with_case_weights = with_case_weights)
    
    Stroma_imp_ngrams <- calc_imp_ngrams(ngram_matrix[data_df[["fold"]] != ith_fold, ],
                                         data_df[["S_target"]][data_df[["fold"]] != ith_fold])
    Stroma_model <- train_rf(ngram_matrix[data_df[["fold"]] != ith_fold, ],
                             data_df[["S_target"]][data_df[["fold"]] != ith_fold], 
                             Stroma_imp_ngrams,
                             with_case_weights = with_case_weights)
    
    
    train_tat <- N_TL_TAT_seqs[names(N_TL_TAT_seqs) %in% filter(data_df, fold != ith_fold & Tat_target == TRUE)[["seq_name"]]]
    train_sec <- N_TL_SEC_seqs[names(N_TL_SEC_seqs) %in% filter(data_df, fold != ith_fold & Sec_target == TRUE)[["seq_name"]]]
    Tat_model <- train_profileHMM(train_tat, "tat_model") 
    Sec_model <- train_profileHMM(train_sec, "sec_model") 
    
    # Predict test sequences with 1st tier models
    NP_res <- predict(NP_model, ngram_matrix)[["predictions"]][, "TRUE"]
    Membrane_res <- predict(Membrane_model, ngram_matrix)[["predictions"]][, "TRUE"]
    Stroma_res <- predict(Stroma_model, ngram_matrix)[["predictions"]][, "TRUE"]
    
    NP_Membrane_results <- data_df %>% 
      mutate(Nuclear = NP_res,
             Membrane = Membrane_res,
             Stroma = Stroma_res) 
    
    # Predict test sequences with 2nd tier models
    #N_notmembrane <- c(N_seqs, P_seqs)[filter(NP_Membrane_results, Nuclear > 0.5 & Membrane <= 0.5)[["seq_name"]]]
    Tat_res <- predict_profileHMM(c(N_seqs, P_seqs), "tat_model") %>% 
      mutate(Tat = (2^sequence_score) / (1+2^sequence_score)) %>% 
      select(c("domain_name", "Tat")) %>% 
      setNames(c("seq_name", "Tat"))
    Sec_res <- predict_profileHMM(c(N_seqs, P_seqs), "sec_model") %>% 
      mutate(Sec = (2^sequence_score) / (1+2^sequence_score)) %>% 
      select(c("domain_name", "Sec")) %>% 
      setNames(c("seq_name", "Sec"))
    hmm_res <- left_join(data.frame(seq_name = names(c(N_seqs, P_seqs))), 
                         full_join(Tat_res, Sec_res, by = "seq_name"),
                         by = "seq_name")
    hmm_res[is.na(hmm_res)] <- 0
    
    #N_Membrane <- filter(NP_Membrane_results, Nuclear > 0.4 & Membrane > 0.3)[["seq_name"]]
    #N_Membrane_ngrams <- ngram_matrix[which(data_df[["seq_name"]] %in% N_Membrane), ]
    N_Membrane_res <- cbind(data.frame(seq_name = data_df[["seq_name"]]),
                            predict(Nuclear_membrane_model, ngram_matrix)[["predictions"]])
    
    #P_Membrane <- filter(NP_Membrane_results, Nuclear <= 0.5 & Membrane > 0.5)[["seq_name"]]
    #P_Membrane_ngrams <- ngram_matrix[which(data_df[["seq_name"]] %in% P_Membrane), ]
    P_Membrane_res <- cbind(data.frame(seq_name = data_df[["seq_name"]]),
                            P_IM = predict(Plastid_membrane_model, ngram_matrix)[["predictions"]][, "TRUE"])
    all_Membrane_res <- full_join(N_Membrane_res, 
                                  P_Membrane_res,
                                  by = "seq_name")
    all_results <- left_join(NP_Membrane_results, full_join(hmm_res, all_Membrane_res, by = "seq_name"), by = "seq_name") 
    all_results[is.na(all_results)] <- 0
    
    decision_model_train_dat <- all_results %>% 
      filter(fold != ith_fold) %>% 
      select(c(Nuclear, Membrane, Stroma, Tat, Sec, IM, OM, TM, P_IM, dataset)) %>% 
      mutate(dataset = as.factor(dataset))
    decision_model <- ranger(data = decision_model_train_dat, dependent.variable.name = "dataset",
                             write.forest = TRUE, probability = TRUE, verbose = FALSE,
                             class.weights = sapply(levels(decision_model_train_dat[["dataset"]]), function(i) 1/sum(decision_model_train_dat[["dataset"]] == i)/nrow(decision_model_train_dat), USE.NAMES = FALSE))
    decision_preds <- cbind(all_results[all_results[["fold"]] == ith_fold, c("seq_name", "dataset")],
                            predict(decision_model, filter(all_results, fold == ith_fold))[["predictions"]])
    mutate(decision_preds, prediction = c("N_IM", "N_OM", "N_S", "N_TL_SEC", "N_TL_TAT", "N_TM", "P_IM", "P_S", "P_TM")
           [max.col(decision_preds[c("N_IM", "N_OM", "N_S", "N_TL_SEC", "N_TL_TAT", "N_TM", "P_IM", "P_S", "P_TM")])])
    
    
  }) %>% bind_rows()
  
}
