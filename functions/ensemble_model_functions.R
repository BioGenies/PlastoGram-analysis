create_folds <- function(target_df, n_fold) {
  folds <- lapply(unique(target_df[["dataset"]]), function(ith_dataset) {
    selected <- filter(target_df, dataset == ith_dataset)
    folded <- cvFolds(nrow(selected), K = n_fold)
    data.frame(seq_name = selected[["seq_name"]][folded[["subsets"]]], 
               fold = folded[["which"]],
               stringsAsFactors = FALSE)
  }) %>% bind_rows()
  left_join(target_df, folds, by = "seq_name")
}

train_ngram_models <- function(model_df, ngram_matrix, data_df, ith_fold) {
  df <- filter(model_df, `input.data` == "ngram_matrix")
  lapply(1:nrow(df), function(i) {
    if(df[["input.data.filtering"]][i] != "") {
      filtering <- paste0("fold != ", ith_fold, " & ", df[["input.data.filtering"]][i])
    } else {
      filtering <- paste0("fold != ", ith_fold)
    }
    selected <- filter(data_df, eval(parse(text = filtering)))[["seq_name"]]
    target <- data_df[[df[["target.name"]][i]]][which(data_df[["seq_name"]] %in% selected)]
    train_dat <- ngram_matrix[which(data_df[["seq_name"]] %in% selected), ]
    if(df[["multiclass"]][i] == TRUE) {
      train_rf(train_dat,
               target,
               unique(
                 unlist(
                   unname(
                     get_imp_ngrams_mc(train_dat, 
                                       data_df[which(data_df[["seq_name"]] %in% selected), ],
                                       df[["target.name"]][i])))),
               with_class_weights = df[["class_weights"]][i])
    } else {
      train_rf(train_dat,
               target,
               calc_imp_ngrams(train_dat, as.logical(target)),
               with_class_weights = df[["class_weights"]][i])
    }
  }) %>% setNames(df[["model.name"]])
}

train_profile_HMM_models <- function(model_df, sequences, data_df, ith_fold) {
  df <- filter(model_df, `input.data` == "sequences")
  lapply(1:nrow(df), function(i) {
    if(df[["input.data.filtering"]][i] != "") {
      filtering <- paste0("fold != ", ith_fold, " & ", df[["input.data.filtering"]][i])
    } else {
      filtering <- paste0("fold != ", ith_fold)
    }
    train_seq <- sequences[names(sequences) %in% filter(data_df, eval(parse(text = filtering)))[["seq_name"]]]
    train_profileHMM(train_seq, df[["model.name"]][i])
    df[["model.name"]][i]
  }) %>% setNames(df[["model.name"]])
}

is_empty <- function(x) {
  x == "" | is.na(x)
}

calculate_evaluation_order <- function(model_df) {
  model_names <- model_df[["model.name"]]
  target_names <- filter(df, multiclass %in% c(FALSE, NA))[["target.name"]]
  res_names <- sapply(target_names, function(i) gsub("_target", "", i), USE.NAMES = FALSE)
  
  orders <- rep(NA, length(model_names))
  for (i in seq_along(model_names)) {
    if(!(any(sapply(res_names, function(j) grepl(j, model_df[["test.filtering"]][i]))))) {
      orders[i] <- 1
    } 
  }
  calculated_res <- res_names[which(!(is.na(orders)))]
  current_order <- 2
  while(!(all(res_names %in% calculated_res))) {
    for (i in which(is.na(orders))) {
      needed_res <- res_names[sapply(res_names, function(j) grepl(j, model_df[["test.filtering"]][i]))]
      if(all(needed_res %in% calculated_res)) {
        orders[i] <- current_order
      } 
    }
    current_order <- current_order + 1
    calculated_res <- res_names[which(!(is.na(orders)))]
  }
  orders
}

predict_with_models <- function(model_df, ngram_matrix, sequences, data_df, ith_fold, ngram_models) {
  model_df <- mutate(model_df, `evaluation.order` = calculate_evaluation_order(model_df))
  max_order <- max(model_df[["evaluation.order"]])
  current_order <- 1
  data_df_results <- data_df
  while(current_order < max_order) {
    df <- filter(model_df, `evaluation.order` == current_order)
    results <- lapply(1:nrow(df), function(i) {
      if(df[["input.data"]][i] == "ngram_matrix") {
        selected <- if(is_empty(df[["test.filtering"]][i])) {
          data_df_results[["seq_name"]]
        } else {
          filter(data_df_results, eval(parse(text = gsub("ith_fold", "get('ith_fold')", df[["test.filtering"]][i]))))[["seq_name"]]
        }
        test_dat <- ngram_matrix[which(data_df_results[["seq_name"]] %in% selected), ]
        if(df[["multiclass"]][i] == TRUE) {
          cbind(data.frame(seq_name = selected),
                         predict(ngram_models[[df[["model.name"]][i]]], test_dat)[["predictions"]])
        } else {
          data.frame(seq_name = selected,
                              pred = predict(ngram_models[[df[["model.name"]][i]]], test_dat)[["predictions"]][, "TRUE"]) %>% 
            setNames(c("seq_name", gsub("_target", "", df[["target.name"]][i])))
        }
        
      } else if(df[["input.data"]][i] == "sequences") {
        test_seqs <- if(is_empty(df[["test.filtering"]][i])) {
          sequences
        } else {
          sequences[names(sequences) %in% filter(data_df_results, eval(parse(text = gsub("ith_fold", "get('ith_fold')", df[["test.filtering"]][i]))))[["seq_name"]]]
        }
        predict_profileHMM(test_seqs, df[["model.name"]][i]) %>% 
          mutate(pred = (2^sequence_score) / (1+2^sequence_score)) %>% 
          select(c("domain_name", "pred")) %>% 
          setNames(c("seq_name", gsub("_target", "", df[["target.name"]][i])))
      }
      
    }) %>% reduce(., full_join, by = "seq_name")
    data_df_results <- reduce(list(data_df_results, results), left_join, by = "seq_name")
    current_order <- current_order + 1
  } 
  data_df_results
}

get_decision <- function(results) {
  mutate(results, Decision = case_when(
    NP > 0.5 & membrane_all <= 0.5 & Tat > Sec ~ "N_TL_TAT",
    NP > 0.5 & membrane_all <= 0.5 & Sec > Tat ~ "N_TL_SEC",
    NP > 0.5 & membrane_all <= 0.5 & Sec == 0 & Tat == 0 ~ "N_S",
    NP <= 0.5 & membrane_all <= 0.5 ~ "P_S",
    NP <= 0.5 & membrane_all > 0.5 & membrane_IM > 0.5 ~ "P_IM",
    NP <= 0.5 & membrane_all > 0.5 & membrane_IM <= 0.5 ~ "P_TM",
    NP > 0.5 & membrane_all > 0.5 & OM > IM & OM > TM ~ "N_OM",
    NP > 0.5 & membrane_all > 0.5 & IM > OM & IM > TM ~ "N_IM",
    NP > 0.5 & membrane_all > 0.5 & TM > IM & TM > IM ~ "N_TM"
  ))
}


evaluate_plastogram <- function(ngram_matrix, data_df, sequences, desc_file) {
  
  model_df <- read.csv(desc_file)
  
  plastogram_predictions <- lapply(unique(data_df[["fold"]]), function(ith_fold) {
    
    # Train models
    ngram_models <- train_ngram_models(model_df, ngram_matrix, data_df, ith_fold)
    profile_HMM_models <- train_profile_HMM_models(model_df, sequences, data_df, ith_fold)
    
    # Predict test sequences
    results <- predict_with_models(model_df, ngram_matrix, c(N_seqs, P_seqs), data_df, ith_fold, ngram_models)
    results[is.na(results)] <- 0
    # If no decision model is used, get final prediction based on our decision rules
    if(!("predictions" %in% model_df[["input.data"]])) {
      final_res <- get_decision(filter(results, fold == ith_fold))
    } else { 
      # Train additional model if a decision model is used
      train_dat <- filter(results, fold != ith_fold) %>% 
        select(-c("seq_name", "NP_target", "Sec_target", "Tat_target", "membrane_target", 
                  "S_target", "membrane_all_target", "membrane_OM_target", "membrane_IM_target", 
                  "membrane_TM_target", "fold")) %>% 
        mutate(dataset = as.factor(dataset))
      if(model_df[["class_weights"]][which(model_df[["input.data"]] == "predictions")] == TRUE) {
        cw <- sapply(levels(train_dat[["dataset"]]), function(i) 1/sum(train_dat[["dataset"]] == i)/nrow(train_dat), USE.NAMES = FALSE)
      } else {
        cw <- NULL
      }
      decision_model <- ranger(data = train_dat, dependent.variable.name = "dataset",
                               write.forest = TRUE, probability = TRUE, verbose = FALSE,
                               class.weights = cw)
      test_dat <- filter(results, fold == ith_fold) %>% 
        select(c("seq_name", "dataset", "NP_target", "Sec_target", "Tat_target", "membrane_target", 
                 "S_target", "membrane_all_target", "membrane_OM_target", "membrane_IM_target", 
                 "membrane_TM_target", "fold"))
      final_res <- cbind(test_dat, predict(decision_model, filter(results, fold == ith_fold))[["predictions"]])
      mutate(final_res, Decision = c("N_IM", "N_OM", "N_S", "N_TL_SEC", "N_TL_TAT", "N_TM", "P_IM", "P_S", "P_TM")
             [max.col(final_res[c("N_IM", "N_OM", "N_S", "N_TL_SEC", "N_TL_TAT", "N_TM", "P_IM", "P_S", "P_TM")])])
    }
  }) %>% bind_rows()
  
  write.csv(plastogram_predictions, gsub(".csv", "_predictions.csv", desc_file), row.names = FALSE)
  
  plastogram_predictions
  
}


test_ensembles <- function(desc_files_dir, ngram_matrix, data_df, sequences) {
  desc_files <- list.files(desc_files_dir, full.names = TRUE)
  lapply(desc_files, function(ith_file) {
    res <- evaluate_plastogram(ngram_matrix, data_df, sequences, ith_file)
    data.frame(ensemble = gsub(".csv", "", last(strsplit(ith_file, "/")[[1]])),
               kappa = KAPPA(res[["dataset"]], res[["Decision"]]))
  }) %>% bind_rows()
}
