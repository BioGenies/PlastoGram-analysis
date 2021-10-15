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

train_ngram_models <- function(model_df, ngram_matrix, data_df, filtering_colname, filtering_term) {
  df <- filter(model_df, Input_data == "ngrams")
  lapply(1:nrow(df), function(i) {
    selected <- if(is.null(filtering_colname) & is.null(filtering_term)) {
      if(df[["Input_filtering"]][i] != "") {
        filtering <- paste0(df[["Input_filtering"]][i])
        filter(data_df, eval(parse(text = filtering)))[["seq_name"]]
      } else {
        data_df[["seq_name"]]
      }
    } else { 
      filtering <- if(df[["Input_filtering"]][i] != "") {
        paste0(filtering_colname, " != ", filtering_term, " & ", df[["Input_filtering"]][i])
      } else {
        paste0(filtering_colname, " != ", filtering_term)
      }
      filter(data_df, eval(parse(text = filtering)))[["seq_name"]]
    }
    target <- data_df[[df[["Target_name"]][i]]][which(data_df[["seq_name"]] %in% selected)]
    train_dat <- ngram_matrix[which(data_df[["seq_name"]] %in% selected), ]
    imp_ngrams <- if(df[["Multiclass"]][i] == TRUE) {
      unique(
        unlist(
          unname(
            get_imp_ngrams_mc(train_dat, 
                              data_df[which(data_df[["seq_name"]] %in% selected), ],
                              df[["Target_name"]][i], cutoff = 0.01))))
    } else{
      calc_imp_ngrams(train_dat, as.logical(target), cutoff = 0.01)
    }
    if(grepl("SMOTE", df[["Model_name"]][i])) {
      train_rf(train_dat, target, imp_ngrams, with_class_weights = FALSE, smote = TRUE)
    } else {
      train_rf(train_dat, target, imp_ngrams, with_class_weights = FALSE, smote = FALSE)
    }
    
  }) %>% setNames(df[["Model_name"]])
}

train_profile_HMM_models <- function(model_df, sequences, data_df, filtering_colname, filtering_term) {
  df <- filter(model_df, Input_data == "sequences")
  lapply(1:nrow(df), function(i) {
    train_seq <- if(is.null(filtering_colname) & is.null(filtering_term)) {
      if(df[["Input_filtering"]][i] != "") {
        filtering <- paste0(df[["Input_filtering"]][i])
        sequences[names(sequences) %in% filter(data_df, eval(parse(text = filtering)))[["seq_name"]]]
      } else {
        sequences
      }
    } else { 
      filtering <- if(df[["Input_filtering"]][i] != "") {
        paste0(filtering_colname, " != ", filtering_term, " & ", df[["Input_filtering"]][i])
      } else {
        paste0(filtering_colname, " != ", filtering_term)
      }
      sequences[names(sequences) %in% filter(data_df, eval(parse(text = filtering)))[["seq_name"]]]
    }
    
    train_profileHMM(train_seq, df[["Model_name"]][i])
    df[["Model_name"]][i]
  }) %>% setNames(df[["Model_name"]])
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


predict_with_all_models <- function(model_df, test_ngram_matrix, test_df, test_sequences, ngram_models, 
                                    hmm_models = list("Sec_model" = "Sec_model", "Tat_model" = "Tat_model"), ith_fold, remove_hmm_files = TRUE) {
  all_preds <- lapply(1:nrow(model_df), function(i) {
    results <- if(model_df[["Input_data"]][i] == "ngrams") {
      if(model_df[["Multiclass"]][i] == TRUE) {
        x <- cbind(data.frame(seq_name = test_df[["seq_name"]]),
                   predict(ngram_models[[model_df[["Model_name"]][i]]], test_ngram_matrix)[["predictions"]])
        colnames(x) <- sapply(colnames(x), function(j) ifelse(j == "seq_name", "seq_name", paste0(model_df[["Model_name"]][i], "_", j)))
        x
      } else {
        data.frame(seq_name = test_df[["seq_name"]],
                   pred = predict(ngram_models[[model_df[["Model_name"]][i]]], test_ngram_matrix)[["predictions"]][, "TRUE"]) %>% 
          setNames(c("seq_name", model_df[["Model_name"]][i]))
      }
      
    } else if(model_df[["Input_data"]][i] == "sequences") {
      predict_profileHMM(test_sequences, hmm_models[model_df[["Model_name"]][i]], remove_files = remove_hmm_files) %>% 
        mutate(pred = (2^sequence_score) / (1+2^sequence_score)) %>% 
        select(c("domain_name", "pred")) %>% 
        setNames(c("seq_name", model_df[["Model_name"]][i]))
    }
    #write.csv(results, paste0(model_df[["Model_name"]][i], "_predictions_fold", ith_fold, ".csv"), row.names = FALSE)
  }) %>% reduce(., full_join, by = "seq_name")
  
  if(is.null(ith_fold)) {
    all_preds
  } else {
    mutate(all_preds, fold = ith_fold)
  }
}

get_all_models_predictions <- function(ngram_matrix, sequences, data_df, model_df, data_path, remove_hmm_files = FALSE) {
  
  res <- lapply(unique(data_df[["fold"]]), function(ith_fold) {
    dat <- ngram_matrix[data_df[["fold"]] != ith_fold, ]
    test_dat <- ngram_matrix[data_df[["fold"]] == ith_fold, ]
    test_df <- filter(data_df, fold == ith_fold)
    test_seqs <- sequences[which(names(sequences) %in% test_df[["seq_name"]])]
    
    ngram_models <- train_ngram_models(model_df, ngram_matrix, data_df, filtering_colname = "fold", filtering_term =  ith_fold)
    hmm_models <- train_profile_HMM_models(model_df, sequences, data_df, filtering_colname = "fold", filtering_term = ith_fold)
    
    predict_with_all_models(model_df, test_dat, test_df, test_seqs, ngram_models, hmm_models, ith_fold, remove_hmm_files = remove_hmm_files)
  }) %>% bind_rows()
  write.csv(res, paste0(data_path, "All_models_predictions.csv"), row.names = FALSE)
  res
} 


filter_results_for_single_architecture <- function(architecture_file, all_models_results) {
  all_models_results[is.na(all_models_results)] <- 0
  arch <- read.csv(architecture_file)
  res_colnames <- unlist(sapply(arch[["Model_name"]], function(ith_name) {
    if(grepl("Nuclear_membrane_model", ith_name)) 
    {paste0(ith_name, c("_OM", "_IM", "_TM"))
    } else {
      ith_name
    }
  }, USE.NAMES = FALSE))
  res <- lapply(1:nrow(arch), function(i) {
    if(arch[["Multiclass"]][i] == FALSE) {
      if(is.na(arch[["Test_filtering"]][i]) | arch[["Test_filtering"]][i] == "") {
        if(arch[["SMOTE"]][i] == FALSE) {
          select(all_models_results, c(seq_name, fold, arch[["Model_name"]][i]))
        } else {
          select(all_models_results, c(seq_name, fold, paste0(arch[["Model_name"]][i], "_SMOTE")))
        }
      } else {
        filtering <- arch[["Test_filtering"]][i]
        if(arch[["SMOTE"]][i] == FALSE) {
          select(filter(all_models_results, eval(parse(text = filtering))), c(seq_name, fold, arch[["Model_name"]][i]))
        } else {
          select(filter(all_models_results, eval(parse(text = filtering))), c(seq_name, fold, paste0(arch[["Model_name"]][i], "_SMOTE")))
        }
      }
    } else {
      if(is.na(arch[["Test_filtering"]][i]) | arch[["Test_filtering"]][i] == "") {
        if(arch[["SMOTE"]][i] == FALSE) {
          select(all_models_results, c(seq_name, fold, Nuclear_membrane_model_OM, Nuclear_membrane_model_IM, Nuclear_membrane_model_TM))
        } else {
          select(all_models_results, c(seq_name, fold, Nuclear_membrane_model_SMOTE_OM, Nuclear_membrane_model_SMOTE_IM, Nuclear_membrane_model_SMOTE_TM))
        }
      } else {
        filtering <- arch[["Test_filtering"]][i]
        if(arch[["SMOTE"]][i] == FALSE) {
          select(filter(all_models_results, eval(parse(text = filtering))), 
                 c(seq_name, fold, Nuclear_membrane_model_OM, Nuclear_membrane_model_IM, Nuclear_membrane_model_TM))
        } else {
          select(filter(all_models_results, eval(parse(text = filtering))),  
                 c(seq_name, fold, Nuclear_membrane_model_SMOTE_OM, Nuclear_membrane_model_SMOTE_IM, Nuclear_membrane_model_SMOTE_TM))
        }
      }
    }
  }) %>% reduce(., left_join, by = c("seq_name", "fold"))
  res[is.na(res)] <- 0
  res
}


generate_results_for_architectures <- function(architecture_file_list, all_models_results, outdir, data_df) {
  lapply(architecture_file_list, function(ith_file) {
    model <- gsub(".csv", "", last(strsplit(ith_file, "/")[[1]]))
    res <- left_join(data_df[, c("seq_name", "dataset")], 
                     filter_results_for_single_architecture(ith_file, all_models_results),
                     by = "seq_name")
    full_results <- lapply(unique(res[["fold"]]), function(ith_fold) {
      train_dat <- select(filter(res, fold != ith_fold), -seq_name)
      test_dat <- filter(res, fold == ith_fold)
      lm_model <- multinom(dataset ~ ., train_dat, model = TRUE)
      preds <- test_dat %>%
        select(c("dataset", "fold")) %>%
        mutate(Prediction = predict(lm_model, test_dat))
      probs <- predict(lm_model, test_dat, type = "probs") %>% 
        as.data.frame() %>% 
        mutate(Probability = sapply(1:nrow(.), function(i) max(.[i,])))
      cbind(preds, probs) %>% 
        mutate(fold = ith_fold)
    }) %>% bind_rows()
    write.csv(full_results, paste0(outdir, model, "_results.csv"), row.names = FALSE)
  })
}

evaluate_all_architectures <- function(res_files, outfile) {
  performance_results <- lapply(res_files, function(ith_file) {
    res <- read.csv(ith_file)
    lapply(unique(res[["fold"]]), function(ith_fold) {
      dat <- filter(res, fold == ith_fold)
      data.frame(
        model = gsub("_results.csv", "", last(strsplit(ith_file, "/")[[1]])),
        fold = ith_fold,
        AU1U = multiclass.AU1U(dat[, 4:(ncol(dat)-1)], dat[["dataset"]]),
        kappa = KAPPA(dat[["dataset"]], dat[["Prediction"]]),
        N_IM_sensitivity = TPR(ifelse(dat[["dataset"]] == "N_IM", TRUE, FALSE),
                               ifelse(dat[["Prediction"]] == "N_IM", TRUE, FALSE), TRUE),
        N_OM_sensitivity = TPR(ifelse(dat[["dataset"]] == "N_OM", TRUE, FALSE),
                               ifelse(dat[["Prediction"]] == "N_OM", TRUE, FALSE), TRUE),
        N_TM_sensitivity = TPR(ifelse(dat[["dataset"]] == "N_TM", TRUE, FALSE),
                               ifelse(dat[["Prediction"]] == "N_TM", TRUE, FALSE), TRUE),
        N_S_sensitivity = TPR(ifelse(dat[["dataset"]] == "N_S", TRUE, FALSE),
                              ifelse(dat[["Prediction"]] == "N_S", TRUE, FALSE), TRUE),
        N_TL_SEC_sensitivity = TPR(ifelse(dat[["dataset"]] == "N_TL_SEC", TRUE, FALSE),
                                   ifelse(dat[["Prediction"]] == "N_TL_SEC", TRUE, FALSE), TRUE),
        N_TL_TAT_sensitivity = TPR(ifelse(dat[["dataset"]] == "N_TL_TAT", TRUE, FALSE),
                                   ifelse(dat[["Prediction"]] == "N_TL_TAT", TRUE, FALSE), TRUE),
        P_IM_sensitivity = TPR(ifelse(dat[["dataset"]] == "P_IM", TRUE, FALSE),
                               ifelse(dat[["Prediction"]] == "P_IM", TRUE, FALSE), TRUE),
        P_TM_sensitivity = TPR(ifelse(dat[["dataset"]] == "P_TM", TRUE, FALSE),
                               ifelse(dat[["Prediction"]] == "P_TM", TRUE, FALSE), TRUE),
        P_S_sensitivity = TPR(ifelse(dat[["dataset"]] == "P_S", TRUE, FALSE),
                              ifelse(dat[["Prediction"]] == "P_S", TRUE, FALSE), TRUE)
      )
    }) %>% bind_rows()
  }) %>% bind_rows()
  write.csv(performance_results, outfile, row.names = FALSE)
  performance_results
}


get_mean_performance_of_architectures <- function(performance_results, outfile) {
  summ <- performance_results %>% 
    group_by(model) %>% 
    summarise(mean_AU1U = mean(AU1U),
              mean_kappa = mean(kappa),
              mean_N_IM_sensitivity = mean(N_IM_sensitivity),
              mean_N_OM_sensitivity = mean(N_OM_sensitivity),
              mean_N_TM_sensitivity = mean(N_TM_sensitivity),
              mean_N_S_sensitivity = mean(N_S_sensitivity),
              mean_N_TL_SEC_sensitivity = mean(N_TL_SEC_sensitivity),
              mean_N_TL_TAT_sensitivity = mean(N_TL_TAT_sensitivity),
              mean_P_IM_sensitivity = mean(P_IM_sensitivity),
              mean_P_TM_sensitivity = mean(P_TM_sensitivity),
              mean_P_S_sensitivity = mean(P_S_sensitivity))
  write.csv(summ, outfile, row.names = FALSE)
  summ
}

