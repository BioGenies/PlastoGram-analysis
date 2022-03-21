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

create_cv_folds <- function(target_df, n_fold, seeds) {
  lapply(seeds, function(ith_seed) {
    set.seed(ith_seed)
    create_folds(target_df, n_fold)
  })
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
    print(paste0("Now training ", df[["Model_name"]][[i]]))
    imp_ngrams <- if(df[["Multiclass"]][i] == TRUE) {
      unique(
        unlist(
          unname(
            get_imp_ngrams_mc(train_dat, 
                              data_df[which(data_df[["seq_name"]] %in% selected), ],
                              df[["Target_name"]][i], cutoff = 0.01))))
    } else {
      calc_imp_ngrams(train_dat, as.logical(target), cutoff = 0.01)
    }
    print(paste0("Number of informative ngrams for ", df[["Model_name"]][[i]], ": ", length(imp_ngrams)))
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
    print(paste0("Training models for fold ", ith_fold))
    ngram_models <- train_ngram_models(model_df, ngram_matrix, data_df, filtering_colname = "fold", filtering_term =  ith_fold)
    hmm_models <- train_profile_HMM_models(model_df, sequences, data_df, filtering_colname = "fold", filtering_term = ith_fold)
    print(paste0("Finished training models for fold ", ith_fold))
    predict_with_all_models(model_df, test_dat, test_df, test_seqs, ngram_models, hmm_models, ith_fold, remove_hmm_files = remove_hmm_files)
  }) %>% bind_rows()
  res
} 

get_all_models_predictions_cv <- function(ngram_matrix, sequences, data_dfs_cv, model_df, data_path, remove_hmm_files = FALSE) {
  lapply(1:length(data_dfs_cv), function(i) {
    res <- get_all_models_predictions(ngram_matrix, sequences, data_dfs_cv[[i]], model_df, data_path, remove_hmm_files)
    write.csv(res, paste0(data_path, "All_models_predictions_rep", i, ".csv"), row.names = FALSE)
    res
  })
}


filter_results_for_single_architecture <- function(architecture_file, all_models_results) {
  all_models_results[is.na(all_models_results)] <- 0
  arch <- read.csv(architecture_file)
  res_colnames <- unlist(sapply(arch[["Model_name"]], function(ith_name) {
    if(grepl("Nuclear_membrane_mc_model", ith_name)) {
      paste0(ith_name, c("_OM", "_IM", "_TM"))
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
          select(all_models_results, c(seq_name, fold, Nuclear_membrane_mc_model_OM, Nuclear_membrane_mc_model_IM, Nuclear_membrane_mc_model_TM))
        } else {
          select(all_models_results, c(seq_name, fold, Nuclear_membrane_mc_model_SMOTE_OM, Nuclear_membrane_mc_model_SMOTE_IM, Nuclear_membrane_mc_model_SMOTE_TM))
        }
      } else {
        filtering <- arch[["Test_filtering"]][i]
        if(arch[["SMOTE"]][i] == FALSE) {
          select(filter(all_models_results, eval(parse(text = filtering))), 
                 c(seq_name, fold, Nuclear_membrane_mc_model_OM, Nuclear_membrane_mc_model_IM, Nuclear_membrane_mc_model_TM))
        } else {
          select(filter(all_models_results, eval(parse(text = filtering))),  
                 c(seq_name, fold, Nuclear_membrane_mc_model_SMOTE_OM, Nuclear_membrane_mc_model_SMOTE_IM, Nuclear_membrane_mc_model_SMOTE_TM))
        }
      }
    }
  }) %>% reduce(., left_join, by = c("seq_name", "fold"))
  res[is.na(res)] <- 0
  res
}


generate_results_for_architectures <- function(architecture_file_list, all_models_results, outdir, data_df, higher_level_models = c("RF", "GLM")) {
  lapply(1:length(all_models_results), function(ith_rep) {
    lapply(higher_level_models, function(ith_hl_model) {
      lapply(architecture_file_list, function(ith_file) {
        model <- gsub(".csv", "", last(strsplit(ith_file, "/")[[1]]))
        res <- left_join(data_df[, c("seq_name", "dataset")], 
                         filter_results_for_single_architecture(ith_file, all_models_results[[ith_rep]]),
                         by = "seq_name")
        full_results <- lapply(unique(res[["fold"]]), function(ith_fold) {
          train_dat <- select(filter(res, fold != ith_fold), -c(seq_name, fold))
          test_dat <- filter(res, fold == ith_fold)
          if(ith_hl_model == "GLM") {
            lm_model <- multinom(dataset ~ ., train_dat, model = TRUE)
            preds <- test_dat %>%
              select(c("seq_name", "dataset", "fold")) %>%
              mutate(Prediction = predict(lm_model, test_dat))
            probs <- predict(lm_model, test_dat, type = "probs") %>% 
              as.data.frame() %>% 
              mutate(Probability = sapply(1:nrow(.), function(i) max(.[i,])))
            cbind(preds, probs)
          } else {
            rf_model <- ranger(dataset ~ ., data = mutate(train_dat, dataset = as.factor(dataset)),
                               write.forest = TRUE, probability = TRUE, num.trees = 500, 
                               verbose = FALSE, seed = 427244)
            preds <- predict(rf_model, test_dat)[["predictions"]]
            cbind(select(test_dat, c("seq_name", "dataset", "fold")),
                  Prediction = c(colnames(preds)[max.col(preds[, c(colnames(preds))])]),
                  preds,
                  Probability = sapply(1:nrow(preds), function(i) max(preds[i,])))
          }
        }) %>% bind_rows()
        write.csv(full_results, paste0(outdir, model, "_", ith_hl_model, "_rep", ith_rep, "_results.csv"), row.names = FALSE)
      })
    })
  })
}

evaluate_all_architectures <- function(res_files, outfile, data_df) {
  performance_results <- lapply(res_files, function(ith_file) {
    res <- read.csv(ith_file)
    x <- table(data_df[["dataset"]])
    lapply(unique(res[["fold"]]), function(ith_fold) {
      dat <- filter(res, fold == ith_fold)
      data.frame(
        model = gsub("_rep._results.csv|_rep.._results.csv", "", last(strsplit(ith_file, "/")[[1]])),
        rep = last(strsplit(gsub("_results.csv", "", last(strsplit(ith_file, "/")[[1]])), "_")[[1]]),
        fold = ith_fold,
        AU1U = multiclass.AU1U(dat[, 5:(ncol(dat)-1)], dat[["dataset"]]),
        kappa = KAPPA(dat[["dataset"]], dat[["Prediction"]]),
        weighted_kappa = ckap(dat[, c("dataset", "Prediction")], 
                              weight = as.matrix(bind_rows(lapply(x, function(i) ifelse(i == x, as.vector(sum(x)/(i+x)), 1)))))[["est"]],
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
    group_by(model, rep) %>% 
    summarise(mean_rep_kappa = mean(kappa),
              mean_rep_AU1U = mean(AU1U),
              mean_rep_N_IM_sens = mean(N_IM_sensitivity),
              mean_rep_N_OM_sens = mean(N_OM_sensitivity),
              mean_rep_N_TM_sens = mean(N_TM_sensitivity),
              mean_rep_N_S_sens = mean(N_S_sensitivity),
              mean_rep_N_TL_SEC_sens = mean(N_TL_SEC_sensitivity),
              mean_rep_N_TL_TAT_sens = mean(N_TL_TAT_sensitivity),
              mean_rep_P_IM_sens = mean(P_IM_sensitivity),
              mean_rep_P_TM_sens = mean(P_TM_sensitivity),
              mean_rep_P_S_sens = mean(P_S_sensitivity),
              sd_rep_N_IM_sens = sd(N_IM_sensitivity),
              sd_rep_N_OM_sens = sd(N_OM_sensitivity),
              sd_rep_N_TM_sens = sd(N_TM_sensitivity),
              sd_rep_N_S_sens = sd(N_S_sensitivity),
              sd_rep_N_TL_SEC_sens = sd(N_TL_SEC_sensitivity),
              sd_rep_N_TL_TAT_sens = sd(N_TL_TAT_sensitivity),
              sd_rep_P_IM_sens = sd(P_IM_sensitivity),
              sd_rep_P_TM_sens = sd(P_TM_sensitivity),
              sd_rep_P_S_sens = sd(P_S_sensitivity)) %>% 
    ungroup() %>% 
    group_by(model) %>% 
    summarise(mean_kappa = mean(mean_rep_kappa),
              mean_AU1U = mean(mean_rep_AU1U),
              mean_N_IM_sens = mean(mean_rep_N_IM_sens),
              mean_N_OM_sens = mean(mean_rep_N_OM_sens),
              mean_N_TM_sens = mean(mean_rep_N_TM_sens),
              mean_N_S_sens = mean(mean_rep_N_S_sens),
              mean_N_TL_SEC_sens = mean(mean_rep_N_TL_SEC_sens),
              mean_N_TL_TAT_sens = mean(mean_rep_N_TL_TAT_sens),
              mean_P_IM_sens = mean(mean_rep_P_IM_sens),
              mean_P_TM_sens = mean(mean_rep_P_TM_sens),
              mean_P_S_sens = mean(mean_rep_P_S_sens),
              sd_N_IM_sens = mean(sd_rep_N_IM_sens),
              sd_N_OM_sens = mean(sd_rep_N_OM_sens),
              sd_N_TM_sens = mean(sd_rep_N_TM_sens),
              sd_N_S_sens = mean(sd_rep_N_S_sens),
              sd_N_TL_SEC_sens = mean(sd_rep_N_TL_SEC_sens),
              sd_N_TL_TAT_sens = mean(sd_rep_N_TL_TAT_sens),
              sd_P_IM_sens = mean(sd_rep_P_IM_sens),
              sd_P_TM_sens = mean(sd_rep_P_TM_sens),
              sd_P_S_sens = mean(sd_rep_P_S_sens)) 
  write.csv(summ, outfile, row.names = FALSE)
  summ
}

rank_architectures <- function(mean_architecture_performance) {
  ranking_columns <- colnames(mean_architecture_performance)[which(!(colnames(mean_architecture_performance) %in% c("model", "mean_AU1U")))]
  ranks_df <- lapply(ranking_columns, function(ith_col) {
    ranks <- if(grepl("mean", ith_col)) {
      rank(-mean_architecture_performance[[ith_col]], ties.method = "min")
    } else {
      rank(mean_architecture_performance[[ith_col]], ties.method = "min")
    } 
    df <- data.frame("ranks" = ranks)
    colnames(df) <- paste0(ith_col, "_ranks")
    df
  }) %>% bind_cols() 
  cbind(data.frame(model = mean_architecture_performance[["model"]]), 
        ranks_df) %>% 
    mutate(rank_sum = rowSums(ranks_df))
}

train_multinom <- function(prediction_results, data_df) {
  train_dat <- left_join(prediction_results, data_df[, c("seq_name", "dataset")], by = "seq_name") %>% 
    select(-seq_name)
  train_dat[is.na(train_dat)] <- 0
  multinom(dataset ~ ., data = train_dat, model = TRUE)
}

scaled_train_multinom <- function(train_df, data_df) {
  train_df[is.na(train_df)] <- 0
  scaled_train_df <- scale(select(train_df, -seq_name))
  glm_scale_scale <- attr(scaled_train_df, "scaled:scale")
  glm_scale_center <- attr(scaled_train_df, "scaled:center")
  scaled_train_df <- data.frame(scaled_train_df)
  scaled_train_df[["dataset"]] <- data_df[["dataset"]]
  
  hl_model <- invisible(multinom(dataset ~ ., data = scaled_train_df, model = TRUE))
  
  list("model" = hl_model,
       "scale" = glm_scale_scale,
       "center" = glm_scale_center)
 
}

train_higher_level_model <- function(PlastoGram_predictions, data_df, best_architecture_name) {
  if(grepl("GLM", best_architecture_name)) {
    train_multinom(PlastoGram_predictions, data_df)
  } else {
    train_dat <- left_join(PlastoGram_predictions, data_df[, c("seq_name", "dataset")], by = "seq_name") %>% 
      select(-seq_name)
    train_dat[is.na(train_dat)] <- 0
    ranger(dataset ~ ., data = mutate(train_dat, dataset = as.factor(dataset)),
           write.forest = TRUE, probability = TRUE, num.trees = 500, 
           verbose = FALSE, seed = 427244)
  }
}


generate_and_test_architectures <- function(model_variants, smote_models, sequence_models, model_dat, filtering_df, architectures_output_dir,
                                            all_models_predictions, architecture_res_output_dir, data_df_final, performance_outfile) {
  
  generate_all_architectures(model_variants = model_variants,
                             smote_models = smote_models,
                             sequence_models = sequence_models,
                             model_dat = model_dat,
                             filtering_df = filtering_df,
                             output_dir = architectures_output_dir)
  
  architecture_files <- list.files(architectures_output_dir, full.names = TRUE)
  
  architecture_results <- generate_results_for_architectures(architecture_files,
                                                             all_models_predictions,
                                                             architecture_res_output_dir,
                                                             data_df_final)
  
  architecture_results_files <- list.files(architecture_res_output_dir, full.names = TRUE)
  
  evaluate_all_architectures(architecture_results_files,
                             performance_outfile,
                             data_df_final)
}
