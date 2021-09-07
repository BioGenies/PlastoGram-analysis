create_ngram_matrix <- function(seqs, target_df) {
  ngrams <- lapply(target_df[["seq_name"]], function(ith_seq) {
    count_multimers(seqs[ith_seq],
                    k_vector = c(1, rep(2,4), rep(3,4)),
                    kmer_gaps_list = list(NULL, NULL, 1, 2, 3, c(0,0), c(0,1), c(1,0), c(1,1)),
                    alphabet = toupper(colnames(biogram::aaprop)),
                    with_kmer_counts = FALSE) %>% 
      as.matrix() %>% 
      as.data.frame()
  })  %>% bind_rows() 
  
  ngrams[is.na(ngrams)] <- 0
  ngrams
}


calc_imp_ngrams <- function(ngram_matrix, target, cutoff = 0.001) {
  test_bis <- test_features(target = target, 
                            features = ngram_matrix)
  cut(test_bis, breaks = c(0, cutoff, 1))[[1]]
}


train_rf <- function(ngram_matrix, target, imp_ngrams, with_class_weights = FALSE) {
  ranger_train_data <- data.frame(ngram_matrix[, imp_ngrams],
                                  tar = as.factor(target))
  if(with_class_weights == FALSE) {
    class_weights <- NULL
  } else {
    class_weights <- sapply(levels(ranger_train_data[["tar"]]), function(i) 1/sum(ranger_train_data[["tar"]] == i)/nrow(ranger_train_data), USE.NAMES = FALSE)
  }
  ranger(dependent.variable.name = "tar", data =  ranger_train_data, 
         write.forest = TRUE, probability = TRUE, num.trees = 500, 
         verbose = FALSE, class.weights = class_weights)
}


do_cv <- function(ngram_matrix, data_df, target_df, target_col, n_fold, cutoff, mc = FALSE, with_class_weights = FALSE, fcbf = FALSE) {
  lapply(1:n_fold, function(ith_fold) {
    dat <- ngram_matrix[data_df[["fold"]] != ith_fold, ]
    test_dat <- ngram_matrix[data_df[["fold"]] == ith_fold, ]
    imp_ngrams <- if(mc == TRUE) {
      if(fcbf == TRUE) {
        do_fcbf(dat, data_df[[target_col]][which(data_df[["fold"]] != ith_fold)], min_su = cutoff, multiclass = TRUE)
      } else {
        imp <- get_imp_ngrams_mc(dat, data_df[which(data_df[["fold"]] != ith_fold), ], target_col, cutoff)
        unique(unlist(unname(imp)))
      }
    } else {
      if(fcbf == TRUE) {
        do_fcbf(dat, data_df[[target_col]][which(data_df[["fold"]] != ith_fold)], min_su = cutoff)
      } else {
        calc_imp_ngrams(dat, data_df[[target_col]][data_df[["fold"]] != ith_fold], cutoff) 
      }
    }
    trained_model <- train_rf(dat, data_df[[target_col]][data_df[["fold"]] != ith_fold], imp_ngrams, with_class_weights)
    
    full_res <- if(mc == TRUE) {
      classes <- unique(target_df[[target_col]])
      res <- data_df %>% 
        filter(fold == ith_fold) %>% 
        select(seq_name, fold, target_col) %>% 
        setNames(c("seq_name", "fold", "target")) %>% 
        bind_cols(as.data.frame(predict(trained_model, test_dat)[["predictions"]]))
      mutate(res, pred = c(classes)[max.col(res[, c(classes)])],
             imp_ngrams = length(imp_ngrams))
    } else {
      data_df %>% 
        filter(fold == ith_fold) %>% 
        select(seq_name, fold, target_col) %>% 
        setNames(c("seq_name", "fold", "target")) %>% 
        mutate(prob = predict(trained_model, test_dat)[["predictions"]][, "TRUE"],
               pred = ifelse(prob > 0.5, TRUE, FALSE),
               imp_ngrams = length(imp_ngrams))
    }
    full_res
  }) %>% bind_rows()
} 


get_cv_res_summary <- function(cv_res, positive, ngrams = TRUE) {
  if(ngrams == FALSE) {
    cv_res %>% 
      group_by(fold) 
  } else {
    cv_res %>% 
      group_by(fold, imp_ngrams) 
  } %>% 
    summarise(TP = mlr3measures::tp(factor(target), factor(pred), positive),
              TN = mlr3measures::tn(factor(target), factor(pred), positive),
              FP = mlr3measures::fp(factor(target), factor(pred), positive),
              FN = mlr3measures::fn(factor(target), factor(pred), positive),
              AUC = mlr3measures::auc(factor(target), prob, positive), 
              ACC = mlr3measures::acc(factor(target), factor(pred)))
} 


get_imp_ngrams_mc <- function(ngram_matrix, target_df, target_col, cutoff = 0.001) {
  combns <- combn(unique(target_df[[target_col]]), 2, simplify = FALSE)
  lapply(combns, function(ith_cmbn) {
    tar <- filter(target_df, get(target_col) %in% ith_cmbn) %>% 
      mutate(target = ifelse(get(target_col) == ith_cmbn[1], 1, 0))
    features <- ngram_matrix[target_df[[target_col]] %in% ith_cmbn, ]
    test_bis <- test_features(tar[["target"]], features)
    res <- cut(test_bis, breaks = c(0, cutoff, 1))[[1]]
  }) %>% setNames(sapply(combns, function(ith_combn) paste(ith_combn, collapse = "_")))
}


get_cv_res_summary_mc <- function(mc_cv_res) {
  classes <- unique(mc_cv_res[["target"]])
  lapply(unique(mc_cv_res[["fold"]]), function(ith_fold) {
    dat <- filter(mc_cv_res, fold == ith_fold)
    data.frame(
      fold = ith_fold,
      imp_ngrams = dat[["imp_ngrams"]][1],
      Accuracy = ACC(dat[["target"]], dat[["pred"]]),
      AU1U = multiclass.AU1U(dat[, c(classes)], dat[["target"]]),
      KapS = KAPPA(dat[["target"]], dat[["pred"]])) 
  }) %>% bind_rows()
}


calc_case_weights <- function(target) {
  lvls <- levels(as.factor(target))
  weights <- sapply(lvls, function(i) 1/(sum(target == i)/length(target)))
  sapply(target, function(i) weights[which(names(weights) == i)], USE.NAMES = FALSE)
}


do_fcbf <- function(ngrams, target, multiclass = FALSE, min_su = 0.01) {
  ngrams[sapply(ngrams, is.numeric)] <- lapply(ngrams[sapply(ngrams, is.numeric)], 
                                               as.factor)
  if(multiclass == TRUE) {
    combns <- combn(unique(target), 2, simplify = FALSE)
    lapply(combns, function(ith_cmbn) {
      tar <- sapply(target[which(target %in% ith_cmbn)], function(i) ifelse(i == ith_cmbn[1], 1, 0))
      features <- ngrams[which(target %in% ith_cmbn), ]
      res <- fcbf(features, as.factor(tar), min_su, samples_in_rows = TRUE)
      row.names(res)
    }) %>% unlist() %>% 
      unique()
  } else {
    res <- fcbf(ngrams, as.factor(target), min_su, samples_in_rows = TRUE)
    row.names(res)
  }
}


test_fcbf <- function(ngram_matrix, data_df, target_df, cutoff_vec, target_col, class_weights = FALSE, mc = FALSE) {
  lapply(cutoff_vec, function(ith_cutoff) {
    res <- do_cv(ngram_matrix, data_df, target_df, target_col, 5, cutoff = ith_cutoff, 
                 mc = mc, with_class_weights = class_weights, fcbf = TRUE)
    stats <- if(mc == TRUE) {
      get_cv_res_summary_mc(res)
    } else {
      get_cv_res_summary(res, "TRUE")
    }
   mutate(stats, cutoff = ith_cutoff)
  }) %>% bind_rows()
}
