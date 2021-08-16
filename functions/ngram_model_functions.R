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


do_cv <- function(ngram_matrix, target_df, target_col, n_fold, cutoff, mc = FALSE, with_class_weights = FALSE) {
  
  folds <- lapply(unique(target_df[[target_col]]), function(ith_target) {
    selected <- filter(target_df, get(target_col) == ith_target)
    folded <- cvFolds(nrow(selected), K = n_fold)
    fold_df <- data.frame(seq_name = selected[["seq_name"]][folded[["subsets"]]], 
                          fold = folded[["which"]],
                          stringsAsFactors = FALSE)
  }) %>% bind_rows()
  
  data_df <- left_join(target_df, folds, by = "seq_name")
  
  lapply(1:n_fold, function(ith_fold) {
    dat <- ngram_matrix[data_df[["fold"]] != ith_fold, ]
   # filtered_cw <- case_weights[which(data_df[["fold"]] != ith_fold)]
    test_dat <- ngram_matrix[data_df[["fold"]] == ith_fold, ]
    if(mc == TRUE) {
      imp_ngrams <- get_imp_ngrams_mc(dat, data_df[which(data_df[["fold"]] != ith_fold), ], target_col, cutoff)
      imp_ngrams <- unique(unlist(unname(imp_ngrams)))
    } else {
      imp_ngrams <- calc_imp_ngrams(dat, data_df[[target_col]][data_df[["fold"]] != ith_fold], cutoff)
    }
    trained_model <- train_rf(dat, data_df[[target_col]][data_df[["fold"]] != ith_fold], imp_ngrams, with_class_weights)
    
    if(mc == TRUE) {
      classes <- unique(target_df[[target_col]])
      res <- data_df %>% 
        filter(fold == ith_fold) %>% 
        select(seq_name, fold, target_col) %>% 
        setNames(c("seq_name", "fold", "target")) %>% 
        bind_cols(as.data.frame(predict(trained_model, test_dat)[["predictions"]]))
      full_res <- mutate(res, pred = c(classes)[max.col(res[, c(classes)])])
    } else {
      full_res <- data_df %>% 
        filter(fold == ith_fold) %>% 
        select(seq_name, fold, target_col) %>% 
        setNames(c("seq_name", "fold", "target")) %>% 
        mutate(prob = predict(trained_model, test_dat)[["predictions"]][, "TRUE"],
               pred = ifelse(prob > 0.5, TRUE, FALSE))
    }
    full_res
  }) %>% bind_rows()
} 


get_cv_res_summary <- function(cv_res, positive) {
  cv_res %>% 
    group_by(fold) %>% 
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
