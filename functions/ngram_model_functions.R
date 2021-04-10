create_ngram_matrix <- function(seqs, target_df) {
  ngrams <- lapply(target_df[["seq_name"]], function(ith_seq) {
    count_multimers(seqs[[ith_seq]],
                    k_vector = c(1, rep(2,2), rep(3,4)),
                    kmer_gaps_list = list(NULL, NULL, 1, c(0,0), c(0,1), c(1,0), c(1,1)),
                    alphabet = toupper(colnames(biogram::aaprop))) %>% 
      binarize() %>% 
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


train_rf <- function(ngram_matrix, target, imp_ngrams) {
  ranger_train_data <- data.frame(ngram_matrix[, imp_ngrams],
                                  tar = target)
  
  ranger(dependent.variable.name = "tar", data =  ranger_train_data, 
         write.forest = TRUE, probability = TRUE, num.trees = 500, 
         verbose = FALSE, classification = TRUE)
}


do_cv <- function(ngram_matrix, target_df, target_col, n_fold, cutoff) {
  
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
    test_dat <- ngram_matrix[data_df[["fold"]] == ith_fold, ]
    imp_ngrams <- calc_imp_ngrams(dat, data_df[[target_col]][data_df[["fold"]] != ith_fold], cutoff)
    trained_model <- train_rf(dat, data_df[[target_col]][data_df[["fold"]] != ith_fold], imp_ngrams)
    data_df %>% 
      filter(fold == ith_fold) %>% 
      select(seq_name, fold, target_col) %>% 
      setNames(c("seq_name", "fold", "target")) %>% 
      mutate(prob = predict(trained_model, test_dat)[["predictions"]][, 1],
             pred = ifelse(prob > 0.5, 1, 0))
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
