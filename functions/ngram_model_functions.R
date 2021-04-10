create_ngram_matrix <- function(dataset_list) {
  ngrams <- lapply(names(dataset_list), function(ith_set) {
    seqs <- dataset_list[[ith_set]]
    lapply(seq_along(seqs), function(i) {
      count_multimers(seqs[[i]],
                      k_vector = c(1, rep(2,2), rep(3,4)),
                      kmer_gaps_list = list(NULL, NULL, 1, c(0,0), c(0,1), c(1,0), c(1,1)),
                      alphabet = toupper(colnames(biogram::aaprop))) %>% 
        binarize() %>% 
        as.matrix() %>% 
        as.data.frame() %>% 
        mutate(seq_name = names(seqs[i]),
               dataset = ith_set)
    })  %>% bind_rows() 
  }) %>% bind_rows()
  
  ngrams[is.na(ngrams)] <- 0
  ngrams
  }


calc_imp_ngrams <- function(ngram_matrix, target1, cutoff = 0.001) {
  tar <- ifelse(ngram_matrix[["dataset"]] == target1, 1,0)
  test_bis <- test_features(target = tar, 
                            features = select(ngram_matrix, -c("dataset", "seq_name")))
  cut(test_bis, breaks = c(0, cutoff, 1))[[1]]
}


train_rf <- function(ngram_matrix, imp_ngrams) {
  ranger_train_data <- data.frame(ngram_matrix[, imp_ngrams],
                                  tar = as.factor(ngram_matrix[["dataset"]]))
  
  ranger(dependent.variable.name = "tar", data =  ranger_train_data, 
         write.forest = TRUE, probability = TRUE, num.trees = 500, 
         verbose = FALSE)
}


do_cv <- function(dataset_list, K, target1, cutoff) {
  
  folds <- lapply(names(dataset_list), function(ith_set) {
    folded <- cvFolds(length(dataset_list[[ith_set]]), K = 5)
    fold_df <- data.frame(seq_name = names(dataset_list[[ith_set]])[folded[["subsets"]]], 
                          fold = folded[["which"]],
                          stringsAsFactors = FALSE)
  }) %>% bind_rows()
  
  ngram_matrix <- left_join(create_ngram_matrix(dataset_list),
                            folds)
  
  lapply(1:K, function(ith_fold) {
    dat <- filter(ngram_matrix, fold != ith_fold)
    test_dat <- filter(ngram_matrix, fold == ith_fold)
    imp_ngrams <- calc_imp_ngrams(dat, target1, cutoff)
    trained_model <- train_rf(dat, imp_ngrams)
    mutate(select(test_dat, c("seq_name", "dataset", "fold")), 
           prob = predict(trained_model, select(test_dat, -c("seq_name", "dataset", "fold")))[["predictions"]][, 1],
           pred = ifelse(prob > 0.5, trained_model[["forest"]][["levels"]][1], trained_model[["forest"]][["levels"]][2]))
  }) %>% bind_rows()
} 


get_cv_res_summary <- function(cv_res, positive) {
  cv_res %>% 
    group_by(fold) %>% 
    summarise(TP = mlr3measures::tp(factor(dataset), factor(pred), positive),
              TN = mlr3measures::tn(factor(dataset), factor(pred), positive),
              FP = mlr3measures::fp(factor(dataset), factor(pred), positive),
              FN = mlr3measures::fn(factor(dataset), factor(pred), positive),
              AUC = mlr3measures::auc(factor(dataset), prob, positive), 
              ACC = mlr3measures::acc(factor(dataset), factor(pred)))
} 
