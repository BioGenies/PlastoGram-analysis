library(UBL)
library(dplyr)
library(targets)
library(seqR)
library(biogram)
library(ranger)
library(tidyr)
library(measures)
library(ggplot2)

tar_load(c(ngram_matrix, target_df, data_df))

do_binary_smote_cv <- function(ngram_matrix, data_df, cutoff, smote = TRUE, C.perc = "balance", k = k, name, rep) {
  lapply(unique(data_df[["fold"]]), function(ith_fold) {
    dat <- ngram_matrix[which(data_df[["Nuclear_target"]] == TRUE & (data_df[["Membrane_mc_target"]] == 'OM' | data_df[["Stroma_target"]] == TRUE) & data_df[["fold"]] != ith_fold), ]
    tar <- as.factor(filter(data_df, Nuclear_target == TRUE & (Membrane_mc_target == 'OM' | Stroma_target == TRUE) & fold != ith_fold)[["OM_target"]])
    
    imp_ngrams <- calc_imp_ngrams(dat, as.logical(tar), cutoff)
    
    train_dat <- if(smote == FALSE) {
      x <- cbind(dat[, imp_ngrams], target = tar)
    } else {
      x <- cbind(dat[, imp_ngrams], target = tar)
      x[sapply(x, is.numeric)] <- lapply(x[sapply(x, is.numeric)], 
                                         as.factor)
      x <- SmoteClassif(target ~ ., x, C.perc = C.perc, dist = "Overlap")
      modifyList(x, lapply(x[, which(colnames(x) != "target")], function(i) as.numeric(as.character(i))))
    }
    
    rf <- ranger(data = train_dat, dependent.variable.name = "target", write.forest = TRUE,
                 probability = TRUE, num.trees = 500, verbose = FALSE, seed = 123456, classification = TRUE)
    test_dat <- ngram_matrix[which(data_df[["Nuclear_target"]] == TRUE & (data_df[["Membrane_mc_target"]] == 'OM' | data_df[["Stroma_target"]] == TRUE) & data_df[["fold"]] == ith_fold), ]
    
    filter(data_df, Nuclear_target == TRUE & (Membrane_mc_target == 'OM' | Stroma_target == TRUE) & fold == ith_fold) %>% 
      select(c(seq_name, OM_target, fold)) %>% 
      mutate(Prediction = predict(rf, test_dat[, imp_ngrams])[["predictions"]][, "TRUE"])
  }) %>% bind_rows() %>% 
    mutate(Model = name,
           rep = rep)
}

models <- list(
  "Without_SMOTE_0.001" = c(0.001, FALSE, NULL, NULL),
  "Without_SMOTE_0.01" = c(0.01, FALSE, NULL, NULL),
  "SMOTE_0.001_balance_5" = c(0.001, TRUE, "balance", 5),
  "SMOTE_0.001_balance_3" = c(0.001, TRUE, "balance", 3),
  "SMOTE_0.01_balance_5" = c(0.01, TRUE, "balance", 5),
  "SMOTE_0.01_balance_3" = c(0.01, TRUE, "balance", 3),
  "SMOTE_0.001_(5.57,1)_5" = c(0.001, TRUE, list("TRUE" = 5.57, "FALSE" = 1), 5),
  "SMOTE_0.01_(5.57,1)_5" = c(0.01, TRUE, list("TRUE" = 5.57, "FALSE" = 1), 5),
  "SMOTE_0.001_(5.57,1)_3" = c(0.001, TRUE, list("TRUE" = 5.57, "FALSE" = 1), 3),
  "SMOTE_0.01_(5.57,1)_3" = c(0.01, TRUE, list("TRUE" = 5.57, "FALSE" = 1), 3)
)

smote_res <- lapply(1:10, function(i) {
  lapply(names(models), function(ith_model) {
    do_binary_smote_cv(ngram_matrix, data_df, cutoff = models[[ith_model]][1], smote = models[[ith_model]][2], 
                       C.perc = models[[ith_model]][3], k = models[[ith_model]][4], name = ith_model,
                       rep = i)
  }) %>% bind_rows()
}) %>% bind_rows()



stats <- smote_res %>% 
  mutate(Decision = ifelse(Prediction > 0.5, TRUE, FALSE)) %>% 
  group_by(rep, Model, fold) %>% 
  summarise(TP = TP(OM_target, Decision, TRUE),
            TN = TN(OM_target, Decision, FALSE),
            FP = FP(OM_target, Decision, TRUE),
            FN = FN(OM_target, Decision, FALSE),
            Sensitivity = mlr3measures::sensitivity(as.factor(OM_target), as.factor(Decision), "TRUE"),
            Specificity = mlr3measures::specificity(as.factor(OM_target), as.factor(Decision), "TRUE"),
            Accuracy = ACC(OM_target, Decision),
            AUC = AUC(Prediction, OM_target, FALSE, TRUE))

mean_stats <- stats %>% 
  group_by(Model, rep) %>% 
  summarise(mean_accuracy = mean(Accuracy),
            mean_AUC = mean(AUC),
            mean_sensitivity = mean(Sensitivity),
            mean_specificity = mean(Specificity))


ggplot(mean_stats, aes(x = mean_sensitivity, y = mean_specificity, color = Model)) +
  geom_point(size = 3) 

pivot_longer(mean_stats, mean_accuracy:mean_specificity, names_to = "Measure", values_to = "Value") %>% 
  ggplot(aes(x = Model, y = Value, fill = Measure)) +
  geom_boxplot() +
  facet_wrap(~Measure, scales = "free_y", ncol = 1)

pivot_longer(stats, TP:AUC, names_to = "Measure", values_to = "Value") %>% 
  filter(Measure %in% c("Accuracy", "AUC", "Sensitivity", "Specificity")) %>% 
  ggplot(aes(x = Model, y = Value, fill = Measure)) +
  geom_boxplot() +
  facet_wrap(~Measure, scales = "free_y", ncol = 1)
