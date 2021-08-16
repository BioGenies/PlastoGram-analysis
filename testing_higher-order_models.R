library(ranger)
library(nnet)
library(dplyr)
library(measures)

all_models_results <- readRDS("~/RProjects/PlastoGram-analysis/data/all_models_results.rds")

dat <- select(all_models_results, c("dataset", "fold", "Nuclear", "Membrane", "Tat", "Sec", "IM", "OM", "TM", "P_IM"))
dat[is.na(dat)] <- 0

lm_res <- lapply(unique(dat[["fold"]]), function(ith_fold) {
  train_dat <- filter(dat, fold != ith_fold) %>%
    select(-fold)
  test_dat <- filter(dat, fold == ith_fold)
  lm_model <- multinom(dataset ~ ., dat, model = TRUE)
  test_dat %>%
    select(c("dataset", "fold")) %>%
    mutate(Prediction = predict(lm_model, test_dat),
           Probability = predict(lm_model, test_dat, type = "probs"))
}) %>% bind_rows()


rf_res <- lapply(unique(dat[["fold"]]), function(ith_fold) {
  train_dat <- filter(dat, fold != ith_fold) %>%
    select(-fold) %>% 
    mutate(dataset = as.factor(dataset))
  test_dat <- filter(dat, fold == ith_fold)
  rf_model <- ranger(data = train_dat, dependent.variable.name = "dataset",
               write.forest = TRUE, probability = TRUE, verbose = FALSE)
  res <- predict(rf_model, test_dat)[["predictions"]]

  test_dat %>%
    select(c("dataset", "fold")) %>%
    mutate(Prediction = colnames(res)[max.col(res)],
           Probability = sapply(1:nrow(res), function(i) max(res[i,])))
}) %>% bind_rows()


group_by(lm_res, fold) %>%
  summarise(Accuracy = ACC(dataset, Prediction),
            Kappa = KAPPA(dataset, Prediction),
            AU1U = multiclass.AU1U(Probability, dataset))

group_by(rf_res, fold) %>%
  summarise(Accuracy = ACC(dataset, Prediction),
            Kappa = KAPPA(dataset, Prediction),
            AU1U = multiclass.AU1U(Probability, dataset))

table(lm_res[["dataset"]], lm_res[["Prediction"]])
table(rf_res[["dataset"]], rf_res[["Prediction"]])


