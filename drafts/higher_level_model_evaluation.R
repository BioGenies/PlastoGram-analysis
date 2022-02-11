library(nnet)
library(dplyr)
library(targets)
library(pbapply)
library(irr)

tar_load(data_df)
#PATH <- "lower_level_models_results_debug/"
PATH <- "lower_level_models_results/"
list.files(PATH)
Jackknife_results <- read.csv("/home/jakubkala/Downloads/Jackknife_results.csv")

train_model <- function (train_df, test_df)
{
  hl_model <- invisible(multinom(dataset ~ ., data = train_df,
                                 model = TRUE))
  y_pred <- predict(hl_model, test_df)
}

scaled_train_model <- function (train_df, test_df) {
  #browser()
  scaled_train_df <- scale(select(train_df, -dataset))
  scaled_test_df <- scale(test_df, center=attr(scaled_train_df, "scaled:center"),
                          scale=attr(scaled_train_df, "scaled:scale"))
  scaled_train_df <- data.frame(scaled_train_df)
  scaled_train_df$dataset <- train_df$dataset

  hl_model <- invisible(multinom(dataset ~ ., data = scaled_train_df, model = TRUE))

  y_pred <- predict(hl_model, scaled_test_df)

}


jacknife_higher_level_model <- function(path, data_df, higher_level_model) {

  pred_cols <- c("Nuclear_model", "Membrane_model",
                 "Nuclear_OM_model", "Nuclear_IM_model",
                 "Nuclear_TM_model", "Plastid_membrane_model",
                 "Tat_model", "Sec_model",
                 "Nuclear_OM_stroma_model")

  out <- pblapply(list.files(path), function(filepath){

    seqname <- unlist(strsplit(filepath, "[.]"))[1]
    df <- read.csv(file.path(PATH, filepath))

    train_df <-  left_join(df,
                           data_df[, c("seq_name", "dataset")],
                           by = "seq_name") %>%
      select(-c(seq_name, fold))

    test_df <- Jackknife_results %>% filter(fold == seqname) %>% select(pred_cols)
    y_test <- Jackknife_results %>% filter(fold == seqname) %>% select("dataset")


    y_pred <- higher_level_model(train_df, test_df)



    list(seq_name = seqname, y_test = y_test[1, 1], y_pred = as.character(y_pred))

  })

  out <- data.frame(do.call(rbind, out))
  out <- data.frame(lapply(out, as.character))
  data.frame(lapply(out, function(x) unlist(x)))
}

accuracy <- function(results) {
  sum(results$y_test == results$y_pred) / nrow(results)
}

jackknife_results <- jacknife_higher_level_model(PATH, data_df, train_model)



scaled_jackknife_results <- jacknife_higher_level_model(PATH, data_df, scaled_train_model)

accuracy(jackknife_results)
kappa2(jackknife_results[, c("y_test", "y_pred")])
accuracy(scaled_jackknife_results)
kappa2(scaled_jackknife_results[, c("y_test", "y_pred")])



