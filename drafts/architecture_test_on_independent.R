library(targets)
library(dplyr)
library(biogram)
library(readxl)
library(ranger)
library(seqR)
library(cvTools)
library(measures)
library(rhmmer)
library(purrr)
library(combinat)
library(FCBF)
library(UBL)
library(nnet)
library(rel)
library(tidyr)


source("./functions/create_datasets.R")
source("./functions/filter_sequences.R")
source("./functions/ngram_model_functions.R")
source("./functions/create_target_df.R")
source("./functions/profileHMM_functions.R")
source("./functions/evaluate_full_model.R")
source("./functions/ensemble_model_functions.R")
source("./functions/generate_architectures.R")


run_graphpart_new <- function(graphpart_input, out_dir, ith_seqlist, n_threads = 24) {
  write_fasta(graphpart_input, paste0(out_dir, "Datasets_for_graph-part_", ith_seqlist, ".fa"))
  system(paste0("graphpart needle --fasta-file ", out_dir, "Datasets_for_graph-part_", ith_seqlist, ".fa --threshold 0.5 --out-file ", 
                out_dir, "Graph-part_results_", ith_seqlist, ".csv --labels-name label --partitions 4 --threads 24 --no-moving"))
  label_df <- data.frame(dataset = unique(sapply(names(graphpart_input), function(i) gsub("label=", "", strsplit(i, "|", fixed = TRUE)[[1]][2]))),
                         `label.val` = 0:8)
  res <- read.csv(paste0(out_dir, "Graph-part_results_", ith_seqlist, ".csv")) 
  left_join(res, label_df)
}

get_all_models_predictions_new <- function(ngram_matrix, sequences, data_df, model_df, remove_hmm_files = FALSE) {
  
  dat <- ngram_matrix[data_df[["fold"]] == "train", ]
  train_df <- filter(data_df, fold == "train")
  train_seqs <- sequences[which(names(sequences) %in% train_df[["seq_name"]])]
  
  ngram_models <- train_ngram_models(model_df, dat, train_df, filtering_colname = NULL, filtering_term = NULL)
  hmm_models <- train_profile_HMM_models(model_df, train_seqs, train_df, filtering_colname = NULL, filtering_term = NULL)
  
  predict_with_all_models(model_df, ngram_matrix, data_df, sequences, ngram_models, hmm_models, ith_fold = "independent", remove_hmm_files = remove_hmm_files)
  
}

create_target_df_old <- function(annotations_file, sequences) {
  read_xlsx(annotations_file) %>% 
    select(Entry, Final_dataset) %>% 
    unique() %>% 
    setNames(c("seq_name", "dataset")) %>% 
    filter(seq_name %in% names(sequences),
           dataset %in% c("N_OM", "N_IM", "P_IM", "N_S", "P_S", "N_TM", "P_TM", "N_TL_SEC", "N_TL_TAT")) %>% 
    mutate(Nuclear_target = ifelse(grepl("N_", dataset), TRUE, FALSE),
           Sec_target = ifelse(dataset == "N_TL_SEC", TRUE, FALSE),
           Tat_target = ifelse(dataset == "N_TL_TAT", TRUE, FALSE),
           Membrane_mc_target = case_when(dataset == "N_OM" ~ "OM",
                                          dataset %in% c("N_IM", "P_IM") ~ "IM",
                                          dataset %in% c("N_TM", "P_TM") ~ "TM",
                                          !(dataset %in% c("N_OM", "N_IM", "N_TM", "P_IM", "P_TM")) ~ "other"),
           Stroma_target = ifelse(grepl("_S$", dataset), TRUE, FALSE),
           Membrane_target = ifelse(dataset %in% c("N_OM", "N_IM", "P_IM", "N_TM", "P_TM"), TRUE, FALSE),
           OM_target = ifelse(dataset == "N_OM", TRUE, FALSE),
           IM_target = ifelse(dataset %in% c("N_IM", "P_IM"), TRUE, FALSE),
           TM_target = ifelse(dataset %in% c("N_TM", "P_TM"), TRUE, FALSE))
}

generate_results_for_architectures_new <- function(architecture_file_list, all_models_results, outdir, data_df, higher_level_models = c("RF", "GLM")) {
  lapply(higher_level_models, function(ith_hl_model) {
    lapply(architecture_file_list, function(ith_file) {
      model <- gsub(".csv", "", last(strsplit(ith_file, "/")[[1]]))
      res <- left_join(data_df[, c("seq_name", "dataset")], 
                       filter_results_for_single_architecture(ith_file, all_models_results),
                       by = "seq_name")
      
      train_dat <- select(res[which(data_df[["fold"]] == "train"), ], -c(seq_name, fold))
      test_dat <- res[which(data_df[["fold"]] == "test"), ]
      full_results <- if(ith_hl_model == "GLM") {
        lm_model <- multinom(dataset ~ ., train_dat, model = TRUE)
        preds <- test_dat %>%
          select("dataset") %>%
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
        cbind(select(test_dat, "dataset"),
              Prediction = c(colnames(preds)[max.col(preds[, c(colnames(preds))])]),
              preds,
              Probability = sapply(1:nrow(preds), function(i) max(preds[i,])))
      }
      write.csv(full_results, paste0(outdir, model, "_", ith_hl_model, "_results.csv"), row.names = FALSE)
    })
  })
}

evaluate_all_architectures_new <- function(res_files, outfile, data_df) {
  performance_results <- lapply(res_files, function(ith_file) {
    res <- read.csv(ith_file)
    x <- table(data_df[["dataset"]])
    data.frame(
      model = gsub("_results.csv", "", last(strsplit(ith_file, "/")[[1]])),
      AU1U = multiclass.AU1U(res[, 3:(ncol(res)-1)], res[["dataset"]]),
      kappa = KAPPA(res[["dataset"]], res[["Prediction"]]),
      weighted_kappa = ckap(res[, c("dataset", "Prediction")], 
                            weight = as.matrix(bind_rows(lapply(x, function(i) ifelse(i == x, as.vector(sum(x)/(i+x)), 1)))))[["est"]],
      N_IM_sensitivity = TPR(ifelse(res[["dataset"]] == "N_IM", TRUE, FALSE),
                             ifelse(res[["Prediction"]] == "N_IM", TRUE, FALSE), TRUE),
      N_OM_sensitivity = TPR(ifelse(res[["dataset"]] == "N_OM", TRUE, FALSE),
                             ifelse(res[["Prediction"]] == "N_OM", TRUE, FALSE), TRUE),
      N_TM_sensitivity = TPR(ifelse(res[["dataset"]] == "N_TM", TRUE, FALSE),
                             ifelse(res[["Prediction"]] == "N_TM", TRUE, FALSE), TRUE),
      N_S_sensitivity = TPR(ifelse(res[["dataset"]] == "N_S", TRUE, FALSE),
                            ifelse(res[["Prediction"]] == "N_S", TRUE, FALSE), TRUE),
      N_TL_SEC_sensitivity = TPR(ifelse(res[["dataset"]] == "N_TL_SEC", TRUE, FALSE),
                                 ifelse(res[["Prediction"]] == "N_TL_SEC", TRUE, FALSE), TRUE),
      N_TL_TAT_sensitivity = TPR(ifelse(res[["dataset"]] == "N_TL_TAT", TRUE, FALSE),
                                 ifelse(res[["Prediction"]] == "N_TL_TAT", TRUE, FALSE), TRUE),
      P_IM_sensitivity = TPR(ifelse(res[["dataset"]] == "P_IM", TRUE, FALSE),
                             ifelse(res[["Prediction"]] == "P_IM", TRUE, FALSE), TRUE),
      P_TM_sensitivity = TPR(ifelse(res[["dataset"]] == "P_TM", TRUE, FALSE),
                             ifelse(res[["Prediction"]] == "P_TM", TRUE, FALSE), TRUE),
      P_S_sensitivity = TPR(ifelse(res[["dataset"]] == "P_S", TRUE, FALSE),
                            ifelse(res[["Prediction"]] == "P_S", TRUE, FALSE), TRUE)
    )
  }) %>% bind_rows()
  write.csv(performance_results, outfile, row.names = FALSE)
  performance_results
}



withr::with_dir("./", 
                tar_load(c(N_OM_seqs, N_IM_seqs, P_IM_seqs, N_S_seqs, P_S_seqs,
                           N_TM_seqs, P_TM_seqs, N_TL_SEC_seqs, N_TL_TAT_seqs,
                           sequences, annotations_file, model_dat)))

architecture_files <- list.files("~/Dropbox/Projekty/BioNgramProjects/PlastoGram/Model_architectures/", 
                                 full.names = TRUE)

seqlist1 <- list("N_OM" = N_OM_seqs, "N_IM" = N_IM_seqs, "P_IM" = P_IM_seqs,
                 "N_S" = N_S_seqs, "P_S" = P_S_seqs, "N_TM" = N_TM_seqs, "P_TM" = P_TM_seqs, 
                 "N_TL_SEC" = N_TL_SEC_seqs, "N_TL_TAT" = N_TL_TAT_seqs)
seqlist2 <- list("N_IM" = N_IM_seqs, "P_IM" = P_IM_seqs, "N_TM" = N_TM_seqs,
                 "N_S" = N_S_seqs, "N_TL_SEC" = N_TL_SEC_seqs, "P_S" = P_S_seqs,  
                 "P_TM" = P_TM_seqs, "N_OM" = N_OM_seqs, "N_TL_TAT" = N_TL_TAT_seqs)
seqlist3 <- list("N_TL_SEC" = N_TL_SEC_seqs, "P_S" = P_S_seqs, "P_TM" = P_TM_seqs,
                 "N_IM" = N_IM_seqs, "P_IM" = P_IM_seqs, "N_TM" = N_TM_seqs,
                 "N_S" = N_S_seqs, "N_TL_TAT" = N_TL_TAT_seqs, "N_OM" = N_OM_seqs)
seqlist4 <- list("N_TL_TAT" = N_TL_TAT_seqs, "P_IM" = P_IM_seqs, "P_S" = P_S_seqs, 
                 "N_TL_SEC" = N_TL_SEC_seqs,  "P_TM" = P_TM_seqs, "N_S" = N_S_seqs, 
                 "N_OM" = N_OM_seqs, "N_IM" = N_IM_seqs,  "N_TM" = N_TM_seqs)
seqlist5 <- list("P_IM" = P_IM_seqs, "N_OM" = N_OM_seqs, "N_TL_TAT" = N_TL_TAT_seqs,
                 "P_S" = P_S_seqs, "P_TM" = P_TM_seqs, "N_S" = N_S_seqs, 
                 "N_TL_SEC" = N_TL_SEC_seqs, "N_IM" = N_IM_seqs,  "N_TM" = N_TM_seqs)

seqlist6 <- list("N_OM" = filter_nonstandard_aa(sequences[["N_OM"]]), "N_IM" = filter_nonstandard_aa(sequences[["N_IM"]]), "P_IM" = filter_nonstandard_aa(sequences[["P_IM"]]),
                 "N_S" = filter_nonstandard_aa(sequences[["N_S"]]), "P_S" = filter_nonstandard_aa(sequences[["P_S"]]), "N_TM" = filter_nonstandard_aa(sequences[["N_TM"]]), 
                 "P_TM" = filter_nonstandard_aa(sequences[["P_TM"]]), "N_TL_SEC" = filter_nonstandard_aa(sequences[["N_TL_SEC"]]), "N_TL_TAT" = filter_nonstandard_aa(sequences[["N_TL_TAT"]]))
seqlist7 <- list("N_IM" = filter_nonstandard_aa(sequences[["N_IM"]]), "P_IM" = filter_nonstandard_aa(sequences[["P_IM"]]), "N_TM" = filter_nonstandard_aa(sequences[["N_TM"]]),
                 "N_S" = filter_nonstandard_aa(sequences[["N_S"]]), "N_TL_SEC" = filter_nonstandard_aa(sequences[["N_TL_SEC"]]), "P_S" = filter_nonstandard_aa(sequences[["P_S"]]),  
                 "P_TM" = filter_nonstandard_aa(sequences[["P_TM"]]), "N_OM" = filter_nonstandard_aa(sequences[["N_OM"]]), "N_TL_TAT" = filter_nonstandard_aa(sequences[["N_TL_TAT"]]))
seqlist8 <- list("N_TL_SEC" = filter_nonstandard_aa(sequences[["N_TL_SEC"]]), "P_S" = filter_nonstandard_aa(sequences[["P_S"]]), "P_TM" = filter_nonstandard_aa(sequences[["P_TM"]]),
                 "N_IM" = filter_nonstandard_aa(sequences[["N_IM"]]), "P_IM" = filter_nonstandard_aa(sequences[["P_IM"]]), "N_TM" = filter_nonstandard_aa(sequences[["N_TM"]]),
                 "N_S" = filter_nonstandard_aa(sequences[["N_S"]]), "N_TL_TAT" = filter_nonstandard_aa(sequences[["N_TL_TAT"]]), "N_OM" = filter_nonstandard_aa(sequences[["N_OM"]]))
seqlist9 <- list("N_TL_TAT" = filter_nonstandard_aa(sequences[["N_TL_TAT"]]), "P_IM" = filter_nonstandard_aa(sequences[["P_IM"]]), "P_S" = filter_nonstandard_aa(sequences[["P_S"]]), 
                 "N_TL_SEC" = filter_nonstandard_aa(sequences[["N_TL_SEC"]]),  "P_TM" = filter_nonstandard_aa(sequences[["P_TM"]]), "N_S" = filter_nonstandard_aa(sequences[["N_S"]]), 
                 "N_OM" = filter_nonstandard_aa(sequences[["N_OM"]]), "N_IM" = filter_nonstandard_aa(sequences[["N_IM"]]),  "N_TM" = filter_nonstandard_aa(sequences[["N_TM"]]))
seqlist10 <- list("P_IM" = filter_nonstandard_aa(sequences[["P_IM"]]), "N_OM" = filter_nonstandard_aa(sequences[["N_OM"]]), "N_TL_TAT" = filter_nonstandard_aa(sequences[["N_TL_TAT"]]),
                  "P_S" = filter_nonstandard_aa(sequences[["P_S"]]), "P_TM" = filter_nonstandard_aa(sequences[["P_TM"]]), "N_S" = filter_nonstandard_aa(sequences[["N_S"]]), 
                  "N_TL_SEC" = filter_nonstandard_aa(sequences[["N_TL_SEC"]]), "N_IM" = filter_nonstandard_aa(sequences[["N_IM"]]),  "N_TM" = filter_nonstandard_aa(sequences[["N_TM"]]))

seqlists <- list("s1-cd-graphpart" = seqlist1, "s2-cd-graphpart" = seqlist2, "s3-cd-graphpart" = seqlist3, 
                 "s4-cd-graphpart" = seqlist4, "s5-cd-graphpart" = seqlist5, "s1-cd" = seqlist1, 
                 "s2-cd" = seqlist2, "s3-cd" = seqlist3, "s4-cd" = seqlist4, "s5-cd" = seqlist5, 
                 "s1-graphpart" = seqlist6, "s2-graphpart" = seqlist7, "s3-graphpart" = seqlist8, 
                 "s4-graphpart" = seqlist9, "s5-graphpart" = seqlist10)

lapply(names(seqlists), function(ith_seqlist) {
  all_sequences <- unlist(unname(seqlists[[ith_seqlist]]), recursive = FALSE)
  if(grepl("graphpart", ith_seqlist)) {
    graphpart_input <- process_for_graphpart(seqlists[[ith_seqlist]])
    graphpart_res <- run_graphpart_new(graphpart_input, "./drafts/Datasets/", ith_seqlist)
    print(paste0("Done graph-part for ", ith_seqlist, " (dataset ", which(names(seqlists) == ith_seqlist), "/15)"))
    target_df <- create_target_df(annotations_file,
                                  all_sequences,
                                  graphpart_res) %>% 
      mutate(fold = ifelse(fold == 3, "test", "train"))
  } else {
    target_df <- create_target_df_old(annotations_file, all_sequences) %>% 
      create_folds(4) %>% 
      mutate(fold = ifelse(fold == 4, "test", "train")) 
  }
  
  ngram_matrix <- create_ngram_matrix(all_sequences, target_df)
  all_models_predictions <- get_all_models_predictions_new(ngram_matrix, all_sequences, target_df, model_dat, remove_hmm_files = TRUE)
  print(paste0("Done all models predictions for ", ith_seqlist, " (dataset ", which(names(seqlists) == ith_seqlist), "/15)"))
  architecture_results <- generate_results_for_architectures_new(architecture_files, all_models_predictions, 
                                                                 paste0("./drafts/Architecture_results_", ith_seqlist, "/"), 
                                                                 target_df, higher_level_models = c("RF", "GLM"))
  print(paste0("Done generating architecture results for ", ith_seqlist, " (dataset ", which(names(seqlists) == ith_seqlist), "/15)"))
  architecture_results_files <- list.files(paste0("./drafts/Architecture_results_", ith_seqlist), full.names = TRUE)
  
  evaluate_all_architectures_new(architecture_results_files,
                                 paste0("./drafts/Architectures_evaluation_", ith_seqlist, ".csv"),
                                 target_df)
  print(paste0("Done evaluation of all architectures for ", ith_seqlist, " (dataset ", which(names(seqlists) == ith_seqlist), "/15)"))
  
})



### Results

all_res_files <- list.files("~/RProjects/PlastoGram-analysis/drafts/", pattern = "Architectures_evaluation",
                            full.names = TRUE)

all_res_independent <- lapply(all_res_files, function(i) {
  name <- last(strsplit(gsub(".csv", "", last(strsplit(i, "/")[[1]])), "_")[[1]])
  read.csv(i) %>% 
    mutate(type = name,
           reduction = gsub("cd", "cdhit", gsub("s.-", "", name)),
           rep = gsub("s", "", strsplit(name, "-")[[1]][1]),
           cdhit = ifelse(grepl("cd", name), TRUE, FALSE),
           graphpart = ifelse(grepl("graphpart", name), TRUE, FALSE))
}) %>% bind_rows() %>% 
  mutate(arch_variant = sapply(.[["model"]], function(i) gsub("-Sec|-Tat", "", strsplit(i, "_")[[1]][2])),
         hierarchical = grepl(pattern = "Filtering", x = model, ignore.case = FALSE),
         algorithm = sapply(.[["model"]], function(i) last(strsplit(i, "_")[[1]])),
         sec = !grepl(pattern = "-Sec", x = model, ignore.case = FALSE),
         tat = !grepl(pattern = "-Tat", x = model, ignore.case = FALSE),
         smote = !grepl(pattern = "0-1", x = model),
         smote_variant = sapply(.[["model"]], function(i) strsplit(i, "_")[[1]][3])) %>% 
  pivot_longer(cols = !c(type, reduction, rep, cdhit, graphpart, model, arch_variant, hierarchical, algorithm, sec, tat, smote, smote_variant), names_to = "measure") 

library(ggplot2)
library(ggbeeswarm)

all_res_independent %>% 
  filter(measure == "kappa") %>% 
  ggplot(aes(x = rep, y = value)) +
  geom_quasirandom() +
  facet_wrap(~reduction, scales = "free_x") +
  theme_bw(base_size = 6)

p1 <- all_res_independent %>% 
  filter(measure %in% c("kappa", "AU1U")) %>% 
  ggplot(aes(x = rep, y = value)) +
  geom_quasirandom(size = 0.2) +
  facet_grid(measure~reduction, scales = "free") +
  ggtitle("Performance of all architectures on the independent dataset") +
  theme_bw(base_size = 8)

p2 <- all_res_independent %>% 
  filter(measure %in% c("kappa", "AU1U")) %>% 
  ggplot(aes(x = rep, y = value, color = hierarchical)) +
  geom_quasirandom(size = 0.2) +
  facet_grid(measure~reduction, scales = "free") +
  ggtitle("Performance of all architectures on the independent dataset (hierarchy)") +
  theme_bw(base_size = 8)

p3 <- all_res_independent %>% 
  filter(measure %in% c("kappa", "AU1U")) %>% 
  mutate(higher_level_model = ifelse(grepl("GLM", model), "GLM", "RF")) %>% 
  ggplot(aes(x = rep, y = value, color = higher_level_model)) +
  geom_quasirandom(size = 0.2) +
  facet_grid(measure~reduction, scales = "free") +
  ggtitle("Performance of all architectures on the independent dataset (higher-level model)") +
  theme_bw(base_size = 8)

p4 <- all_res_independent %>%
  filter(measure %in% c("kappa", "AU1U")) %>% 
  filter(grepl("GLM", model)) %>% 
  ggplot(aes(x = rep, y = value, color = hierarchical)) +
  geom_quasirandom(size = 0.2) +
  facet_grid(measure~reduction, scales = "free") +
  ggtitle("Performance of architectures with GLM higher-level model on the independent dataset") +
  theme_bw(base_size = 8)

p5 <- all_res_independent %>%
  filter(measure %in% c("kappa", "AU1U")) %>% 
  filter(grepl("RF", model)) %>% 
  ggplot(aes(x = rep, y = value, color = hierarchical)) +
  geom_quasirandom(size = 0.2) +
  facet_grid(measure~reduction, scales = "free") +
  ggtitle("Performance of architectures with RF higher-level model on the independent dataset") +
  theme_bw(base_size = 8)

library(patchwork)
full_plot_of_results <- p1/p2/p3/p4/p5 
ggsave("Independent_dataset_results.png", full_plot_of_results, width = 9, height = 25)

p6 <- all_res_independent %>% 
  filter(measure %in% c("N_OM_sensitivity", "N_IM_sensitivity", "N_TL_SEC_sensitivity")) %>% 
  ggplot(aes(x = rep, y = value, color = hierarchical)) +
  geom_quasirandom(size = 0.2) +
  facet_grid(measure~reduction, scales = "free") +
  ggtitle("Performance for most difficult classes on the independent dataset (hierarchy)") +
  theme_bw(base_size = 8)

p7 <- all_res_independent %>% 
  filter(measure %in% c("N_OM_sensitivity", "N_IM_sensitivity", "N_TL_SEC_sensitivity")) %>% 
  mutate(higher_level_model = ifelse(grepl("GLM", model), "GLM", "RF")) %>% 
  ggplot(aes(x = rep, y = value, color = higher_level_model)) +
  geom_quasirandom(size = 0.2) +
  facet_grid(measure~reduction, scales = "free") +
  ggtitle("Performance for most difficult classes on the independent dataset (higher-level model)") +
  theme_bw(base_size = 8)

difficult_classes <- p6/p7

ggsave("Independent_dataset_results_difficult_classes.png", difficult_classes, width = 9, height = 12)
