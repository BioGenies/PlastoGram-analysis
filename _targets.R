library(dplyr)
library(biogram)
library(targets)
library(seqR)
library(cvTools)
library(ranger)
library(readxl)
library(measures)
library(rhmmer)
library(purrr)
library(combinat)
library(FCBF)
library(UBL)
library(nnet)
library(rel)

data_path <- "~/Dropbox/Projekty/BioNgramProjects/PlastoGram/"

source("./functions/create_datasets.R")
source("./functions/filter_sequences.R")
source("./functions/ngram_model_functions.R")
source("./functions/create_target_df.R")
source("./functions/profileHMM_functions.R")
source("./functions/evaluate_full_model.R")
source("./functions/ensemble_model_functions.R")
source("./functions/generate_architectures.R")
source("./functions/baseline_model.R")

set.seed(427244)

list(
  tar_target(
    all_seqs_file,
    paste0(data_path, "All_sequences.fasta"),
    format = "file"
  ),
  tar_target(
    annotations_file,
    paste0(data_path, "Dataset_annotations_references.xlsx"),
    format = "file"
  ),
  tar_target(
    dataset_names,
    c("N_OM", "N_IM", "P_IM", "N_S", "P_S", "N_TM", "P_TM", "N_TL_SEC", "N_TL_TAT")
  ),
  tar_target(
    sequences,
    create_datasets(all_seqs_file,
                    annotations_file,
                    paste0(data_path, "Sequences/"),
                    dataset_names),
  ),
  tar_target(
    N_OM_seqs,
    filter_with_cdhit(
      filter_nonstandard_aa(sequences[["N_OM"]]), 
      threshold = 0.9)
  ),
  tar_target(
    N_IM_seqs,
    filter_with_cdhit(
      filter_nonstandard_aa(sequences[["N_IM"]]), 
      threshold = 0.9)
  ),
  tar_target(
    P_IM_seqs,
    filter_with_cdhit(
      filter_nonstandard_aa(sequences[["P_IM"]]), 
      threshold = 0.9)
  ),
  tar_target(
    N_S_seqs, 
    filter_with_cdhit(
      filter_nonstandard_aa(sequences[["N_S"]]), 
      threshold = 0.9)
  ),
  tar_target(
    P_S_seqs, 
    filter_with_cdhit(
      filter_nonstandard_aa(sequences[["P_S"]]),
      threshold = 0.9)
  ),
  tar_target(
    N_TM_seqs, 
    filter_with_cdhit(
      filter_nonstandard_aa(sequences[["N_TM"]]),
      threshold = 0.9)
  ),
  tar_target(
    P_TM_seqs, 
    filter_with_cdhit(
      filter_nonstandard_aa(sequences[["P_TM"]]),
      threshold = 0.9)
  ),
  tar_target(
    N_TL_SEC_seqs, 
    filter_with_cdhit(
      filter_nonstandard_aa(sequences[["N_TL_SEC"]]),
      threshold = 0.9)
  ),
  tar_target(
    N_TL_TAT_seqs, 
    filter_with_cdhit(
      filter_nonstandard_aa(sequences[["N_TL_TAT"]]),
      threshold = 0.9)
  ),
  tar_target(
    sequence_list,
    list("N_OM" = N_OM_seqs, "N_IM" = N_IM_seqs, "P_IM" = P_IM_seqs,
         "N_S" = N_S_seqs, "P_S" = P_S_seqs, "N_TM" = N_TM_seqs, "P_TM" = P_TM_seqs, 
         "N_TL_SEC" = N_TL_SEC_seqs, "N_TL_TAT" = N_TL_TAT_seqs)
  ),
  tar_target(
    all_sequences,
    unlist(unname(sequence_list), recursive = FALSE)
  ),
  tar_target(
    graphpart_input,
    process_for_graphpart(sequence_list)
  ),
  tar_target(
    graphpart_res,
    run_graphpart(graphpart_input, 
                  paste0(data_path, "Sequences/"))
  ),
  tar_target(
    target_df, 
    create_target_df(annotations_file,
                     all_sequences,
                     graphpart_res)),
  tar_target(
    data_df,
    filter(target_df, cluster != 1)
  ),
  tar_target(
    data_df_independent,
    filter(target_df, cluster == 1)
  ),
  tar_target(
    ngram_matrix, 
    create_ngram_matrix(all_sequences, data_df)
  ),
  tar_target(
    ngram_matrix_independent, 
    create_ngram_matrix(all_sequences, data_df_independent)
  ),
  tar_target(
    sequences_cv,
    all_sequences[which(names(all_sequences) %in% data_df[["seq_name"]])]
  ),
  tar_target(
    sequences_independent,
    all_sequences[which(names(all_sequences) %in% data_df_independent[["seq_name"]])]
  ),
  tar_target(
    model_variants,
    get_model_variants()
  ),
  tar_target(
    model_dat_file,
    "./data/PlastoGram_models_info.csv",
    format = "file"
  ),
  tar_target(
    model_dat,
    read.csv(model_dat_file)
  ),
  tar_target(
    test_filtering_options_file,
    "./data/PlastoGram_models_results_filtering.csv",
    format = "file"
  ),
  tar_target(
    filtering_df,
    read.csv(test_filtering_options_file)
  ),
  tar_target(
    target_dfs_cv,
    create_cv_folds(data_df, 5, c(427244, 58713042, 6513, 901374, 388123648, 43671, 71234, 209147, 1847820, 248114))
  ),
  tar_target(
    all_models_predictions,
    get_all_models_predictions_cv(ngram_matrix, sequences_cv, target_dfs_cv, model_dat, data_path, remove_hmm_files = TRUE)
  ),
  tar_target(
    architectures_performance,
    generate_and_test_architectures(model_variants = model_variants,
                                    smote_models = c("OM_Stroma_model", "Nuclear_membrane_model",
                                                     "N_OM_model", "N_IM_model", "N_TM_model"),
                                    sequence_models = c("Sec_model", "Tat_model"),
                                    model_dat = model_dat,
                                    filtering_df = filtering_df,
                                    architectures_output_dir = paste0(data_path, "Model_architectures/"),
                                    all_models_predictions = all_models_predictions, 
                                    architecture_res_output_dir = paste0(data_path, "Model_architectures_results/"), 
                                    data_df_final = data_df,
                                    performance_outfile = paste0(data_path, "Architectures_performance.csv"))
  ),
  tar_target(
    mean_architecture_performance,
    get_mean_performance_of_architectures(architectures_performance,
                                          paste0(data_path, "Architectures_mean_performance.csv"))
  # ),
  # tar_target(
  #   PlastoGram_best_architecture_glm,
  #   read.csv(paste0(data_path, "Model_architectures/Architecture_v8_0-1_No_filtering.csv"))
  # ),
  # tar_target(
  #   PlastoGram_best_architecture_rf,
  #   read.csv(paste0(data_path, "Model_architectures/Architecture_v8_1-2_No_filtering.csv"))
  # ),
  # tar_target(
  #   PlastoGram_ngram_models,
  #   train_ngram_models(PlastoGram_final_architecture, ngram_matrix, data_df_final, filtering_colname = NULL, filtering_term = NULL)
  # ),
  # tar_target(
  #   PlastoGram_hmm_Sec,
  #   paste0(train_profileHMM(N_TL_SEC_seqs, "PlastoGram_Sec_model", remove_files = TRUE), ".hmm"),
  #   format = "file"
  # ),
  # tar_target(
  #   PlastoGram_hmm_Tat,
  #   paste0(train_profileHMM(N_TL_TAT_seqs, "PlastoGram_Tat_model", remove_files = TRUE), ".hmm"),
  #   format = "file"
  # ),
  # tar_target(
  #   PlastoGram_predictions,
  #   predict_with_all_models(PlastoGram_final_architecture, ngram_matrix, data_df_final, c(N_seqs, P_seqs), PlastoGram_ngram_models, 
  #                           list("Sec_model" = gsub(".hmm", "", PlastoGram_hmm_Sec), "Tat_model" = gsub(".hmm", "", PlastoGram_hmm_Tat)), 
  #                           ith_fold = NULL, remove_hmm_files = FALSE)
  # ),
  # tar_target(
  #   PlastoGram_multinom_model,
  #   train_multinom(PlastoGram_predictions, data_df_final)
  ),
  tar_target(
    baseline_model_cv_res,
    do_baseline_cv(ngram_matrix, target_dfs_cv)
  # ),
  # tar_target(
  #   jackknife_results_glm,
  #   do_jackknife(ngram_matrix, c(N_seqs, P_seqs), data_df_final, PlastoGram_best_architecture_glm, data_path, higher_level_model = "GLM")
  # ),
  # tar_target(
  #   jackknife_results_rf,
  #   do_jackknife(ngram_matrix, c(N_seqs, P_seqs), data_df_final, PlastoGram_best_architecture_rf, data_path, higher_level_model = "RF")
  )
)

