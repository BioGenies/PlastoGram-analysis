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
library(ggplot2)
library(tidyr)
library(patchwork)

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
source("./functions/publication_results.R")

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
    model_variants,
    get_updated_model_variants()
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
  ### Holdout approach with OM and IM 
  tar_target(
    holdouts,
    generate_holdout(sequence_list)
  ),
  tar_target(
    traintest,
    unlist(unname(holdouts[["traintest"]]), recursive = FALSE)
  ),
  tar_target(
    benchmark,
    unlist(unname(holdouts[["benchmark"]]), recursive = FALSE)
  ),
  tar_target(
    holdout_target_df,
    create_target_df_holdout(annotations_file,
                             all_sequences)
  ),
  tar_target(
    traintest_data_df,
    filter(holdout_target_df, seq_name %in% names(traintest))
  ),
  tar_target(
    benchmark_data_df,
    filter(holdout_target_df, seq_name %in% names(benchmark))
  ),
  tar_target(
    traintest_ngram_matrix, 
    create_ngram_matrix(traintest, traintest_data_df)
  ),
  tar_target(
    benchmark_ngram_matrix, 
    create_ngram_matrix(benchmark, benchmark_data_df)
  ),
  tar_target(
    holdout_target_dfs_cv,
    create_cv_folds(traintest_data_df, 5, c(427244, 58713042, 6513, 901374, 388123648)),
  ),
  tar_target(
    holdout_all_models_predictions,
    get_all_models_predictions_cv(traintest_ngram_matrix, traintest, holdout_target_dfs_cv, 
                                  model_dat, data_path, type = "holdout", remove_hmm_files = TRUE)
  ),
  tar_target(
    holdout_architectures_performance,
    generate_and_test_architectures(model_variants = model_variants,
                                    smote_models = c("OM_Stroma_model", "Nuclear_membrane_model",
                                                     "Nuclear_membrane_mc_model",
                                                     "N_OM_model", "N_IM_model", "N_TM_model"),
                                    sequence_models = c("Sec_model", "Tat_model"),
                                    model_dat = model_dat,
                                    filtering_df = filtering_df,
                                    architectures_output_dir = paste0(data_path, "Model_architectures/"),
                                    all_models_predictions = holdout_all_models_predictions, 
                                    architecture_res_output_dir = paste0(data_path, "Model_architectures_holdout_results/"), 
                                    data_df_final = traintest_data_df,
                                    performance_outfile = paste0(data_path, "Architectures_performance_holdout.csv"),
                                    type = "holdout")
  ),
  tar_target(
    holdout_mean_architecture_performance,
    get_mean_performance_of_architectures(holdout_architectures_performance,
                                          paste0(data_path, "Architectures_mean_performance_holdout.csv"))
  ),
  tar_target(
    holdout_ranked_architecture_performance,
    rank_architectures(holdout_mean_architecture_performance)
  ),
  ### OM and IM together as envelope
  tar_target(
    envelope_target_df,
    create_envelope_target_df(annotations_file,
                              all_sequences)
  ),
  tar_target(
    envelope_traintest_data_df,
    filter(envelope_target_df, seq_name %in% names(traintest))
  ),
  tar_target(
    envelope_benchmark_data_df,
    filter(envelope_target_df, seq_name %in% names(benchmark))
  ),
  tar_target(
    envelope_traintest_ngram_matrix, 
    create_ngram_matrix(traintest, envelope_traintest_data_df)
  ),
  tar_target(
    envelope_benchmark_ngram_matrix, 
    create_ngram_matrix(benchmark, envelope_benchmark_data_df)
  ),
  tar_target(
    envelope_target_dfs_cv,
    create_cv_folds(envelope_traintest_data_df, 5, c(427244, 58713042, 6513, 901374, 388123648)),
  ),
  tar_target(
    envelope_model_variants,
    get_envelope_model_variants()
  ),
  tar_target(
    envelope_model_dat_file,
    "./data/PlastoGram_envelope_models_info.csv",
    format = "file"
  ),
  tar_target(
    envelope_model_dat,
    read.csv(envelope_model_dat_file)
  ),
  tar_target(
    envelope_test_filtering_options_file,
    "./data/PlastoGram_envelope_models_results_filtering.csv",
    format = "file"
  ),
  tar_target(
    envelope_filtering_df,
    read.csv(envelope_test_filtering_options_file)
  ),
  tar_target(
    envelope_all_models_predictions,
    get_all_models_predictions_cv(envelope_traintest_ngram_matrix, traintest, envelope_target_dfs_cv, 
                                  envelope_model_dat, data_path, type = "envelope", remove_hmm_files = TRUE)
  ),
  tar_target(
    envelope_architectures_performance,
    generate_and_test_architectures(model_variants = envelope_model_variants,
                                    smote_models = c("N_E_all_model", "N_TM_all_model", "TM_all_model"),
                                    sequence_models = c("Sec_model", "Tat_model"),
                                    model_dat = envelope_model_dat,
                                    filtering_df = envelope_filtering_df,
                                    architectures_output_dir = paste0(data_path, "Model_architectures_envelope/"),
                                    all_models_predictions = envelope_all_models_predictions, 
                                    architecture_res_output_dir = paste0(data_path, "Model_architectures_envelope_results/"), 
                                    data_df_final = envelope_traintest_data_df,
                                    performance_outfile = paste0(data_path, "Architectures_envelope_performance.csv"),
                                    type = "envelope")
  ),
  tar_target(
    envelope_mean_architecture_performance,
    get_mean_performance_of_envelope_architectures(envelope_architectures_performance,
                                                   paste0(data_path, "Architectures_envelope_mean_performance.csv"))
  ),
  tar_target(
    envelope_ranked_architecture_performance,
    rank_architectures(envelope_mean_architecture_performance)
  ),
  # Results and final models
  tar_target(
    PlastoGram_best_architecture_name,
    arrange(envelope_mean_architecture_performance, desc(mean_kappa))[["model"]][1]
  ),
  tar_target(
    PlastoGram_best_architecture,
    read.csv(paste0(data_path, "Model_architectures_envelope/", gsub("_GLM|_RF", ".csv", PlastoGram_best_architecture_name)))
  ),
  tar_target(
    PlastoGram_ngram_models,
    train_ngram_models(PlastoGram_best_architecture, envelope_traintest_ngram_matrix, envelope_traintest_data_df, filtering_colname = NULL, filtering_term = NULL)
  ),
  tar_target(
    PlastoGram_hmm_Sec,
    paste0(train_profileHMM(holdouts[["traintest"]][["N_TL_SEC"]], "PlastoGram_Sec_model", remove_files = TRUE), ".hmm"),
    format = "file"
  ),
  tar_target(
    PlastoGram_hmm_Tat,
    paste0(train_profileHMM(holdouts[["traintest"]][["N_TL_TAT"]], "PlastoGram_Tat_model", remove_files = TRUE), ".hmm"),
    format = "file"
  ),
  tar_target(
    PlastoGram_predictions,
    predict_with_all_models(PlastoGram_best_architecture, envelope_traintest_ngram_matrix, envelope_traintest_data_df, traintest, PlastoGram_ngram_models,
                            list("Sec_model" = gsub(".hmm", "", PlastoGram_hmm_Sec), "Tat_model" = gsub(".hmm", "", PlastoGram_hmm_Tat)),
                            ith_fold = NULL, remove_hmm_files = FALSE)
  ),
  tar_target(
    PlastoGram_higher_level_model,
    train_higher_level_model(PlastoGram_predictions, envelope_traintest_data_df, PlastoGram_best_architecture_name, PlastoGram_best_architecture)
  ),
  tar_target(
    PlastoGram_informative_ngrams,
    get_all_imp_ngrams(PlastoGram_ngram_models)
  ),
  tar_target(
    PlastoGram_evaluation,
    predict_with_PlastoGram(PlastoGram_ngram_models, 
                            list("Sec_model" = gsub(".hmm", "", PlastoGram_hmm_Sec), "Tat_model" = gsub(".hmm", "", PlastoGram_hmm_Tat)), 
                            PlastoGram_higher_level_model, envelope_benchmark_ngram_matrix, benchmark, envelope_benchmark_data_df,
                            PlastoGram_informative_ngrams, PlastoGram_best_architecture, PlastoGram_best_architecture_name)
  ),
  tar_target(
    PlastoGram_OM_IM_model,
    train_om_im_model(envelope_traintest_ngram_matrix, envelope_traintest_data_df)
  ),
  tar_target(
    PlastoGram_evaluation_OM_IM,
    evaluate_om_im_model(PlastoGram_OM_IM_model, envelope_benchmark_ngram_matrix, envelope_benchmark_data_df, PlastoGram_evaluation)
  ),
  tar_target(
    baseline_model_cv_res,
    do_baseline_cv(envelope_traintest_ngram_matrix, envelope_target_dfs_cv)
  ),
  tar_target(
    baseline_mean_performance,
    get_mean_performance_of_baseline(baseline_model_cv_res, paste0(data_path, "Baseline_envelope_mean_performance.csv"))
  ),
  tar_target(
    baseline_model,
    train_baseline_model(envelope_traintest_ngram_matrix, envelope_traintest_data_df)
  ),
  tar_target(
    baseline_model_evaluation,
    evaluate_baseline_model(baseline_model, envelope_benchmark_ngram_matrix, envelope_benchmark_data_df)
  ),
  tar_target(
    architecture_plot_data,
    get_architecture_plot_data(envelope_mean_architecture_performance)
  ),
  tar_target(
    performance_distr_plot,
    get_performance_distr_plot(architecture_plot_data, baseline_mean_performance, paste0(data_path, "Publication_results/"))
  ),
  tar_target(
    performance_parameters_table,
    get_parameters_of_performance_distribution_table(envelope_mean_architecture_performance, paste0(data_path, "Publication_results/"))
  ),
  tar_target(
    best_model_cv_plot,
    get_best_model_cv_plot(architecture_plot_data, paste0(data_path, "Publication_results/"))
  ),
  tar_target(
    best_model_cv_table,
    get_best_model_cv_table(architecture_plot_data, paste0(data_path, "Publication_results/"))
  ),
  tar_target(
    independent_dataset_performance_table,
    get_independent_dataset_performance_table(PlastoGram_evaluation, paste0(data_path, "Publication_results/"))
  ),
  tar_target(
    om_im_model_cv_performance_table,
    get_om_im_model_cv_res_table(traintest_ngram_matrix, holdout_target_dfs_cv, paste0(data_path, "Publication_results/")) 
  )
)

