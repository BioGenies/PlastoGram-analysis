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
library(ggpubr)

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
source("./functions/pairwise_identity.R")

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
  
  # CD-HIT + graphpart approach with OM and IM as envelope class
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
    target_df_graphpart, 
    create_target_df(annotations_file,
                     all_sequences,
                     graphpart_res)),
  tar_target(
    data_df_graphpart,
    filter(target_df_graphpart, cluster != 1)
  ),
  tar_target(
    data_df_independent_graphpart,
    filter(target_df_graphpart, cluster == 1)
  ),
  tar_target(
    ngram_matrix_traintest_graphpart, 
    create_ngram_matrix(all_sequences, data_df_graphpart)
  ),
  tar_target(
    ngram_matrix_independent_graphpart, 
    create_ngram_matrix(all_sequences, data_df_independent_graphpart)
  ),
  tar_target(
    sequences_cv_graphpart,
    all_sequences[which(names(all_sequences) %in% data_df_graphpart[["seq_name"]])]
  ),
  tar_target(
    sequences_independent_graphpart,
    all_sequences[which(names(all_sequences) %in% data_df_independent_graphpart[["seq_name"]])]
  ),
  tar_target(
    filtering_df_graphpart,
    read.csv(test_filtering_options_file)
  ),
  tar_target(
    target_dfs_cv_graphpart,
    create_cv_folds(data_df_graphpart, 5, c(427244, 58713042, 6513, 901374, 388123648)),#, 43671, 71234, 209147, 1847820, 248114))
  ),
  tar_target(
    envelope_all_models_predictions_graphpart,
    get_all_models_predictions_cv(ngram_matrix_traintest_graphpart, sequences_cv_graphpart, target_dfs_cv_graphpart, 
                                  envelope_model_dat, data_path, type = "envelope", remove_hmm_files = TRUE)
  ),
  tar_target(
    architectures_performance_graphpart,
    generate_and_test_architectures(model_variants = envelope_model_variants,
                                    smote_models = c("N_E_all_model", "N_TM_all_model", "TM_all_model"),
                                    sequence_models = c("Sec_model", "Tat_model"),
                                    model_dat = envelope_model_dat,
                                    filtering_df = envelope_filtering_df,
                                    architectures_output_dir = paste0(data_path, "Model_architectures_envelope_graphpart/"),
                                    all_models_predictions = envelope_all_models_predictions_graphpart, 
                                    architecture_res_output_dir = paste0(data_path, "Model_architectures_envelope_results_graphpart/"), 
                                    data_df_final = data_df_graphpart,
                                    performance_outfile = paste0(data_path, "Architectures_envelope_performance_graphpart.csv"),
                                    type = "envelope")
  ),
  tar_target(
    mean_architecture_performance_graphpart,
    get_mean_performance_of_envelope_architectures(architectures_performance_graphpart,
                                                   paste0(data_path, "Architectures_mean_performance_graphpart.csv"))
  ),
  tar_target(
    PlastoGram_best_architecture_name_graphpart,
    arrange(mean_architecture_performance_graphpart, desc(mean_kappa))[["model"]][1]
  ),
  tar_target(
    PlastoGram_best_architecture_graphpart,
    change_model_names(read.csv(paste0(data_path, "Model_architectures_envelope_graphpart/", gsub("_GLM|_RF", ".csv", PlastoGram_best_architecture_name_graphpart))))
  ),
  tar_target(
    PlastoGram_ngram_models_graphpart,
    train_ngram_models(PlastoGram_best_architecture_graphpart, ngram_matrix_traintest_graphpart, data_df_graphpart, filtering_colname = NULL, filtering_term = NULL)
  ),
  tar_target(
    PlastoGram_hmm_Sec_graphpart,
    paste0(train_profileHMM(sequences_cv_graphpart[which(names(sequences_cv_graphpart) %in% filter(data_df_graphpart, dataset == "N_TL_SEC")[["seq_name"]])], 
                            "PlastoGram_graphpart_Sec_model", remove_files = TRUE), ".hmm"),
    format = "file"
  ),
  tar_target(
    PlastoGram_hmm_Tat_graphpart,
    paste0(train_profileHMM(sequences_cv_graphpart[which(names(sequences_cv_graphpart) %in% filter(data_df_graphpart, dataset == "N_TL_TAT")[["seq_name"]])], 
                            "PlastoGram_graphpart_Tat_model", remove_files = TRUE), ".hmm"),
    format = "file"
  ),
  tar_target(
    PlastoGram_predictions_graphpart,
    predict_with_all_models(PlastoGram_best_architecture_graphpart, ngram_matrix_traintest_graphpart, data_df_graphpart, sequences_cv_graphpart, PlastoGram_ngram_models_graphpart,
                            list("Sec_model" = gsub(".hmm", "", PlastoGram_hmm_Sec_graphpart), "Tat_model" = gsub(".hmm", "", PlastoGram_hmm_Tat_graphpart)),
                            ith_fold = NULL, remove_hmm_files = FALSE)
  ),
  tar_target(
    PlastoGram_higher_level_model_graphpart,
    train_higher_level_model(PlastoGram_predictions_graphpart, data_df_graphpart, PlastoGram_best_architecture_name_graphpart, PlastoGram_best_architecture_graphpart)
  ),
  tar_target(
    PlastoGram_informative_ngrams_graphpart,
    get_all_imp_ngrams(PlastoGram_ngram_models_graphpart)
  ),
  tar_target(
    PlastoGram_evaluation_graphpart,
    predict_with_PlastoGram(PlastoGram_ngram_models_graphpart, 
                            list("Sec_model" = gsub(".hmm", "", PlastoGram_hmm_Sec_graphpart), "Tat_model" = gsub(".hmm", "", PlastoGram_hmm_Tat_graphpart)), 
                            PlastoGram_higher_level_model_graphpart, ngram_matrix_independent_graphpart, sequences_independent_graphpart, data_df_independent_graphpart,
                            PlastoGram_informative_ngrams_graphpart, PlastoGram_best_architecture_graphpart, PlastoGram_best_architecture_name_graphpart)
  ),
  tar_target(
    PlastoGram_OM_IM_model_graphpart,
    train_om_im_model(ngram_matrix_traintest_graphpart, data_df_graphpart)
  ),
  tar_target(
    PlastoGram_evaluation_OM_IM_graphpart,
    evaluate_om_im_model(PlastoGram_OM_IM_model_graphpart, ngram_matrix_independent_graphpart, data_df_independent_graphpart, PlastoGram_evaluation_graphpart)
  ),
  # Results and final models
  tar_target(
    PlastoGram_best_architecture_name,
    arrange(envelope_mean_architecture_performance, desc(mean_kappa))[["model"]][1]
  ),
  tar_target(
    PlastoGram_best_architecture,
    change_model_names(read.csv(paste0(data_path, "Model_architectures_envelope/", gsub("_GLM|_RF", ".csv", PlastoGram_best_architecture_name))))
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
    graphpart_baseline_model_cv_res,
    do_baseline_cv(ngram_matrix_traintest_graphpart, target_dfs_cv_graphpart)
  ),
  tar_target(
    graphpart_baseline_mean_performance,
    get_mean_performance_of_baseline(graphpart_baseline_model_cv_res, paste0(data_path, "Baseline_graphpart_envelope_mean_performance.csv"))
  ),
  tar_target(
    graphpart_baseline_model,
    train_baseline_model(ngram_matrix_traintest_graphpart, data_df_graphpart)
  ),
  tar_target(
    graphpart_baseline_model_evaluation,
    evaluate_baseline_model(graphpart_baseline_model, ngram_matrix_independent_graphpart, data_df_independent_graphpart)
  ),
  # Publication results
  tar_target(
    architecture_plot_data,
    get_architecture_plot_data(envelope_mean_architecture_performance)
  ),
  tar_target(
    graphpart_architecture_plot_data,
    get_architecture_plot_data(mean_architecture_performance_graphpart)
  ),
  tar_target(
    performance_distr_plot,
    get_performance_distr_plot(architecture_plot_data, graphpart_architecture_plot_data, baseline_mean_performance, 
                               graphpart_baseline_mean_performance, paste0(data_path, "Publication_results/"))
  ),
  tar_target(
    performance_parameters_table,
    get_parameters_of_performance_distribution_table(envelope_mean_architecture_performance, paste0(data_path, "Publication_results/"),
                                                     "Performance_distribution_parameters_holdout.csv")
  ),
  tar_target(
    performance_parameters_table_graphpart,
    get_parameters_of_performance_distribution_table(mean_architecture_performance_graphpart, paste0(data_path, "Publication_results/"),
                                                     "Performance_distribution_parameters_partitioning.csv")
  ),
  tar_target(
    best_model_cv_plot,
    get_best_model_cv_plot(architecture_plot_data, graphpart_architecture_plot_data, paste0(data_path, "Publication_results/"))
  ),
  tar_target(
    cv_and_independent_res_table,
    get_cv_and_independent_res_table(architecture_plot_data, PlastoGram_evaluation, paste0(data_path, "Publication_results/"),
                                     "CV+independent_res_holdout.csv")
  ),
  tar_target(
    cv_and_independent_res_table_graphpart,
    get_cv_and_independent_res_table(graphpart_architecture_plot_data, PlastoGram_evaluation_graphpart, paste0(data_path, "Publication_results/"),
                                     "CV+independent_res_partitioning.csv")
  ),
  tar_target(
    om_im_model_cv_performance_table,
    get_om_im_model_cv_res_table(envelope_traintest_ngram_matrix, envelope_target_dfs_cv, paste0(data_path, "Publication_results/"),
                                 "OM_IM_model_cv_res_holdout.csv") 
  ),
  tar_target(
    om_im_model_cv_performance_table_graphpart,
    get_om_im_model_cv_res_table(ngram_matrix_traintest_graphpart, target_dfs_cv_graphpart, paste0(data_path, "Publication_results/"),
                                 "OM_IM_model_cv_res_partitioning.csv") 
  ),
  tar_target(
    om_im_model_mean_cv_performance_table,
    get_mean_performance_of_om_im_models(om_im_model_cv_performance_table,
                                         paste0(data_path, "Publication_results/"),
                                         "OM_IM_model_mean_cv_res_holdout.csv")
  ),
  tar_target(
    om_im_model_mean_cv_performance_table_graphpart,
    get_mean_performance_of_om_im_models(om_im_model_cv_performance_table_graphpart,
                                         paste0(data_path, "Publication_results/"), 
                                         "OM_IM_model_mean_cv_res_partitioning.csv")
  ),
  tar_target(
    physicochemical_prop_plot,
    ggsave(paste0(data_path, "Publication_results/OM_stroma_properties.eps"), 
           get_physicochemical_properties_plot(all_sequences, holdout_target_df, 
                                               colors = c("N_OM" = "#b172d8", "N_IM" = "#7281d8", "N_TM" = "#d87272", "N_S" = "#76d872")), 
           width = 9, height = 2.5)
  ),
  tar_target(
    physicochemical_prop_test,
    get_physicochemical_properties_test(all_sequences, holdout_target_df, paste0(data_path, "Publication_results/")) 
  ),
  tar_target(
    pca_props_plot,
    get_pca_and_props_plot(sequence_list, traintest, traintest_data_df, paste0(data_path, "Publication_results/"))
  ),
  tar_target(
    benchmark_file,
    write_fasta(benchmark, "./data/Independent_dataset.fa")
  ),
  tar_target(
    benchmark_file_graphpart,
    write_fasta(sequences_independent_graphpart, "./data/Independent_dataset_partitioning.fa")
  ),
  tar_target(
    SChloro_benchmark_res_file,
    "./data/schloro_independent_dataset_holdout_res.gff",
    format = "file"
  ),
  tar_target(
    SChloro_benchmark_res_file_graphpart,
    "./data/schloro_independent_dataset_partitioning_res.gff",
    format = "file"
  ),
  tar_target(
    SChloro_benchmark_res,
    read.delim(SChloro_benchmark_res_file, skip = 1, header = FALSE)
  ),
  tar_target(
    SChloro_benchmark_res_graphpart,
    read.delim(SChloro_benchmark_res_file_graphpart, skip = 1, header = FALSE)
  ),
  tar_target(
    benchmark_table,
    write.csv(get_benchmark_res_table(PlastoGram_evaluation, SChloro_benchmark_res), 
              paste0(data_path, "Publication_results/Benchmark_results_holdout.csv"), row.names = FALSE)
  ),
  tar_target(
    benchmark_table_graphpart,
    write.csv(get_benchmark_res_table(PlastoGram_evaluation_graphpart, SChloro_benchmark_res_graphpart), 
              paste0(data_path, "Publication_results/Benchmark_results_partitioning.csv"), row.names = FALSE)
  ),
  tar_target(
    benchmark_plot,
    get_benchmark_res_plot(PlastoGram_evaluation, PlastoGram_evaluation_graphpart, SChloro_benchmark_res, SChloro_benchmark_res_graphpart, 
                           paste0(data_path, "Publication_results/")) 
  ),
  tar_target(
    PlastoGram_model_holdout,
    get_final_plastogram_model(PlastoGram_ngram_models, PlastoGram_higher_level_model, PlastoGram_OM_IM_model, PlastoGram_informative_ngrams)
  ),
  tar_target(
    PlastoGram_model_partitioning,
    get_final_plastogram_model(PlastoGram_ngram_models_graphpart, PlastoGram_higher_level_model_graphpart, PlastoGram_OM_IM_model_graphpart, PlastoGram_informative_ngrams_graphpart)
  ),
  tar_target(
    pairwise_identity,
    calculate_pairwise_identity(traintest, benchmark, sequences_cv_graphpart, sequences_independent_graphpart,
                                envelope_target_df, target_df_graphpart, data_path)
  ),
  tar_target(
    mean_pairwise_identity,
    calculate_mean_pairwise_identity(pairwise_identity, paste0(data_path, "Publication_results/"))
  ),
  tar_target(
    cv_res_statistical_test,
    get_cv_statistical_test(envelope_architectures_performance, PlastoGram_best_architecture_name, "Holdout",
                            architectures_performance_graphpart, PlastoGram_best_architecture_name_graphpart, "Partitioning",
                            paste0(data_path, "Publication_results/"))
  ),
  tar_target(
    feature_size_table,
    get_feature_size_table(PlastoGram_ngram_models, PlastoGram_ngram_models_graphpart, 
                           paste0(data_path, "Publication_results/"))
  ),
  tar_target(
    feature_venn_diagram,
    get_features_venn_diagram(PlastoGram_ngram_models, PlastoGram_ngram_models_graphpart,
                              paste0(data_path, "Publication_results/"))
  )
)

