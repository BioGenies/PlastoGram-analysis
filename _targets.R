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

data_path <- "~/Dropbox/Projekty/BioNgramProjects/PlastoGram/"

source("./functions/create_datasets.R")
source("./functions/filter_sequences.R")
source("./functions/ngram_model_functions.R")
source("./functions/create_target_df.R")
source("./functions/profileHMM_functions.R")
source("./functions/evaluate_full_model.R")
source("./functions/ensemble_model_functions.R")
source("./functions/generate_architectures.R")

set.seed(108567)

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
    sequences,
    create_datasets(all_seqs_file,
                    annotations_file,
                    paste0(data_path, "Sequences/")),
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
    N_seqs,
    c(N_OM_seqs, N_IM_seqs, N_S_seqs, N_TM_seqs, N_TL_SEC_seqs, N_TL_TAT_seqs)
  ),
  tar_target(
    P_seqs, 
    c(P_IM_seqs, P_S_seqs, P_TM_seqs)
  ),
  tar_target(
    target_df, 
    create_target_df(annotations_file,
                     c(N_seqs, P_seqs))),
  tar_target(
    ngram_matrix, 
    create_ngram_matrix(c(N_seqs, P_seqs), target_df)
  ),
  tar_target(
    Plastid_Nuclear_CV, 
    do_cv(ngram_matrix, data_df, "Nuclear_target", 5, 0.001)
  ),
  tar_target(
    Plastid_Nuclear_CV_res_stats, 
    get_cv_res_summary(Plastid_Nuclear_CV, "TRUE")
  ),
  tar_target(
    Plastid_Nuclear_cw_CV,
    do_cv(ngram_matrix, data_df, "Nuclear_target", 5, 0.001, with_class_weights = TRUE)
  ),
  tar_target(
    Plastid_Nuclear_cw_CV_res_stats,
    get_cv_res_summary(Plastid_Nuclear_cw_CV, "TRUE")
  ),
  tar_target(
    Membrane_Notmembrane_CV,
    do_cv(ngram_matrix, data_df, "Membrane_target", 5, 0.001)
  ),
  tar_target(
    ngram_matrix_membrane,
    ngram_matrix[which(target_df[["Membrane_target"]] == TRUE),]
  ),
  tar_target(
    target_df_membrane,
    filter(target_df, Membrane_target == TRUE)
  ),
  tar_target(
    Membrane_Notmembrane_CV_res_stats,
    get_cv_res_summary(Membrane_Notmembrane_CV, "TRUE")
  ),
  tar_target(
    Membrane_Notmembrane_cw_CV,
    do_cv(ngram_matrix, data_df, "Membrane_target", 5, 0.001, with_class_weights = TRUE)
  ),
  tar_target(
    Membrane_Notmembrane_cw_CV_res_stats,
    get_cv_res_summary(Membrane_Notmembrane_cw_CV, "TRUE")
  ),
  tar_target(
    ngram_matrix_membrane_plastid,
    ngram_matrix[which(target_df[["Membrane_target"]] == TRUE & target_df[["Nuclear_target"]] == FALSE),]
  ),
  tar_target(
    ngram_matrix_membrane_nuclear,
    ngram_matrix[which(target_df[["Membrane_target"]] == TRUE & target_df[["Nuclear_target"]] == TRUE),]
  ),
  tar_target(
    Membrane_Plastid_CV,
    do_cv(ngram_matrix_membrane_plastid, 
          data_df[which(target_df[["Membrane_target"]] == TRUE & target_df[["Nuclear_target"]] == FALSE),], "IM_target",
          5, 0.001)
  ),
  tar_target(
    Membrane_Nuclear_CV,
    do_cv(ngram_matrix_membrane_nuclear, 
          data_df[which(target_df[["Membrane_target"]] == TRUE & target_df[["Nuclear_target"]] == TRUE),], "Membrane_mc_target",
          5, 0.001, mc = TRUE)
  ),
  tar_target(
    Membrane_Plastid_cw_CV,
    do_cv(ngram_matrix_membrane_plastid, 
          data_df[which(target_df[["Membrane_target"]] == TRUE & target_df[["Nuclear_target"]] == FALSE),], "IM_target",
          5, 0.001, with_class_weights = TRUE)
  ),
  tar_target(
    Membrane_Nuclear_cw_CV,
    do_cv(ngram_matrix_membrane_nuclear, 
          data_df[which(target_df[["Membrane_target"]] == TRUE & target_df[["Nuclear_target"]] == TRUE),], "Membrane_mc_target",
          5, 0.001, mc = TRUE, with_class_weights = TRUE)
  ),
  tar_target(
    Membrane_Plastid_CV_res_stats,
    get_cv_res_summary(Membrane_Plastid_CV, "TRUE")
  ),
  tar_target(
    Membrane_Nuclear_CV_res_stats,
    get_cv_res_summary_mc(Membrane_Nuclear_CV)
  ),
  tar_target(
    Membrane_Plastid_cw_CV_res_stats,
    get_cv_res_summary(Membrane_Plastid_cw_CV, "TRUE")
  ),
  tar_target(
    Membrane_Nuclear_cw_CV_res_stats,
    get_cv_res_summary_mc(Membrane_Nuclear_cw_CV)
  ),
  tar_target(
    Membrane_Nuclear_CV_2,
    do_cv(ngram_matrix_membrane_nuclear, 
          data_df[which(target_df[["Membrane_target"]] == TRUE & target_df[["Nuclear_target"]] == TRUE),], "Membrane_mc_target",
          5, 0.01, mc = TRUE)
  ),
  tar_target(
    Membrane_Nuclear_cw_CV_2,
    do_cv(ngram_matrix_membrane_nuclear, 
          data_df[which(target_df[["Membrane_target"]] == TRUE & target_df[["Nuclear_target"]] == TRUE),], "Membrane_mc_target",
          5, 0.01, mc = TRUE, with_class_weights = TRUE)
  ),
  tar_target(
    Membrane_Nuclear_CV_3,
    do_cv(ngram_matrix_membrane_nuclear, 
          data_df[which(target_df[["Membrane_target"]] == TRUE & target_df[["Nuclear_target"]] == TRUE),], "Membrane_mc_target",
          5, 0.05, mc = TRUE)
  ),
  tar_target(
    Membrane_Nuclear_cw_CV_3,
    do_cv(ngram_matrix_membrane_nuclear, 
          data_df[which(target_df[["Membrane_target"]] == TRUE & target_df[["Nuclear_target"]] == TRUE),], "Membrane_mc_target",
          5, 0.05, mc = TRUE, with_class_weights = TRUE)
  ),
  tar_target(
    Membrane_Nuclear_CV_2_res_stats,
    get_cv_res_summary_mc(Membrane_Nuclear_CV_2)
  ),
  tar_target(
    Membrane_Nuclear_cw_CV_2_res_stats,
    get_cv_res_summary_mc(Membrane_Nuclear_cw_CV_2)
  ),
  tar_target(
    Membrane_Nuclear_CV_3_res_stats,
    get_cv_res_summary_mc(Membrane_Nuclear_CV_3)
  ),
  tar_target(
    Membrane_Nuclear_cw_CV_3_res_stats,
    get_cv_res_summary_mc(Membrane_Nuclear_cw_CV_3)
  ),
  tar_target(
    Stroma_CV,
    do_cv(ngram_matrix, data_df, "Stroma_target", 5, 0.001)
  ),
  tar_target(
    Stroma_cw_CV,
    do_cv(ngram_matrix, data_df, "Stroma_target", 5, 0.001, with_class_weights = TRUE)
  ),
  tar_target(
    Stroma_CV_res_stats,
    get_cv_res_summary(Stroma_CV, "TRUE")
  ),
  tar_target(
    Stroma_cw_CV_res_stats,
    get_cv_res_summary(Stroma_cw_CV, "TRUE")
  ),
  tar_target(
    OM_Stroma_CV,
    do_cv(ngram_matrix[which(target_df[["Nuclear_target"]] == TRUE & (target_df[["Membrane_mc_target"]] == "OM" | target_df[["Stroma_target"]] == TRUE)),], 
          data_df[which(target_df[["Nuclear_target"]] == TRUE & (target_df[["Membrane_mc_target"]] == "OM" | target_df[["Stroma_target"]] == TRUE)),],
          "OM_target", 5, 0.001)
  ),
  tar_target(
    OM_Stroma_CV_res_stats,
    get_cv_res_summary(OM_Stroma_CV, "TRUE")
  ),
  tar_target(
    OM_Stroma_cw_CV,
    do_cv(ngram_matrix[which(target_df[["Nuclear_target"]] == TRUE & (target_df[["Membrane_mc_target"]] == "OM" | target_df[["Stroma_target"]] == TRUE)), ],
          data_df[which(target_df[["Nuclear_target"]] == TRUE & (target_df[["Membrane_mc_target"]] == "OM" | target_df[["Stroma_target"]] == TRUE)), ],
          "OM_target", 5, 0.001, with_class_weights = TRUE)
  ),
  tar_target(
    OM_Stroma_cw_CV_res_stats,
    get_cv_res_summary(OM_Stroma_cw_CV, "TRUE")
  ),
  tar_target(
    OM_Stroma_CV_2,
    do_cv(ngram_matrix[which(target_df[["Nuclear_target"]] == TRUE & (target_df[["Membrane_mc_target"]] == "OM" | target_df[["Stroma_target"]] == TRUE)),], 
          data_df[which(target_df[["Nuclear_target"]] == TRUE & (target_df[["Membrane_mc_target"]] == "OM" | target_df[["Stroma_target"]] == TRUE)),],
          "OM_target", 5, 0.01)
  ),
  tar_target(
    OM_Stroma_CV_2_res_stats,
    get_cv_res_summary(OM_Stroma_CV_2, "TRUE")
  ),
  tar_target(
    OM_Stroma_cw_CV_2,
    do_cv(ngram_matrix[which(target_df[["Nuclear_target"]] == TRUE & (target_df[["Membrane_mc_target"]] == "OM" | target_df[["Stroma_target"]] == TRUE)), ], 
          data_df[which(target_df[["Nuclear_target"]] == TRUE & (target_df[["Membrane_mc_target"]] == "OM" | target_df[["Stroma_target"]] == TRUE)), ],
          "OM_target", 5, 0.01, with_class_weights = TRUE)
  ),
  tar_target(
    OM_Stroma_cw_CV_2_res_stats,
    get_cv_res_summary(OM_Stroma_cw_CV_2, "TRUE")
  ),
  tar_target(
    P_stroma_CV,
    do_cv(ngram_matrix[which(target_df[["Nuclear_target"]] == FALSE), ], 
          data_df[which(target_df[["Nuclear_target"]] == FALSE), ],
          "Stroma_target", 5, 0.001)
  ),
  tar_target(
    P_stroma_cw_CV,
    do_cv(ngram_matrix[which(target_df[["Nuclear_target"]] == FALSE), ], 
          data_df[which(target_df[["Nuclear_target"]] == FALSE), ],
          "Stroma_target", 5, 0.001, with_class_weights = TRUE)
  ),
  tar_target(
    P_stroma_CV_2,
    do_cv(ngram_matrix[which(target_df[["Nuclear_target"]] == FALSE), ], 
          data_df[which(target_df[["Nuclear_target"]] == FALSE), ],
          "Stroma_target", 5, 0.01)
  ),
  tar_target(
    P_stroma_cw_CV_2,
    do_cv(ngram_matrix[which(target_df[["Nuclear_target"]] == FALSE), ], 
          data_df[which(target_df[["Nuclear_target"]] == FALSE), ],
          "Stroma_target", 5, 0.01, with_class_weights = TRUE)
  ),
  tar_target(
    P_stroma_CV_res_stats,
    get_cv_res_summary(P_stroma_CV, "TRUE")
  ),
  tar_target(
    P_stroma_CV_2_res_stats,
    get_cv_res_summary(P_stroma_CV_2, "TRUE")
  ),
  tar_target(
    P_stroma_cw_CV_res_stats,
    get_cv_res_summary(P_stroma_cw_CV, "TRUE")
  ),
  tar_target(
    P_stroma_cw_CV_2_res_stats,
    get_cv_res_summary(P_stroma_cw_CV_2, "TRUE")
  ),
  tar_target(
    N_stroma_CV,
    do_cv(ngram_matrix[which(target_df[["Nuclear_target"]] == TRUE), ], 
          data_df[which(target_df[["Nuclear_target"]] == TRUE), ],
          "Stroma_target", 5, 0.001)
  ),
  tar_target(
    N_stroma_cw_CV,
    do_cv(ngram_matrix[which(target_df[["Nuclear_target"]] == TRUE), ], 
          data_df[which(target_df[["Nuclear_target"]] == TRUE), ],
          "Stroma_target", 5, 0.001, with_class_weights = TRUE)
  ),
  tar_target(
    N_stroma_CV_2,
    do_cv(ngram_matrix[which(target_df[["Nuclear_target"]] == TRUE), ], 
          data_df[which(target_df[["Nuclear_target"]] == TRUE), ],
          "Stroma_target", 5, 0.01)
  ),
  tar_target(
    N_stroma_cw_CV_2,
    do_cv(ngram_matrix[which(target_df[["Nuclear_target"]] == TRUE), ], 
          data_df[which(target_df[["Nuclear_target"]] == TRUE), ],
          "Stroma_target", 5, 0.01, with_class_weights = TRUE)
  ),
  tar_target(
    N_stroma_CV_res_stats,
    get_cv_res_summary(N_stroma_CV, "TRUE")
  ),
  tar_target(
    N_stroma_CV_2_res_stats,
    get_cv_res_summary(N_stroma_CV_2, "TRUE")
  ),
  tar_target(
    N_stroma_cw_CV_res_stats,
    get_cv_res_summary(N_stroma_cw_CV, "TRUE")
  ),
  tar_target(
    N_stroma_cw_CV_2_res_stats,
    get_cv_res_summary(N_stroma_cw_CV_2, "TRUE")
  ),
  tar_target(
    OM_other_membrane_CV,
    do_cv(ngram_matrix_membrane_nuclear, 
          data_df[which(target_df[["Membrane_target"]] == TRUE & target_df[["Nuclear_target"]] == TRUE), ], "OM_target",
          5, 0.001)
  ),
  tar_target(
    OM_other_membrane_cw_CV,
    do_cv(ngram_matrix_membrane_nuclear, 
          data_df[which(target_df[["Membrane_target"]] == TRUE & target_df[["Nuclear_target"]] == TRUE), ], "OM_target",
          5, 0.001, with_class_weights = TRUE)
  ),
  tar_target(
    OM_other_membrane_CV_res_stats,
    get_cv_res_summary(OM_other_membrane_CV, "TRUE")
  ),
  tar_target(
    OM_other_membrane_cw_CV_res_stats,
    get_cv_res_summary(OM_other_membrane_cw_CV, "TRUE")
  ),
  tar_target(
    OM_other_membrane_CV_2,
    do_cv(ngram_matrix_membrane_nuclear, 
          data_df[which(target_df[["Membrane_target"]] == TRUE & target_df[["Nuclear_target"]] == TRUE), ], "OM_target",
          5, 0.01)
  ),
  tar_target(
    OM_other_membrane_cw_CV_2,
    do_cv(ngram_matrix_membrane_nuclear, 
          data_df[which(target_df[["Membrane_target"]] == TRUE & target_df[["Nuclear_target"]] == TRUE), ], "OM_target",
          5, 0.01, with_class_weights = TRUE)
  ),
  tar_target(
    OM_other_membrane_CV_2_res_stats,
    get_cv_res_summary(OM_other_membrane_CV_2, "TRUE")
  ),
  tar_target(
    OM_other_membrane_cw_CV_2_res_stats,
    get_cv_res_summary(OM_other_membrane_cw_CV_2, "TRUE")
  ),
  tar_target(
    IM_other_membrane_CV,
    do_cv(ngram_matrix_membrane_nuclear, 
          data_df[which(target_df[["Membrane_target"]] == TRUE & target_df[["Nuclear_target"]] == TRUE), ], "IM_target",
          5, 0.001)
  ),
  tar_target(
    IM_other_membrane_cw_CV,
    do_cv(ngram_matrix_membrane_nuclear, 
          data_df[which(target_df[["Membrane_target"]] == TRUE & target_df[["Nuclear_target"]] == TRUE),], "IM_target",
          5, 0.001, with_class_weights = TRUE)
  ),
  tar_target(
    IM_other_membrane_CV_res_stats,
    get_cv_res_summary(IM_other_membrane_CV, "TRUE")
  ),
  tar_target(
    IM_other_membrane_cw_CV_res_stats,
    get_cv_res_summary(IM_other_membrane_cw_CV, "TRUE")
  ),
  tar_target(
    IM_other_membrane_CV_2,
    do_cv(ngram_matrix_membrane_nuclear, 
          data_df[which(target_df[["Membrane_target"]] == TRUE & target_df[["Nuclear_target"]] == TRUE),], "IM_target",
          5, 0.01)
  ),
  tar_target(
    IM_other_membrane_CV_2_res_stats,
    get_cv_res_summary(IM_other_membrane_CV_2, "TRUE")
  ),
  tar_target(
    TM_other_membrane_CV,
    do_cv(ngram_matrix_membrane_nuclear, 
          data_df[which(target_df[["Membrane_target"]] == TRUE & target_df[["Nuclear_target"]] == TRUE),], "TM_target",
          5, 0.001)
  ),
  tar_target(
    TM_other_membrane_cw_CV,
    do_cv(ngram_matrix_membrane_nuclear, 
          data_df[which(target_df[["Membrane_target"]] == TRUE & target_df[["Nuclear_target"]] == TRUE),], "TM_target",
          5, 0.001, with_class_weights = TRUE)
  ),
  tar_target(
    TM_other_membrane_CV_res_stats,
    get_cv_res_summary(TM_other_membrane_CV, "TRUE")
  ),
  tar_target(
    TM_other_membrane_cw_CV_res_stats,
    get_cv_res_summary(TM_other_membrane_cw_CV, "TRUE")
  ),
  tar_target(
    TM_other_membrane_CV_2,
    do_cv(ngram_matrix_membrane_nuclear, 
          data_df[which(target_df[["Membrane_target"]] == TRUE & target_df[["Nuclear_target"]] == TRUE),], "TM_target",
          5, 0.01)
  ),
  tar_target(
    TM_other_membrane_cw_CV_2,
    do_cv(ngram_matrix_membrane_nuclear, 
          data_df[which(target_df[["Membrane_target"]] == TRUE & target_df[["Nuclear_target"]] == TRUE),],  "TM_target",
          5, 0.01, with_class_weights = TRUE)
  ),
  tar_target(
    TM_other_membrane_CV_2_res_stats,
    get_cv_res_summary(TM_other_membrane_CV_2, "TRUE")
  ),
  tar_target(
    TM_other_membrane_cw_CV_2_res_stats,
    get_cv_res_summary(TM_other_membrane_cw_CV_2, "TRUE")
  ),
  tar_target(
    TL_CV,
    do_hmmer_cv(N_TL_TAT_seqs, N_TL_SEC_seqs)
  ),
  tar_target(
    TL_TAT_CV_res_stats,
    summarise_hmmer_results(TL_CV, 'Tat')
  ),
  tar_target(
    TL_SEC_CV_res_stats,
    summarise_hmmer_results(TL_CV, 'Sec')
  ),
  tar_target(
    data_df,
    create_folds(target_df, 5)
  ),
  tar_target(
    model_variants,
    get_model_variants()
  ),
  tar_target(
    model_dat_file,
    paste0(data_path, "PlastoGram_models_info.csv"),
    format = "file"
  ),
  tar_target(
    model_dat,
    read.csv(model_dat_file)
  ),
  tar_target(
    test_filtering_options_file,
    paste0(data_path, "PlastoGram_models_results_filtering.csv"),
    format = "file"
  ),
  tar_target(
    filtering_df,
    read.csv(test_filtering_options_file)
  ),
  tar_target(
    architectures,
    generate_all_architectures(model_variants = model_variants,
                               smote_models = c("OM_Stroma_model", "Nuclear_membrane_model",
                                                "N_OM_model", "N_IM_model", "N_TM_model"),
                               sequence_models = c("Sec_model", "Tat_model"),
                               model_dat = model_dat,
                               filtering_df = filtering_df,
                               output_dir = paste0(data_path, "Model_architectures/"))
  ),
  tar_target(  
    Plastid_Nuclear_FCBF_CV_res,
    test_fcbf(ngram_matrix, data_df, c(0.05, 0.025, 0.01, 0.005, 0.001), "Nuclear_target")
  ),
  tar_target(
    Membrane_Notmembrane_FCBF_CV_res,
    test_fcbf(ngram_matrix, data_df, c(0.05, 0.025, 0.01, 0.005, 0.001), "Membrane_target")
  ),
  tar_target(
    Nuclear_Membrane_FCBF_CV_res,
    test_fcbf(ngram_matrix_membrane_nuclear, 
              data_df[which(target_df[["Membrane_target"]] == TRUE & target_df[["Nuclear_target"]] == TRUE),], 
              c(0.05, 0.025, 0.01, 0.005, 0.001), "Membrane_mc_target", mc = TRUE)
  ),
  tar_target(
    Plastid_Membrane_FCBF_CV_res,
    test_fcbf(ngram_matrix_membrane_plastid, 
              data_df[which(target_df[["Membrane_target"]] == TRUE & target_df[["Nuclear_target"]] == FALSE),], 
              c(0.05, 0.025, 0.01, 0.005, 0.001), "IM_target")
  ),
  tar_target(
    Stroma_FCBF_CV_res,
    test_fcbf(ngram_matrix, data_df, c(0.05, 0.025, 0.01, 0.005, 0.001), "Stroma_target")
  ),
  tar_target(
    OM_Stroma_FCBF_CV_res,
    test_fcbf(ngram_matrix[which(target_df[["Nuclear_target"]] == TRUE & (target_df[["Membrane_mc_target"]] == "OM" | target_df[["Stroma_target"]] == TRUE)),], 
              data_df[which(target_df[["Nuclear_target"]] == TRUE & (target_df[["Membrane_mc_target"]] == "OM" | target_df[["Stroma_target"]] == TRUE)),],
              c(0.05, 0.025, 0.01, 0.005, 0.001), "OM_target")
  ),
  tar_target(
    P_Stroma_FCBF_CV_res,
    test_fcbf(ngram_matrix[which(target_df[["Nuclear_target"]] == FALSE), ], 
              data_df[which(target_df[["Nuclear_target"]] == FALSE), ],
              c(0.05, 0.025, 0.01, 0.005, 0.001), "Stroma_target")
  ),
  tar_target(
    N_Stroma_FCBF_CV_res,
    test_fcbf(ngram_matrix[which(target_df[["Nuclear_target"]] == TRUE), ], 
              data_df[which(target_df[["Nuclear_target"]] == TRUE), ],
              c(0.05, 0.025, 0.01, 0.005, 0.001), "Stroma_target")
  ),
  tar_target(
    OM_other_membrane_FCBF_CV_res,
    test_fcbf(ngram_matrix_membrane_nuclear, 
              data_df[which(target_df[["Membrane_target"]] == TRUE & target_df[["Nuclear_target"]] == TRUE), ], 
              c(0.05, 0.025, 0.01, 0.005, 0.001), "OM_target")
  ),
  tar_target(
    IM_other_membrane_FCBF_CV_res,
    test_fcbf(ngram_matrix_membrane_nuclear, 
              data_df[which(target_df[["Membrane_target"]] == TRUE & target_df[["Nuclear_target"]] == TRUE), ], 
              c(0.05, 0.025, 0.01, 0.005, 0.001), "IM_target")
  ),
  tar_target(
    TM_other_membrane_FCBF_CV_res,
    test_fcbf(ngram_matrix_membrane_nuclear, 
              data_df[which(target_df[["Membrane_target"]] == TRUE & target_df[["Nuclear_target"]] == TRUE), ], 
              c(0.05, 0.025, 0.01, 0.005, 0.001), "TM_target")
  ),
  tar_target(
    data_df_final,
    create_folds(target_df, 10)
  ),
  tar_target(
    all_models_predictions,
    get_all_models_predictions(ngram_matrix, c(N_seqs, P_seqs), data_df_final, model_dat, data_path)
  ),
  tar_target(
    architecture_files,
    list.files(paste0(data_path, "Model_architectures/"), full.names = TRUE)
  ),
  tar_target(
    architecture_results,
    generate_results_for_architectures(architecture_files,
                                       paste0(data_path, "All_models_predictions.csv"),
                                       paste0(data_path, "Model_architectures_results/"),
                                       data_df_final)
  ),
  tar_target(
    architecture_results_files,
    list.files(paste0(data_path, "Model_architectures_results"), full.names = TRUE)
  ),
  tar_target(
    architectures_performance,
    evaluate_all_architectures(architecture_results_files,
                               paste0(data_path, "Architectures_performance.csv"))
  ),
  tar_target(
    mean_architecture_performance,
    get_mean_performance_of_architectures(architectures_performance,
                                          paste0(data_path, "Architectures_mean_performance.csv"))
  )
)

