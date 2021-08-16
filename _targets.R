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

data_path <- "~/Dropbox/Projekty/BioNgramProjects/PlastoGram/"

source("./functions/create_datasets.R")
source("./functions/filter_sequences.R")
source("./functions/ngram_model_functions.R")
source("./functions/create_target_df.R")
source("./functions/profileHMM_functions.R")
source("./functions/evaluate_full_model.R")
source("./functions/ensemble_model_functions.R")

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
    do_cv(ngram_matrix, target_df, "Nuclear_target", 5, 0.001)
  ),
  tar_target(
    Plastid_Nuclear_CV_res_stats, 
    get_cv_res_summary(Plastid_Nuclear_CV, "TRUE")
  ),
  tar_target(
    Plastid_Nuclear_CV_2, 
    do_cv(ngram_matrix, target_df, "Nuclear_target", 5, 0.00001)
  ),
  tar_target(
    Plastid_Nuclear_CV_2_res_stats, 
    get_cv_res_summary(Plastid_Nuclear_CV_2, "TRUE")
  ),
  tar_target(
    Plastid_Nuclear_cw_CV,
    do_cv(ngram_matrix, target_df, "Nuclear_target", 5, 0.001, with_class_weights = TRUE)
  ),
  tar_target(
    Plastid_Nuclear_cw_CV_res_stats,
    get_cv_res_summary(Plastid_Nuclear_cw_CV, "TRUE")
  ),
  tar_target(
    Membrane_Notmembrane_CV,
    do_cv(ngram_matrix, target_df, "Membrane_target", 5, 0.001)
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
    do_cv(ngram_matrix, target_df, "Membrane_target", 5, 0.001, with_class_weights = TRUE)
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
          target_df[which(target_df[["Membrane_target"]] == TRUE & target_df[["Nuclear_target"]] == FALSE),], "IM_target",
          5, 0.001)
  ),
  tar_target(
    Membrane_Nuclear_CV,
    do_cv(ngram_matrix_membrane_nuclear, 
          target_df[which(target_df[["Membrane_target"]] == TRUE & target_df[["Nuclear_target"]] == TRUE),], "Membrane_mc_target",
          5, 0.001, mc = TRUE)
  ),
  tar_target(
    Membrane_Plastid_cw_CV,
    do_cv(ngram_matrix_membrane_plastid, 
          target_df[which(target_df[["Membrane_target"]] == TRUE & target_df[["Nuclear_target"]] == FALSE),], "IM_target",
          5, 0.001, with_class_weights = TRUE)
  ),
  tar_target(
    Membrane_Nuclear_cw_CV,
    do_cv(ngram_matrix_membrane_nuclear, 
          target_df[which(target_df[["Membrane_target"]] == TRUE & target_df[["Nuclear_target"]] == TRUE),], "Membrane_mc_target",
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
          target_df[which(target_df[["Membrane_target"]] == TRUE & target_df[["Nuclear_target"]] == TRUE),], "Membrane_mc_target",
          5, 0.01, mc = TRUE)
  ),
  tar_target(
    Membrane_Nuclear_cw_CV_2,
    do_cv(ngram_matrix_membrane_nuclear, 
          target_df[which(target_df[["Membrane_target"]] == TRUE & target_df[["Nuclear_target"]] == TRUE),], "Membrane_mc_target",
          5, 0.01, mc = TRUE, with_class_weights = TRUE)
  ),
  tar_target(
    Membrane_Nuclear_CV_3,
    do_cv(ngram_matrix_membrane_nuclear, 
          target_df[which(target_df[["Membrane_target"]] == TRUE & target_df[["Nuclear_target"]] == TRUE),], "Membrane_mc_target",
          5, 0.05, mc = TRUE)
  ),
  tar_target(
    Membrane_Nuclear_cw_CV_3,
    do_cv(ngram_matrix_membrane_nuclear, 
          target_df[which(target_df[["Membrane_target"]] == TRUE & target_df[["Nuclear_target"]] == TRUE),], "Membrane_mc_target",
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
    do_cv(ngram_matrix, target_df, "Stroma_target", 5, 0.001)
  ),
  tar_target(
    Stroma_cw_CV,
    do_cv(ngram_matrix, target_df, "Stroma_target", 5, 0.001, with_class_weights = TRUE)
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
          target_df[which(target_df[["Nuclear_target"]] == TRUE & (target_df[["Membrane_mc_target"]] == "OM" | target_df[["Stroma_target"]] == TRUE)),],
          "OM_target", 5, 0.001)
  ),
  tar_target(
    OM_Stroma_CV_res_stats,
    get_cv_res_summary(OM_Stroma_CV, "TRUE")
  ),
  tar_target(
    OM_Stroma_cw_CV,
    do_cv(ngram_matrix[which(target_df[["Nuclear_target"]] == TRUE & (target_df[["Membrane_mc_target"]] == "OM" | target_df[["Stroma_target"]] == TRUE)), ],
          target_df[which(target_df[["Nuclear_target"]] == TRUE & (target_df[["Membrane_mc_target"]] == "OM" | target_df[["Stroma_target"]] == TRUE)), ],
          "OM_target", 5, 0.001, with_class_weights = TRUE)
  ),
  tar_target(
    OM_Stroma_cw_CV_res_stats,
    get_cv_res_summary(OM_Stroma_cw_CV, "TRUE")
  ),
  tar_target(
    OM_Stroma_CV_2,
    do_cv(ngram_matrix[which(target_df[["Nuclear_target"]] == TRUE & (target_df[["Membrane_mc_target"]] == "OM" | target_df[["Stroma_target"]] == TRUE)),],
          target_df[which(target_df[["Nuclear_target"]] == TRUE & (target_df[["Membrane_mc_target"]] == "OM" | target_df[["Stroma_target"]] == TRUE)),],
          "OM_target", 5, 0.01)
  ),
  tar_target(
    OM_Stroma_CV_2_res_stats,
    get_cv_res_summary(OM_Stroma_CV_2, "TRUE")
  ),
  tar_target(
    OM_Stroma_cw_CV_2,
    do_cv(ngram_matrix[which(target_df[["Nuclear_target"]] == TRUE & (target_df[["Membrane_mc_target"]] == "OM" | target_df[["S_target"]] == TRUE)), ],
          target_df[which(target_df[["Nuclear_target"]] == TRUE & (target_df[["Membrane_mc_target"]] == "OM" | target_df[["S_target"]] == TRUE)), ],
          "OM_target", 5, 0.01, with_class_weights = TRUE)
  ),
  tar_target(
    OM_Stroma_cw_CV_2_res_stats,
    get_cv_res_summary(OM_Stroma_cw_CV_2, "TRUE")
  ),
  tar_target(
    P_stroma_CV,
    do_cv(ngram_matrix[which(target_df[["Nuclear_target"]] == FALSE), ],
          target_df[which(target_df[["Nuclear_target"]] == FALSE), ],
          "Stroma_target", 5, 0.001)
  ),
  tar_target(
    P_stroma_cw_CV,
    do_cv(ngram_matrix[which(target_df[["Nuclear_target"]] == FALSE), ],
          target_df[which(target_df[["Nuclear_target"]] == FALSE), ],
          "Stroma_target", 5, 0.001, with_class_weights = TRUE)
  ),
  tar_target(
    P_stroma_CV_2,
    do_cv(ngram_matrix[which(target_df[["Nuclear_target"]] == FALSE), ],
          target_df[which(target_df[["Nuclear_target"]] == FALSE), ],
          "Stroma_target", 5, 0.01)
  ),
  tar_target(
    P_stroma_cw_CV_2,
    do_cv(ngram_matrix[which(target_df[["Nuclear_target"]] == FALSE), ],
          target_df[which(target_df[["Nuclear_target"]] == FALSE ), ],
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
          target_df[which(target_df[["Nuclear_target"]] == TRUE), ],
          "Stroma_target", 5, 0.001)
  ),
  tar_target(
    N_stroma_cw_CV,
    do_cv(ngram_matrix[which(target_df[["Nuclear_target"]] == TRUE), ],
          target_df[which(target_df[["Nuclear_target"]] == TRUE), ],
          "Stroma_target", 5, 0.001, with_class_weights = TRUE)
  ),
  tar_target(
    N_stroma_CV_2,
    do_cv(ngram_matrix[which(target_df[["Nuclear_target"]] == TRUE), ],
          target_df[which(target_df[["Nuclear_target"]] == TRUE), ],
          "Stroma_target", 5, 0.01)
  ),
  tar_target(
    N_stroma_cw_CV_2,
    do_cv(ngram_matrix[which(target_df[["Nuclear_target"]] == TRUE), ],
          target_df[which(target_df[["Nuclear_target"]] == TRUE), ],
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
          target_df[which(target_df[["Membrane_target"]] == TRUE & target_df[["Nuclear_target"]] == TRUE),], "OM_target",
          5, 0.001)
  ),
  tar_target(
    OM_other_membrane_cw_CV,
    do_cv(ngram_matrix_membrane_nuclear, 
          target_df[which(target_df[["Membrane_target"]] == TRUE & target_df[["Nuclear_target"]] == TRUE),], "OM_target",
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
          target_df[which(target_df[["Membrane_target"]] == TRUE & target_df[["Nuclear_target"]] == TRUE),], "OM_target",
          5, 0.01)
  ),
  tar_target(
    OM_other_membrane_cw_CV_2,
    do_cv(ngram_matrix_membrane_nuclear, 
          target_df[which(target_df[["Membrane_target"]] == TRUE & target_df[["Nuclear_target"]] == TRUE),], "OM_target",
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
          target_df[which(target_df[["Membrane_target"]] == TRUE & target_df[["Nuclear_target"]] == TRUE),], "IM_target",
          5, 0.001)
  ),
  tar_target(
    IM_other_membrane_cw_CV,
    do_cv(ngram_matrix_membrane_nuclear, 
          target_df[which(target_df[["Membrane_target"]] == TRUE & target_df[["Nuclear_target"]] == TRUE),], "IM_target",
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
          target_df[which(target_df[["Membrane_target"]] == TRUE & target_df[["Nuclear_target"]] == TRUE),], "IM_target",
          5, 0.01)
  ),
  tar_target(
    IM_other_membrane_cw_CV_2,
    do_cv(ngram_matrix_membrane_nuclear, 
          target_df[which(target_df[["Membrane_target"]] == TRUE & target_df[["Nuclear_target"]] == TRUE),], "IM_target",
          5, 0.01, with_class_weights = TRUE)
  ),
  tar_target(
    IM_other_membrane_CV_2_res_stats,
    get_cv_res_summary(IM_other_membrane_CV_2, "TRUE")
  ),
  tar_target(
    IM_other_membrane_cw_CV_2_res_stats,
    get_cv_res_summary(IM_other_membrane_cw_CV_2, "TRUE")
  ),
  tar_target(
    TM_other_membrane_CV,
    do_cv(ngram_matrix_membrane_nuclear, 
          target_df[which(target_df[["Membrane_target"]] == TRUE & target_df[["Nuclear_target"]] == TRUE),], "TM_target",
          5, 0.001)
  ),
  tar_target(
    TM_other_membrane_cw_CV,
    do_cv(ngram_matrix_membrane_nuclear, 
          target_df[which(target_df[["Membrane_target"]] == TRUE & target_df[["Nuclear_target"]] == TRUE),], "TM_target",
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
          target_df[which(target_df[["Membrane_target"]] == TRUE & target_df[["Nuclear_target"]] == TRUE),], "TM_target",
          5, 0.01)
  ),
  tar_target(
    TM_other_membrane_cw_CV_2,
    do_cv(ngram_matrix_membrane_nuclear, 
          target_df[which(target_df[["Membrane_target"]] == TRUE & target_df[["Nuclear_target"]] == TRUE),], "TM_target",
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
  )
)
