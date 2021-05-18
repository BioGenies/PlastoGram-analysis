library(dplyr)
library(biogram)
library(targets)
library(seqR)
library(cvTools)
library(ranger)
library(readxl)
library(measures)

data_path <- "~/Dropbox/Projekty/BioNgramProjects/PlastoGram/"

source("./functions/filter_sequences.R")
source("./functions/ngram_model_functions.R")
source("./functions/create_target_df.R")

set.seed(108567)

list(
  tar_target(
    N_OM_seqs,
    filter_with_cdhit(
      filter_nonstandard_aa(read_fasta(paste0(data_path, "Sequences/N_OM.fa"))), 
                            threshold = 0.9)
  ),
  tar_target(
    N_IM_seqs,
    filter_with_cdhit(
      filter_nonstandard_aa(read_fasta(paste0(data_path, "Sequences/N_IM.fa"))), 
                            threshold = 0.9)
  ),
  tar_target(
    P_IM_seqs,
    filter_with_cdhit(
      filter_nonstandard_aa(read_fasta(paste0(data_path, "Sequences/P_IM.fa"))), 
                            threshold = 0.9)
  ),
  tar_target(
    N_S_seqs, 
    filter_with_cdhit(
      filter_nonstandard_aa(read_fasta(paste0(data_path, "Sequences/N_S.fa"))), 
                            threshold = 0.9)
  ),
  tar_target(
    P_S_seqs, 
    filter_with_cdhit(
      filter_nonstandard_aa(read_fasta(paste0(data_path, "Sequences/P_S.fa"))),
                            threshold = 0.9)
  ),
  tar_target(
    N_TM_seqs, 
    filter_with_cdhit(
      filter_nonstandard_aa(read_fasta(paste0(data_path, "Sequences/N_TM.fa"))),
                            threshold = 0.9)
  ),
  tar_target(
    P_TM_seqs, 
    filter_with_cdhit(
      filter_nonstandard_aa(read_fasta(paste0(data_path, "Sequences/P_TM.fa"))),
                            threshold = 0.9)
  ),
  tar_target(
    N_TL_SEC_seqs, 
    filter_with_cdhit(
      filter_nonstandard_aa(read_fasta(paste0(data_path, "Sequences/N_TL_SEC.fa"))),
                            threshold = 0.9)
  ),
  tar_target(
    N_TL_TAT_seqs, 
    filter_with_cdhit(
      filter_nonstandard_aa(read_fasta(paste0(data_path, "Sequences/N_TL_TAT.fa"))),
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
    create_target_df(paste0(data_path, "Annotations/All_proteins.xlsx"),
                     c(N_seqs, P_seqs))),
  tar_target(
    ngram_matrix, 
    create_ngram_matrix(c(N_seqs, P_seqs), target_df)
  ),
  tar_target(
    Plastid_Nuclear_CV, 
    do_cv(ngram_matrix, target_df, "NP_target", 5, 0.001)
  ),
  tar_target(
    Plastid_Nuclear_CV_res_stats, 
    get_cv_res_summary(Plastid_Nuclear_CV, "TRUE")
  ),
  tar_target(
    Plastid_Nuclear_CV_2, 
    do_cv(ngram_matrix, target_df, "NP_target", 5, 0.00001)
  ),
  tar_target(
    Plastid_Nuclear_CV_2_res_stats, 
    get_cv_res_summary(Plastid_Nuclear_CV_2, "TRUE")
  ),
  tar_target(
    Stroma_CV,
    do_cv(ngram_matrix, target_df, "S_target", 5, 0.001)
  ),
  tar_target(
    Stroma_CV_res_stats,
    get_cv_res_summary(Stroma_CV, "TRUE")
  ),
  tar_target(
    Stroma_CV_2,
    do_cv(ngram_matrix, target_df, "S_target", 5, 0.00001)
  ),
  tar_target(
    Stroma_CV_2_res_stats,
    get_cv_res_summary(Stroma_CV_2, "TRUE")
  ),
  tar_target(
    Membrane_imp_ngrams,
    get_imp_ngrams_mc(ngram_matrix, target_df, "membrane_target", 0.001)
  ),
  tar_target(
    Membrane_CV,
    do_cv(ngram_matrix, target_df, "membrane_target", 5, 0.001, mc = TRUE)
  ),
  tar_target(
    Membrane_CV_res_stats,
    get_cv_res_summary_mc(Membrane_CV)
  ),
  tar_target(
    case_weights,
    get_case_weights(target_df)
  ),
  tar_target(
    Membrane_CV_class_weights,
    do_cv(ngram_matrix, target_df, "membrane_target", 3, 0.001, mc = TRUE, 
          class_weights = c(0.21, 0.75, 0.01, 0.03))
  ),
  tar_target(
    Membrane_CV_class_and_case_weights,
    do_cv(ngram_matrix, target_df, "membrane_target", 3, 0.001, mc = TRUE, 
          class_weights = c(0.21, 0.75, 0.01, 0.03), case_weights = case_weights)
  ),
  tar_target(
    Membrane_CV_class_weights_res_stats,
    get_cv_res_summary_mc(Membrane_CV_class_weights)
  ),
  tar_target(
    Membrane_CV_class_and_case_weights_res_stats,
    get_cv_res_summary_mc(Membrane_CV_class_and_case_weights)
  )
)
