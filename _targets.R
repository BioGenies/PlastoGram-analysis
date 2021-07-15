library(dplyr)
library(biogram)
library(targets)
library(seqR)
library(cvTools)
library(ranger)
library(readxl)
library(measures)
library(rhmmer)

data_path <- "~/Dropbox/Projekty/BioNgramProjects/PlastoGram/"

source("./functions/filter_sequences.R")
source("./functions/ngram_model_functions.R")
source("./functions/create_target_df.R")
source("./functions/profileHMM_functions.R")
source("./functions/evaluate_full_model.R")

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
    Membrane_Notmembrane_CV,
    do_cv(ngram_matrix, target_df, "membrane_all_target", 5, 0.001)
  ),
  tar_target(
    ngram_matrix_membrane,
    ngram_matrix[which(target_df[["membrane_all_target"]] == TRUE),]
  ),
  tar_target(
    target_df_membrane,
    filter(target_df, membrane_all_target == TRUE)
  ),
  tar_target(
    Membrane_Notmembrane_CV_res_stats,
    get_cv_res_summary(Membrane_Notmembrane_CV, "TRUE")
  ),
  tar_target(
    ngram_matrix_membrane_plastid,
    ngram_matrix[which(target_df[["membrane_all_target"]] == TRUE & target_df[["NP_target"]] == FALSE),]
  ),
  tar_target(
    ngram_matrix_membrane_nuclear,
    ngram_matrix[which(target_df[["membrane_all_target"]] == TRUE & target_df[["NP_target"]] == TRUE),]
  ),
  tar_target(
    Membrane_Plastid_CV,
    do_cv(ngram_matrix_membrane_plastid, 
          target_df[which(target_df[["membrane_all_target"]] == TRUE & target_df[["NP_target"]] == FALSE),], "membrane_IM_target",
          5, 0.001)
  ),
  tar_target(
    Membrane_Nuclear_CV,
    do_cv(ngram_matrix_membrane_nuclear, 
          target_df[which(target_df[["membrane_all_target"]] == TRUE & target_df[["NP_target"]] == TRUE),], "membrane_target",
          5, 0.001, mc = TRUE)
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
    full_model_CV,
    evaluate_full_model(ngram_matrix, target_df, 5)
  )
)
