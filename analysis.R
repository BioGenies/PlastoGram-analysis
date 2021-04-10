library(dplyr)
library(biogram)
library(drake)
library(seqR)
library(cvTools)
library(ranger)

data_path <- "~/Dropbox/Projekty/BioNgramProjects/PlastoGram/"

source("./functions/do_cdhit.R")
source("./functions/ngram_model_functions.R")

analysis_PlastoGram <- drake_plan(N_OM_seqs = filter_with_cdhit(read_fasta(paste0(data_path, "Sequences/N_OM.fa")),
                                                                threshold = 0.9),
                                  N_IM_seqs = filter_with_cdhit(read_fasta(paste0(data_path, "Sequences/N_IM.fa")),
                                                                threshold = 0.9),
                                  P_IM_seqs = filter_with_cdhit(read_fasta(paste0(data_path, "Sequences/P_IM.fa")),
                                                                threshold = 0.9),
                                  N_S_seqs = filter_with_cdhit(read_fasta(paste0(data_path, "Sequences/N_S.fa")),
                                                               threshold = 0.9),
                                  P_S_seqs = filter_with_cdhit(read_fasta(paste0(data_path, "Sequences/P_S.fa")),
                                                               threshold = 0.9),
                                  N_TM_seqs = filter_with_cdhit(read_fasta(paste0(data_path, "Sequences/N_TM.fa")),
                                                                threshold = 0.9),
                                  P_TM_seqs = filter_with_cdhit(read_fasta(paste0(data_path, "Sequences/P_TM.fa")),
                                                                threshold = 0.9),
                                  N_TL_SEC_seqs = filter_with_cdhit(read_fasta(paste0(data_path, "Sequences/N_TL_SEC.fa")),
                                                                    threshold = 0.9),
                                  N_TL_TAT_seqs = filter_with_cdhit(read_fasta(paste0(data_path, "Sequences/N_TL_TAT.fa")),
                                                                    threshold = 0.9),
                                  N_seqs = c(N_OM_seqs, N_IM_seqs, N_S_seqs, N_TM_seqs, N_TL_SEC_seqs, N_TL_TAT_seqs),
                                  P_seqs = c(P_IM_seqs, P_S_seqs, P_TM_seqs),
                                  Plastid_Nuclear_CV = do_cv(list("N" = N_seqs, "P" = P_seqs), 5, "N", 0.001),
                                  Plastid_Nuclear_CV_res_stats = get_cv_res_summary(Plastid_Nuclear_CV, "N"),
                                  Plastid_Nuclear_CV_2 = do_cv(list("N" = N_seqs, "P" = P_seqs), 5, "N", 0.00001),
                                  Plastid_Nuclear_CV_2_res_stats = get_cv_res_summary(Plastid_Nuclear_CV_2)
)

make(analysis_PlastoGram, seed = 12345)