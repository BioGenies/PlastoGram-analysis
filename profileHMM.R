library(targets)
library(biogram)
library(dplyr)
library(rhmmer)
library(cvTools)

tar_load("N_TL_TAT_seqs")
tar_load("N_TL_SEC_seqs")

all_res <- lapply(1:5, function(i) {
  folded_tat <- cvFolds(length(N_TL_TAT_seqs), K = 5)
  fold_df_tat <- data.frame(seq_name = names(N_TL_TAT_seqs)[folded_tat[["subsets"]]], 
                            fold = folded_tat[["which"]],
                            dataset = "TAT", rep = i)
  folded_sec <- cvFolds(length(N_TL_SEC_seqs), K = 5)
  fold_df_sec <- data.frame(seq_name = names(N_TL_SEC_seqs)[folded_sec[["subsets"]]], 
                            fold = folded_sec[["which"]],
                            dataset = "SEC", rep = i)
  
  lapply(1:5, function(ith_fold) {
    tat_train <- N_TL_TAT_seqs[filter(fold_df_tat, fold != ith_fold)[["seq_name"]]]
    tat_test <- N_TL_TAT_seqs[filter(fold_df_tat, fold == ith_fold)[["seq_name"]]]
    sec_train <- N_TL_SEC_seqs[filter(fold_df_sec, fold != ith_fold)[["seq_name"]]]
    sec_test <- N_TL_SEC_seqs[filter(fold_df_sec, fold == ith_fold)[["seq_name"]]]
    
    write_fasta(tat_train, "tat_train.fa")
    write_fasta(tat_test, "tat_test.fa")
    write_fasta(sec_train, "sec_train.fa")
    write_fasta(sec_test, "sec_test.fa")
    
    system('mafft  --localpair  --maxiterate 16 --reorder "tat_train.fa" > "tat_train_aln.fa"')
    system("hmmbuild tat_model.hmm tat_train_aln.fa")
    system("hmmsearch --tblout tat_res_tat tat_model.hmm tat_test.fa")
    system("hmmsearch --tblout sec_res_tat tat_model.hmm sec_test.fa")
    tat_res <- bind_rows(read_tblout("tat_res_tat"),
                         read_tblout("sec_res_tat"))
    
    system('mafft  --localpair  --maxiterate 16 --reorder "sec_train.fa" > "sec_train_aln.fa"')
    system("hmmbuild sec_model.hmm sec_train_aln.fa")
    system("hmmsearch --tblout tat_res_sec sec_model.hmm tat_test.fa")
    system("hmmsearch --tblout sec_res_sec sec_model.hmm sec_test.fa")
    sec_res <- bind_rows(read_tblout("tat_res_sec"),
                         read_tblout("sec_res_sec"))
    
    results <- bind_rows(fold_df_tat, fold_df_sec) %>% 
      filter(fold == ith_fold) %>% 
      mutate(tat = ifelse(seq_name %in% tat_res[["domain_name"]], TRUE, FALSE),
             sec = ifelse(seq_name %in% sec_res[["domain_name"]], TRUE, FALSE))
    
    file.remove(c("tat_model.hmm", "sec_model.hmm", "tat_res_tat", "sec_res_tat", 
                  "tat_res_sec", "sec_res_sec", "tat_train_aln.fa", "sec_train_aln.fa",
                  "tat_train.fa", "sec_train.fa"))
    
    results
  }) %>% bind_rows()
}) %>% bind_rows()


tat_results <- all_res %>% 
  group_by(rep, fold) %>% 
  summarise(TP = sum(dataset == "TAT" & tat == TRUE),
            TN = sum(dataset == "SEC" & tat == FALSE),
            FP = sum(dataset == "SEC" & tat == TRUE),
            FN = sum(dataset == "TAT" & tat == FALSE),
            acc = measures::ACC(ifelse(dataset == "TAT", TRUE, FALSE), tat))


sec_results <- all_res %>% 
  group_by(rep, fold) %>% 
  summarise(TP = sum(dataset == "SEC" & sec == TRUE),
            TN = sum(dataset == "TAT" & sec == FALSE),
            FP = sum(dataset == "TAT" & sec == TRUE),
            FN = sum(dataset == "SEC" & sec == FALSE),
            acc = measures::ACC(ifelse(dataset == "SEC", TRUE, FALSE), sec))
