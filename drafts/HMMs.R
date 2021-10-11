library(dplyr)

get_max_pos <- function(x) max(which(x == max(x)))

find_rr <- function(x, peptide) {
  
  only_peptide <- x[peptide]
  
  r_pos <- only_peptide == "R"
  
  position_of_the_last_R_in_pair <- c(FALSE, sapply(2L:length(r_pos), function(ith_pos) {
    r_pos[ith_pos] + r_pos[ith_pos - 1] == 2
  }))
  
  putative_RR_end <- get_max_pos(position_of_the_last_R_in_pair)
  
  res <- rep(FALSE, length(only_peptide))
  res[(putative_RR_end - 1):putative_RR_end] <- TRUE
  
  c(res, rep(NA, length(x) - length(only_peptide)))
}

find_cs <- function(peptide) {
  cs_pos <- which.max(cumsum(peptide))
  res <- rep(FALSE, length(peptide))
  res[(cs_pos - 2):(cs_pos + 5)] <- TRUE
  res
}

read_tat <- function(file_lbl) {
  all_lines <- readLines(file_lbl)
  
  n_seqs <- length(all_lines)/3
  
  splitted_seqs <- split(all_lines, unlist(lapply(1L:n_seqs, function(i) rep(i, 3))))
  
  setNames(lapply(splitted_seqs, function(ith_seq) {
    ith_seq_list <- list(pos = 1L:length(strsplit(ith_seq[2], "")[[1]]), aa = strsplit(ith_seq[2], "")[[1]], peptide = strsplit(ith_seq[3], "")[[1]] == "T")
    
    ith_seq_list[["RR"]] <- find_rr(ith_seq_list[["aa"]], ith_seq_list[["peptide"]])
    ith_seq_list[["CS"]] <- find_cs(ith_seq_list[["peptide"]])
    
    ith_seq_list
    
  }), sapply(splitted_seqs, function(ith_seq) ith_seq[[1]]))
}

get_n_region <- function(seqs) {
  lapply(seqs, function(ith_seq) {
    ith_seq[["aa"]][1L:(which.max(ith_seq[["RR"]]) - 1)]
  })
}

get_h_region <- function(seqs) {
  lapply(seqs, function(ith_seq) {
    ith_seq[["aa"]][(which.max(ith_seq[["RR"]]) + 2):(which.max(ith_seq[["CS"]]) - 1)]
  })
}

get_rr_region <- function(seqs) {
  lapply(seqs, function(ith_seq) {
    as.vector(na.omit(ith_seq[["aa"]][ith_seq[["RR"]]]))
  })
}

get_cs_region <- function(seqs) {
  lapply(seqs, function(ith_seq) {
    ith_seq[["aa"]][ith_seq[["CS"]]]
  })
}

get_rest_region <- function(seqs) {
  lapply(seqs, function(ith_seq) {
    ith_seq[["aa"]][(which.max(ith_seq[["CS"]]) + 8):length(ith_seq[["aa"]])]
  })
}

as_aa_factor <- function(x) {
  factor(x, levels = toupper(biogram:::return_elements("prot")))
}

region2t <- function(x) {
  setNames(as.vector(table(as_aa_factor(x))), toupper(biogram:::return_elements("prot")))
}

measure_region <- function(region, max_length = 32) {
  lengths_df <- data.frame(table(region))
  lengths <- structure(lengths_df[["Freq"]], names = as.character(lengths_df[["region"]]))
  res <- rep(0, max_length)
  lengths <- lengths[as.numeric(names(lengths))>0] #removing lengths smaller than 1
  
  start_l <- min(as.numeric(names(lengths)))
  end_l <- max(as.numeric(names(lengths)))
  if(prod(start_l:end_l %in% as.numeric(names(lengths)))){
    max_length_real <- length(lengths) #if all lengths are present in training set
  } else{
    max_length_real <- 1
    sl <- sum(lengths)
    while(sum(lengths[1:max_length_real])/sl <= 0.51) {
      max_length_real <- which.min(start_l:end_l %in% as.numeric(names(lengths))) - 1
      start_l <- start_l + 1  #to assure that max_length_real wouldn't be too small
      max_length_real <- ifelse(max_length_real == 0, length(lengths), max_length_real)
    }
  }
  max_length <- min(max_length, max_length_real)
  
  prop_lengths <- lengths[as.character(1L:max_length)]/sum(lengths[as.character(1L:max_length)], na.rm = TRUE)
  # NA are introduced by lengths not present in vector 1:max_length
  prop_lengths <- prop_lengths[!is.na(prop_lengths)]
  res[as.numeric(names(prop_lengths))] <- prop_lengths
  res
}

tat_seqs <- read_tat("./data/PRED-TAT-data/Train_tat.lbl")

train_tat_hsmm <- function(tat_seqs) {
  
  region_list <- list(n_region = get_n_region(tat_seqs),
                      rr_region = get_rr_region(tat_seqs),
                      h_region = get_h_region(tat_seqs),
                      cs_region = get_cs_region(tat_seqs),
                      rest_region = get_rest_region(tat_seqs))
  
  region_lengths <- as.matrix(data.frame(lapply(region_list, lengths)))

  regional_aa_counts <- lapply(region_list, function(ith_region)
    colSums(do.call(rbind, lapply(ith_region, region2t))))
  
  overall_probs_log <- log(regional_aa_counts[["rest_region"]]/sum(regional_aa_counts[["rest_region"]]))
  
  max_length <- 50
  
  params <- apply(region_lengths[, c("n_region", "rr_region", "h_region", "cs_region")], 2, measure_region, max_length = max_length)
  params[2, "rr_region"] <- 1
  params[8, "cs_region"] <- 1
  params <- cbind(params, rep(1/max_length, max_length))

  pipar <- c(1, 0, 0, 0, 0)
  tpmpar <- matrix(c(0, 1, 0, 0, 0,
                     0, 0, 1, 0, 0,
                     0, 0, 0, 1, 0,
                     0, 0, 0, 0, 1,
                     0, 0, 0, 0, 0), 5, byrow = TRUE)
  
  od <- lapply(regional_aa_counts, function(ith_region) {
    ith_region/sum(ith_region)
  }) %>% 
    do.call(cbind, .) %>% 
    t
  
  list(pipar = pipar, tpmpar = tpmpar, od = od, 
       overall_probs_log = overall_probs_log, params = params)
}

tat_hsmm_model <- train_tat_hsmm(tat_seqs)


predict_hsmm <- function(hsmm_model, prot_seq) {
  prot_seq_numeric <- as.numeric(as_aa_factor(prot_seq)) - 1
  viterbi_res <- signalHsmm:::duration_viterbi(prot_seq_numeric, hsmm_model[["pipar"]], hsmm_model[["tpmpar"]], 
                                               hsmm_model[["od"]], hsmm_model[["params"]])

  last_state_id <- ncol(viterbi_res[["viterbi"]])
  viterbi_path <- viterbi_res[["path"]] + 1
  cs_pos <- ifelse(any(viterbi_path == last_state_id), 
                   max(which(viterbi_path == (last_state_id - 1))) + 1, 
                   length(prot_seq))
  
  prob_hsmm <- viterbi_res[["viterbi"]][cs_pos, viterbi_path[cs_pos]]
  prob_non <- Reduce(function(x, y) x + hsmm_model[["overall_probs_log"]][y + 1], 
                     prot_seq_numeric[1L:cs_pos], 0)
  list(unname(1 - 1/(1 + exp(prob_hsmm - prob_non))),
       data.frame(seq = prot_seq, viterbi_path))
}

lapply(1L:length(tat_seqs), function(i) {
  data.frame(predict_hsmm(tat_hsmm_model, tat_seqs[[i]][["aa"]])[[2]], RR = tat_seqs[[i]][["RR"]])
})

table(sapply(1L:length(tat_seqs), function(i) {
  predict_hsmm(tat_hsmm_model, tat_seqs[[i]][["aa"]])[[1]]
}) > 0.5)
