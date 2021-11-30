get_composition_matrix <- function(data_path, groups, length) {
  aa <- c("A", "Q", "W", "E", "R", "T", "Y", "I", "P", "S", "D", "F", "G", "H", "K", "L", "C", "V", "N", "M")
  
  ngram_composition <- lapply(groups, function(ith_group) {
    prots <- read_fasta(paste0(data_path, ith_group, ".fasta"))
    
    prot_matrix <- lapply(1:length(prots), function(ith_prot) {
      as.data.frame(rbind(unname(prots[[ith_prot]][1:length])))
    }) %>% bind_rows() %>% 
      as.matrix()
    
    unigrams <- count_multimers(prot_matrix, k_vector = 1, kmer_alphabet = aa)
    bigrams <- count_multimers(prot_matrix, k_vector = c(2,2), kmer_gaps_list = list(NULL, 1), kmer_alphabet = aa) 
    trigrams <- count_multimers(prot_matrix, k_vector = c(rep(3,4)), kmer_gaps_list = list(c(0,0), c(0,1), c(1,0), c(1,1)), kmer_alphabet = aa)
    
    freqs <- cbind(
      cbind(as.matrix(unigrams)/20, as.matrix(bigrams)/ncol(as.matrix(bigrams))),
      as.matrix(trigrams)/ncol(as.matrix(trigrams)))
    
    mutate(as.data.frame(freqs),
           Dataset = ith_group)
  }) %>% bind_rows()
  
  ngram_composition[is.na(ngram_composition)] <- 0
  
  ngram_composition
}
