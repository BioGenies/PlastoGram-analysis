get_model_variants <- function() {
  base_models <- list(
    v1 = c("Nuclear_model", "Membrane_model", "Nuclear_membrane_model", "Plastid_membrane_model", "Tat_model", "Sec_model"),
    v2 = c("Nuclear_model", "Membrane_model", "Nuclear_membrane_model", "Plastid_membrane_model", "Tat_model", "Sec_model", "Stroma_model"),
    v3 = c("Nuclear_model", "Membrane_model", "Nuclear_membrane_model", "Plastid_membrane_model", "Tat_model", "Sec_model", "OM_Stroma_model"),
    v4 = c("Nuclear_model", "Membrane_model", "N_OM_model", "N_IM_model", "N_TM_model", "Plastid_membrane_model", "Tat_model", "Sec_model"),
    v5 = c("Nuclear_model", "Membrane_model", "N_OM_model", "N_IM_model", "N_TM_model", "Plastid_membrane_model", "Tat_model", "Sec_model", "Stroma_model"),
    v6 = c("Nuclear_model", "Membrane_model", "N_OM_model", "N_IM_model", "N_TM_model", "Plastid_membrane_model", "Tat_model", "Sec_model", "OM_Stroma_model"),
    v7 = c("Nuclear_model", "Membrane_model", "Nuclear_membrane_model", "Plastid_membrane_model", "Tat_model", "Sec_model", "P_stroma_model", "N_stroma_model"),
    v8 = c("Nuclear_model", "Membrane_model", "N_OM_model", "N_IM_model", "N_TM_model", "Plastid_membrane_model", "Tat_model", "Sec_model", "P_stroma_model", "N_stroma_model")
  )
  without_sec <- lapply(base_models, function(i) i[which(i != "Sec_model")]) %>% 
    setNames(paste0(names(base_models), "-Sec"))
  without_tat <- lapply(base_models, function(i) i[which(i != "Tat_model")]) %>% 
    setNames(paste0(names(base_models), "-Tat"))
  
  c(base_models, without_sec, without_tat)
}


generate_all_architectures <- function(model_variants, smote_models, sequence_models, model_dat, filtering_df, output_dir) {
  model_dat <- model_dat[which(!(grepl("SMOTE", model_dat[["Model_name"]]))), ]
  filtering_options <- colnames(filtering_df)[2:ncol(filtering_df)]
  lapply(filtering_options, function(ith_filtering) {
    lapply(seq_along(model_variants), function(ith_variant) {
      n_models <- length(model_variants[[ith_variant]])
      model_names <- model_variants[[ith_variant]]
      n_smote_models <- sum(model_names %in% smote_models)
      pos_smote_models <- which(model_names %in% smote_models)
      res_filtering <- sapply(model_names, function(i) filtering_df[[ith_filtering]][which(filtering_df[["Model_name"]] == i)])
      df <- left_join(
        data.frame(
          Model_name = model_names,
          Test_filtering = sapply(res_filtering, function(i) ifelse(is.na(i), "", i))),
        model_dat)
      if(n_smote_models == 0) {
        architecture <- mutate(df, SMOTE = FALSE)
        write.csv(architecture, file = paste0(output_dir, "Architecture_", names(model_variants)[ith_variant], ".csv"), row.names = FALSE)
      } else {
        lapply(0:n_smote_models, function(i) {
          smote_permutations <- unique(permn(c(rep(TRUE, i), rep(FALSE, n_smote_models-i))))
          lapply(seq_along(smote_permutations), function(ith_perm) {
            smote <- rep(FALSE, n_models)
            smote[pos_smote_models] <- smote_permutations[[ith_perm]]
            architecture <- mutate(df, SMOTE = smote)
            write.csv(architecture, file = paste0(output_dir, "Architecture_", names(model_variants)[ith_variant], "_", i, "-", ith_perm, "_", ith_filtering, ".csv"), row.names = FALSE)
          })
        })
      }
    }) 
  })
}



