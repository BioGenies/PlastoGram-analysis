create_target_df <- function(annotations_file, sequences) {
  read_xlsx(annotations_file) %>% 
    select(Entry, Final_dataset) %>% 
    unique() %>% 
    setNames(c("seq_name", "dataset")) %>% 
    filter(seq_name %in% names(sequences),
           dataset %in% c("N_OM", "N_IM", "P_IM", "N_S", "P_S", "N_TM", "P_TM", "N_TL_SEC", "N_TL_TAT")) %>% 
    mutate(NP_target = ifelse(grepl("N_", dataset), TRUE, FALSE),
           Sec_target = ifelse(dataset == "N_TL_SEC", TRUE, FALSE),
           Tat_target = ifelse(dataset == "N_TL_TAT", TRUE, FALSE),
           membrane_target = case_when(dataset == "N_OM" ~ "OM",
                                       dataset %in% c("N_IM", "P_IM") ~ "IM",
                                       dataset %in% c("N_TM", "P_TM") ~ "TM",
                                       !(dataset %in% c("N_OM", "N_IM", "N_TM", "P_IM", "P_TM")) ~ "other"),
           S_target = ifelse(grepl("_S$", dataset), TRUE, FALSE),
           membrane_all_target = ifelse(dataset %in% c("N_OM", "N_IM", "P_IM", "N_TM", "P_TM"), TRUE, FALSE),
           membrane_OM_target = ifelse(dataset == "N_OM", TRUE, FALSE),
           membrane_IM_target = ifelse(dataset %in% c("N_IM", "P_IM"), TRUE, FALSE),
           membrane_TM_target = ifelse(dataset %in% c("N_TM", "P_TM"), TRUE, FALSE))
}
