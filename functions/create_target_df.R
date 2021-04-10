create_target_df <- function(annotations_file, sequences) {
  read_xlsx(annotations_file) %>% 
    select(Entry, Final_dataset) %>% 
    unique() %>% 
    setNames(c("seq_name", "dataset")) %>% 
    filter(seq_name %in% names(sequences),
           dataset %in% c("N_OM", "N_IM", "P_IM", "N_S", "P_S", "N_TM", "P_TM", "N_TL_SEC", "N_TL_TAT")) %>% 
    mutate(NP_target = ifelse(grepl("N_", dataset), 1, 0),
           Sec_target = ifelse(dataset == "N_TL_SEC", 1, 0),
           Tat_target = ifelse(dataset == "N_TL_TAT", 1, 0),
           membrane_target = case_when(dataset == "N_OM" ~ "OM",
                                       dataset %in% c("N_IM", "P_IM") ~ "IM",
                                       dataset %in% c("N_TM", "P_TM") ~ "TM"),
           S_target = ifelse(grepl("_S$", dataset), 1, 0))
}
