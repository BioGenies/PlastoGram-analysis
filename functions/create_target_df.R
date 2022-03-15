create_target_df <- function(annotations_file, sequences, graphpart_res) {
  annot <- read_xlsx(annotations_file) %>% 
    select(Entry, Final_dataset) %>% 
    unique() %>% 
    setNames(c("seq_name", "dataset")) %>% 
    filter(seq_name %in% names(graphpart_res[["AC"]]),
           dataset %in% c("N_OM", "N_IM", "P_IM", "N_S", "P_S", "N_TM", "P_TM", "N_TL_SEC", "N_TL_TAT"))
  
  graphpart_res %>% 
    select(c(AC, dataset, cluster)) %>% 
    setNames(c("seq_name", "dataset", "cluster")) %>% 
    left_join(annot) %>% 
    mutate(Nuclear_target = ifelse(grepl("N_", dataset), TRUE, FALSE),
           Sec_target = ifelse(dataset == "N_TL_SEC", TRUE, FALSE),
           Tat_target = ifelse(dataset == "N_TL_TAT", TRUE, FALSE),
           Membrane_mc_target = case_when(dataset == "N_OM" ~ "OM",
                                          dataset %in% c("N_IM", "P_IM") ~ "IM",
                                          dataset %in% c("N_TM", "P_TM") ~ "TM",
                                          !(dataset %in% c("N_OM", "N_IM", "N_TM", "P_IM", "P_TM")) ~ "other"),
           Stroma_target = ifelse(grepl("_S$", dataset), TRUE, FALSE),
           Membrane_target = ifelse(dataset %in% c("N_OM", "N_IM", "P_IM", "N_TM", "P_TM"), TRUE, FALSE),
           OM_target = ifelse(dataset == "N_OM", TRUE, FALSE),
           IM_target = ifelse(dataset %in% c("N_IM", "P_IM"), TRUE, FALSE),
           TM_target = ifelse(dataset %in% c("N_TM", "P_TM"), TRUE, FALSE))
}
