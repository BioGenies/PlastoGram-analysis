get_mean_performance_of_baseline <- function(baseline_model_cv_res, outfile) {
  lapply(unique(baseline_model_cv_res[["rep"]]), function(ith_rep) {
    lapply(unique(baseline_model_cv_res[["fold"]]), function(ith_fold) {
      dat <- filter(baseline_model_cv_res, rep == ith_rep, fold == ith_fold)
      data.frame(
        model = "baseline",
        rep = ith_rep,
        fold = ith_fold,
        AU1U = multiclass.AU1U(dat[, 4:(ncol(dat)-2)], dat[["dataset"]]),
        kappa = KAPPA(dat[["dataset"]], dat[["Localization"]]),
        N_E_sensitivity = TPR(ifelse(dat[["dataset"]] == "N_E", TRUE, FALSE),
                              ifelse(dat[["Localization"]] == "N_E", TRUE, FALSE), TRUE),
        N_TM_sensitivity = TPR(ifelse(dat[["dataset"]] == "N_TM", TRUE, FALSE),
                               ifelse(dat[["Localization"]] == "N_TM", TRUE, FALSE), TRUE),
        N_S_sensitivity = TPR(ifelse(dat[["dataset"]] == "N_S", TRUE, FALSE),
                              ifelse(dat[["Localization"]] == "N_S", TRUE, FALSE), TRUE),
        N_TL_SEC_sensitivity = TPR(ifelse(dat[["dataset"]] == "N_TL_SEC", TRUE, FALSE),
                                   ifelse(dat[["Localization"]] == "N_TL_SEC", TRUE, FALSE), TRUE),
        N_TL_TAT_sensitivity = TPR(ifelse(dat[["dataset"]] == "N_TL_TAT", TRUE, FALSE),
                                   ifelse(dat[["Localization"]] == "N_TL_TAT", TRUE, FALSE), TRUE),
        P_IM_sensitivity = TPR(ifelse(dat[["dataset"]] == "P_IM", TRUE, FALSE),
                               ifelse(dat[["Localization"]] == "P_IM", TRUE, FALSE), TRUE),
        P_TM_sensitivity = TPR(ifelse(dat[["dataset"]] == "P_TM", TRUE, FALSE),
                               ifelse(dat[["Localization"]] == "P_TM", TRUE, FALSE), TRUE),
        P_S_sensitivity = TPR(ifelse(dat[["dataset"]] == "P_S", TRUE, FALSE),
                              ifelse(dat[["Localization"]] == "P_S", TRUE, FALSE), TRUE)
      )
    }) %>% bind_rows()
  }) %>% bind_rows() %>% 
    get_mean_performance_of_envelope_architectures(outfile)
}


get_architecture_plot_data <- function(mean_architecture_performance) {
  variants <- c(rep(c("NM1", "NM1_S", "NM2", "NM2_S", "NM1_PSNS", "NM2_PSNS", "NM1_ES", "NM2_ES"),
                    41))
  added <- c("TL", "NTL", "TL_NMall", "NTL_NMall", "TL_E", "NTL_E", "TL_TMall", "NTL_TMall", "NMall", "E", 
             "TMall", "NMall_E", "NMall_TMall", "E_TMall", "TL_NMall_E", "TL_NMall_TMall", "TL_TMall_E", 
             "NTL_NMall_E", "NTL_NMall_TMall", "NTL_TMall_E", "TL_NMall_E_TMall", "NTL_NMall_E_TMall", 
             "TL_NPE", "NTL_NPE", "NMall_NPE", "E_NPE", "TMall_NPE", "TL_NPE_NMall", "NTL_NPE_NMall", 
             "TL_NPE_E", "NTL_NPE_E", "TL_NPE_TMall", "NTL_NPE_TMall", "TL_E_NPE_NMall", "NTL_E_NPE_NMall",
             "TL_E_NPE_TMall", "NTL_E_NPE_TMall", "TL_TMall_NPE_NMall", "NTL_TMall_NPE_NMall", 
             "TL_TMall_NPE_NMall_E", "NTL_TMall_NPE_NMall_E")
  variants <- sapply(0:(length(variants)-1), function(i) paste0(variants[(1*i+1):(1*i+8)], "_", added[i+1]), 
                     simplify = FALSE) %>% 
    unlist()
  
  var_names_df <- data.frame(
    variant = paste0("v", 1:length(variants)),
    name = variants
  )
  arch_df <- mutate(mean_architecture_performance,
                    arch_variant = sapply(mean_architecture_performance[["model"]], function(i) gsub("-Sec|-Tat", "", strsplit(i, "_")[[1]][2])),
                    hierarchical = grepl(pattern = "Filtering", x = model, ignore.case = FALSE),
                    algorithm = sapply(mean_architecture_performance[["model"]], function(i) last(strsplit(i, "_")[[1]])),
                    sec = !grepl(pattern = "-Sec", x = model, ignore.case = FALSE),
                    tat = !grepl(pattern = "-Tat", x = model, ignore.case = FALSE),
                    smote = !grepl(pattern = "0-1", x = model),
                    smote_variant = sapply(mean_architecture_performance[["model"]], function(i) strsplit(i, "_")[[1]][3])) %>% 
    left_join(var_names_df, by = c("arch_variant" = "variant")) %>% 
    mutate(NM_type = sapply(.[["name"]], function(i) gsub("NM", "", strsplit(i, "_")[[1]][1])),
           TL_model = case_when(grepl("NTL", name) ~ "NTL",
                                grepl("_TL", name) ~ "TL"),
           E_model = ifelse(grepl("_E", name), TRUE, FALSE),
           NPE_model = ifelse(grepl("NPE", name), TRUE, FALSE),
           TMall = ifelse(grepl("TMall", name), TRUE, FALSE),
           NMall_model = ifelse(grepl("_NMall", name), TRUE, FALSE),
           stroma_model = case_when(grepl("ES", name) ~ "NE_Stroma",
                                    grepl("PSNS", name) ~ "NS+PS",
                                    grepl("_S", name) ~ "Stroma")) 
  pivot_longer(arch_df, cols = !c(model, arch_variant, name, NM_type, TL_model, E_model, NPE_model, TMall,
                                  NMall_model, hierarchical, algorithm, sec, tat, stroma_model, smote,
                                  smote_variant), names_to = "measure") 
  
}


get_performance_distr_plot <- function(architecture_plot_data, baseline_mean_performance, res_path) {
  baseline_plot_data <- pivot_longer(baseline_mean_performance, 2:ncol(baseline_mean_performance), 
                                     values_to = "value", names_to = "measure")
  p <- architecture_plot_data %>% 
    filter(measure %in% c("mean_kappa", "mean_AU1U")) %>% 
    mutate(measure = ifelse(measure == "mean_kappa", "Mean Kappa", "Mean AU1U")) %>% 
    ggplot(aes(x = measure, y = value)) +
    geom_boxplot() +
    theme_bw() +
    geom_point(data = mutate(filter(baseline_plot_data, measure %in% c("mean_kappa", "mean_AU1U")), 
                             measure = ifelse(measure == "mean_kappa", "Mean Kappa", "Mean AU1U")), 
               aes(x = measure, y = value), color = "red", size = 3) +
    xlab("Performance measure") +
    ylab("Value")
  ggsave("Performance_distribution_architectures+baseline.eps", p, width = 6, height = 8,
         path = res_path)
  p
}

get_parameters_of_performance_distribution_table <- function(envelope_mean_architecture_performance, res_path) {
  df <- lapply(select(envelope_mean_architecture_performance, c("mean_kappa", "mean_AU1U")), summary) %>% 
    do.call(cbind, .) %>% 
    data.frame() %>% 
    tibble::rownames_to_column() %>% 
    setNames(c("Measure", "Mean kappa", "Mean AU1U")) %>% 
    mutate(Measure = c("Minimum", "1st quartile", "Median", "Mean", "3rd quartile", "Maximum"))
  write.csv(df, paste0(res_path, "Performance_distribution_parameters.csv"), row.names = FALSE)
  df
}



get_best_model_cv_plot <- function(architecture_plot_data, res_path) {
  p <- architecture_plot_data %>% 
    filter(endsWith(measure, "_sens"),
           model == "Architecture_v71_0-1_No_filtering_RF") %>% 
    select(c("measure", "value")) %>% 
    mutate(type = ifelse(grepl("mean_", measure), "mean", "sd"),
           Class = gsub("_sens", "", gsub("mean_|sd_", "", measure)),
           Origin = ifelse(grepl("N_", Class), "Nuclear-encoded", "Plastid-encoded")) %>% 
    pivot_wider(id_cols = c(Class, Origin), names_from = type, values_from = value) %>% 
    ggplot(aes(x = Class, y = mean, fill = Class)) +
    geom_col() +
    geom_errorbar(aes(ymin = mean-sd, ymax = mean+sd), width = 0.2) +
    theme_bw() +
    ylab("Mean sensitivity value in cross-validation") +
    scale_fill_manual("Dataset", values = c("N_E" = "#7281d8", "N_TM" = "#d87272", "N_S" = "#76d872",
                                            "N_TL_SEC" = "#d8c472", "N_TL_TAT" = "#d8a972", 
                                            "P_IM" = "#a8b1e8", "P_TM" = "#e8a8a8", "P_S" = "#ade8a8")) +
    theme(legend.position = "none")
  ggsave(paste0(res_path, "Best_model_cv_res.eps"), width = 9, height = 6)
}

get_best_model_cv_table <- function(architecture_plot_data, res_path) {
  architecture_plot_data %>% 
    filter(startsWith(measure, "mean_"),
           model == "Architecture_v71_0-1_No_filtering_RF") %>% 
    select(c("measure", "value")) %>% 
    mutate(measure = gsub("_sens", " sensitivity", gsub("mean_", "", measure))) %>% 
    setNames(c("Measure", "Value")) %>% 
    write.csv(paste0(res_path, "Best_model_cv_res.csv"), row.names = FALSE)
}


get_independent_dataset_performance_table <- function(PlastoGram_evaluation, res_path) {
  res <- PlastoGram_evaluation[["Final_results"]]
  hl_res <- PlastoGram_evaluation[["Higher-order_model_preds"]]
  df <- data.frame(
    kappa = KAPPA(res[["dataset"]], res[["Localization"]]),
    AU1U = multiclass.AU1U(hl_res[, 3:10], hl_res[["dataset"]]),
    `N_E accuracy` = ACC(filter(res, dataset == "N_E")[["dataset"]], filter(res, dataset == "N_E")[["Localization"]]),
    `N_TM accuracy` = ACC(filter(res, dataset == "N_TM")[["dataset"]], filter(res, dataset == "N_TM")[["Localization"]]),
    `N_S accuracy` = ACC(filter(res, dataset == "N_S")[["dataset"]], filter(res, dataset == "N_S")[["Localization"]]),
    `N_TL_SEC accuracy` = ACC(filter(res, dataset == "N_TL_SEC")[["dataset"]], filter(res, dataset == "N_TL_SEC")[["Localization"]]),
    `N_TL_TAT accuracy` = ACC(filter(res, dataset == "N_TL_TAT")[["dataset"]], filter(res, dataset == "N_TL_TAT")[["Localization"]]),
    `P_IM accuracy` = ACC(filter(res, dataset == "P_IM")[["dataset"]], filter(res, dataset == "P_IM")[["Localization"]]),
    `P_TM accuracy` = ACC(filter(res, dataset == "P_TM")[["dataset"]], filter(res, dataset == "P_TM")[["Localization"]]),
    `P_S accuracy` = ACC(filter(res, dataset == "P_S")[["dataset"]], filter(res, dataset == "P_S")[["Localization"]]),
    check.names = FALSE
  ) %>% pivot_longer(1:ncol(.), names_to = "Measure", values_to = "Value") 
  write.csv(df, paste0(res_path, "Independent_dataset_results.csv"), row.names = FALSE)
  df
}

encode_seq <- function(x, property) {
  sapply(x, function(ith_seq) {
    mean(aaprop[property, tolower(ith_seq)])
  })
}

calculate_properties <- function(datasets) {
  lapply(names(datasets), function(ith_dataset) {
    ds_neg <- datasets[[ith_dataset]]
    data.frame(prot = names(ds_neg),
               dataset = ith_dataset,
               len = lengths(ds_neg),
               KLEP840101 = encode_seq(ds_neg, "KLEP840101"),
               KYTJ820101 = encode_seq(ds_neg, "KYTJ820101"),
               ZIMJ680103 = encode_seq(ds_neg, "ZIMJ680103"))
  }) %>% bind_rows()
}


get_physicochemical_properties_plot <- function(traintest, traintest_data_df, colors, res_path) {
  props <- data.frame(prop = c("KLEP840101", "KYTJ820101", "ZIMJ680103"),
                      description = c("Net charge (Klein et al., 1984)",
                                      "Hydropathy index (Kyte-Doolittle, 1982)", 
                                      "Polarity (Zimmerman et al., 1968)"))
  
  OM_t <- traintest[which(names(traintest) %in% filter(traintest_data_df, dataset == "N_OM")[["seq_name"]])]
  N_IM_t <- traintest[which(names(traintest) %in% filter(traintest_data_df, dataset == "N_IM")[["seq_name"]])]
  N_TM_t <- traintest[which(names(traintest) %in% filter(traintest_data_df, dataset == "N_TM")[["seq_name"]])]
  N_S_t <- traintest[which(names(traintest) %in% filter(traintest_data_df, dataset == "N_S")[["seq_name"]])]
  
  ds <- list("N_OM" = OM_t, "N_IM" = N_IM_t, "N_TM" = N_TM_t, "N_S" = N_S_t)
  prop_df <- calculate_properties(ds)
  
  p_list <- lapply(colnames(prop_df)[4:ncol(prop_df)], function(ith_prop) {
    plot_dat <- prop_df[, c("prot", "dataset", ith_prop)] %>% 
      setNames(c("prot", "Dataset", "Property"))
    ggplot(plot_dat, aes(x = Dataset, y = Property, group = Dataset, fill = Dataset)) +
      geom_violin() +
      ggtitle(props[["description"]][which(props[["prop"]] == ith_prop)]) + 
      xlab(NULL) +
      ylab(NULL) +
      theme_bw(base_size = 8) +
      theme(plot.title = element_text(hjust = 0.5),
            legend.position = "none") +
      scale_fill_manual("Dataset", values = colors)
  })
  
  full_plot <- p_list[[1]] + p_list[[2]] + p_list[[3]]
  
  ggsave(paste0(res_path, "OM_stroma_properties.eps"), full_plot, width = 9, height = 2.5)
}
  

get_om_im_model_cv_res_table <- function(traintest_ngram_matrix, holdout_target_dfs_cv, res_path) {
  res <- lapply(1:length(holdout_target_dfs_cv), function(ith_rep) {
    dat <- traintest_ngram_matrix[which(holdout_target_dfs_cv[[ith_rep]][["Membrane_mc_target"]] %in% c("OM", "IM")),]
    dat_df <- filter(holdout_target_dfs_cv[[ith_rep]], Membrane_mc_target %in% c("OM", "IM"))
    cv_res <- do_cv(dat, dat_df, "OM_target", 5, 0.1)
    perf <- lapply(unique(cv_res[["fold"]]), function(ith_fold) {
      dat <- filter(cv_res, fold == ith_fold)
      data.frame(
        Replication = ith_rep,
        Fold = ith_fold,
        AUC = AUC(dat[["prob"]], dat[["target"]], FALSE, TRUE),
        Kappa = KAPPA(dat[["target"]], dat[["pred"]]),
        Sensitivity = TPR(dat[["target"]], dat[["pred"]], TRUE),
        TP = TP(dat[["target"]], dat[["pred"]], TRUE),
        TN = TN(dat[["target"]], dat[["pred"]], FALSE),
        FP = FP(dat[["target"]], dat[["pred"]], TRUE),
        FN = FN(dat[["target"]], dat[["pred"]], FALSE)
      )
    }) %>% bind_rows()
  }) %>% bind_rows()
  write.csv(res, paste0(res_path, "OM_IM_model_cv_res.csv"), row.names = FALSE)
  res
}
