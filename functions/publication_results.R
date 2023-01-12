get_mean_performance_of_baseline <- function(baseline_model_cv_res, outfile) {
  lapply(unique(baseline_model_cv_res[["rep"]]), function(ith_rep) {
    lapply(unique(baseline_model_cv_res[["fold"]]), function(ith_fold) {
      dat <- filter(baseline_model_cv_res, rep == ith_rep, fold == ith_fold)
      data.frame(
        model = "baseline",
        rep = ith_rep,
        fold = ith_fold,
        AU1U = multiclass.AU1U(dat[, c("N_E", "N_S", "N_TL_SEC", "N_TL_TAT", "N_TM", "P_IM", "P_S", "P_TM")], dat[["dataset"]]),
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


get_performance_distr_plot <- function(architecture_plot_data, graphpart_architecture_plot_data, baseline_mean_performance, graphpart_baseline_mean_performance, res_path) {
  baseline_plot_data <- mutate(baseline_mean_performance, type = "Holdout") %>% 
    bind_rows(mutate(graphpart_baseline_mean_performance, type = "Partitioning")) %>% 
    pivot_longer(colnames(.)[which(!(colnames(.) %in% c("model", "type")))], values_to = "value", names_to = "measure") %>% 
    filter(measure %in% c("mean_kappa", "mean_AU1U") | (startsWith(measure, "mean") & endsWith(measure, "sens"))) %>% 
    mutate(measure = gsub("_sens", " accuracy", gsub("mean_", "", measure))) 
  all_plot_data <- mutate(architecture_plot_data, type = "Holdout") %>% 
    bind_rows(mutate(graphpart_architecture_plot_data, type = "Partitioning"))  %>% 
    filter(measure %in% c("mean_kappa", "mean_AU1U") | (startsWith(measure, "mean") & endsWith(measure, "sens"))) %>% 
    mutate(measure = gsub("_sens", " accuracy", gsub("mean_", "", measure))) 
  p <- ggplot(all_plot_data, aes(x = measure, y = value)) +
    geom_boxplot() +
    theme_bw() +
    geom_point(data = baseline_plot_data, aes(x = measure, y = value, color = "Baseline"), size = 3) +
    geom_point(data = filter(all_plot_data, model == "Architecture_v71_0-1_No_filtering_RF"), 
               aes(x = measure, y = value, color = "PlastoGram"), size = 3) +
    xlab("Mean performance measure") +
    ylab("Value") +
    scale_color_manual("Model", values = c("PlastoGram" = "#76d872", "Baseline" = "#d87272", Other = "black")) +
    facet_wrap(~type, ncol = 1)
  ggsave("Performance_distribution_architectures+baseline.eps", p, width = 15, height = 10,
         path = res_path)
  p
}

get_parameters_of_performance_distribution_table <- function(envelope_mean_architecture_performance, res_path, outfile) {
  df <- lapply(select(envelope_mean_architecture_performance, 
                      c("mean_kappa", "mean_AU1U", 
                        colnames(envelope_mean_architecture_performance)[startsWith(colnames(envelope_mean_architecture_performance), "mean") 
                                                                         & endsWith(colnames(envelope_mean_architecture_performance), "sens")])), summary) %>% 
    do.call(cbind, .) %>% 
    data.frame() %>% 
    tibble::rownames_to_column() %>% 
    setNames(c("Measure", "Mean kappa", "Mean AU1U", "Mean N_E accuracy", "Mean N_TM accuracy", "Mean N_S accuracy", "Mean N_TL_SEC accuracy",
               "Mean N_TL_TAT accuracy", "Mean P_IM accuracy", "Mean P_TM accuracy", "Mean P_S accuracy")) %>% 
    mutate(Measure = c("Minimum", "1st quartile", "Median", "Mean", "3rd quartile", "Maximum")) %>% 
    t() %>% 
    as.data.frame()
  write.table(df, paste0(res_path, outfile), sep = ",", col.names = FALSE, quote = FALSE)
  df
}



get_best_model_cv_plot <- function(architecture_plot_data, graphpart_architecture_plot_data, res_path) {
  p <- mutate(architecture_plot_data, version = "Holdout") %>% 
    bind_rows(mutate(graphpart_architecture_plot_data, version = "Partitioning")) %>% 
    filter(endsWith(measure, "_sens"),
           model == "Architecture_v71_0-1_No_filtering_RF") %>% 
    select(c("measure", "value", "version")) %>% 
    mutate(type = ifelse(grepl("mean_", measure), "mean", "sd"),
           Class = gsub("_sens", "", gsub("mean_|sd_", "", measure)),
           Origin = ifelse(grepl("N_", Class), "Nuclear-encoded", "Plastid-encoded")) %>% 
    pivot_wider(id_cols = c(Class, Origin, version), names_from = type, values_from = value) %>% 
    ggplot(aes(x = Class, y = mean, fill = version, group = version)) +
    geom_col(position = "dodge") +
    geom_errorbar(aes(ymin = mean-sd, ymax = mean+sd), position = position_dodge(width = 0.9), width = 0.5) +
    theme_bw() +
    #facet_wrap(~version) +
    ylab("Mean class-specific accuracy value in cross-validation") +
    # scale_fill_manual("Dataset", values = c("N_E" = "#7281d8", "N_TM" = "#d87272", "N_S" = "#76d872",
    #                                         "N_TL_SEC" = "#d8c472", "N_TL_TAT" = "#d8a972", 
    #                                         "P_IM" = "#a8b1e8", "P_TM" = "#e8a8a8", "P_S" = "#ade8a8")) +
    scale_fill_manual("Data set version", values = c("Holdout" = "#ade8a8", "Partitioning" = "#76d872"))
  ggsave(paste0(res_path, "Best_model_cv_res.eps"), p, width = 9, height = 5)
}

get_best_model_cv_table <- function(architecture_plot_data) {
  architecture_plot_data %>% 
    filter(startsWith(measure, "mean_") | endsWith(measure, "sens"),
           model == "Architecture_v71_0-1_No_filtering_RF") %>% 
    select(c("measure", "value")) %>% 
    mutate(measure = gsub("_sens", " accuracy", gsub("mean_", "", measure))) %>% 
    setNames(c("Measure", "Value"))# %>% 
  #   write.csv(paste0(res_path, "Best_model_cv_res.csv"), row.names = FALSE)
}


get_independent_dataset_performance_table <- function(PlastoGram_evaluation) {
  res <- PlastoGram_evaluation[["Final_results"]]
  hl_res <- PlastoGram_evaluation[["Higher-order_model_preds"]]
  data.frame(
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
  #  write.csv(df, paste0(res_path, "Independent_dataset_results.csv"), row.names = FALSE)
}

get_cv_and_independent_res_table <- function(architecture_plot_data, PlastoGram_evaluation, res_path, outfile) {
  cv <- get_best_model_cv_table(architecture_plot_data) %>% 
    setNames(c("Measure", "Mean in CV"))
  ind <- get_independent_dataset_performance_table(PlastoGram_evaluation) %>% 
    setNames(c("Measure", "Independent"))
  sd <- cv[11:18,] %>% 
    mutate(Measure = gsub("sd_", "", Measure)) %>% 
    setNames(c("Measure", "SD in CV"))
  res <- left_join(left_join(ind, cv), sd)
  write.csv(res, paste0(res_path, outfile), row.names = FALSE)
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


get_physicochemical_properties_plot <- function(traintest, traintest_data_df, colors) {
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
      theme_bw(base_size = 9) +
      theme(plot.title = element_text(hjust = 0.5),
            legend.position = "none") +
      scale_fill_manual("Dataset", values = colors)
  })
  
  (p_list[[1]] + labs(tag = "B") + theme(plot.tag = element_text(size = '16'))) / p_list[[2]] / p_list[[3]] 
  
}

get_physicochemical_properties_test <- function(traintest, traintest_data_df, data_path) {
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
  combns <- combn(unique(prop_df[["dataset"]]), 2, simplify = FALSE)
  test_res <- lapply(colnames(prop_df)[4:ncol(prop_df)], function(ith_prop) {
    lapply(combns, function(ith_combn) {
      dat <- prop_df[, c("prot", "dataset", ith_prop)] %>% 
        setNames(c("prot", "Dataset", "Property")) %>% 
        filter(Dataset %in% ith_combn)
      data.frame(Property = props[["description"]][which(props[["prop"]] == ith_prop)],
                 Dataset1 = factor(ith_combn[1], levels = unique(prop_df[["dataset"]])),
                 Dataset2 = factor(ith_combn[2], levels = unique(prop_df[["dataset"]])),
                 `p-value` = wilcox.test(x = filter(dat, Dataset == ith_combn[1])[["Property"]],
                                         y = filter(dat, Dataset == ith_combn[2])[["Property"]])[["p.value"]],
                 check.names = FALSE)
    }) %>% bind_rows()
  }) %>% bind_rows() %>% 
    mutate(`Adjusted p-value` = p.adjust(`p-value`, method = "BH"),
           `Is significant` = ifelse(`Adjusted p-value` < 0.05, TRUE, FALSE))
  write.csv(test_res, paste0(data_path, "Physicochemical_properties_statistical_test.csv"), row.names = FALSE)
  test_res
}


get_om_im_model_cv_res_table <- function(traintest_ngram_matrix, holdout_target_dfs_cv, res_path, outfile) {
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
  write.csv(res, paste0(res_path, outfile), row.names = FALSE)
  res
}

get_benchmark_res_table <- function(PlastoGram_evaluation, SChloro_benchmark_res) {
  schloro_df <- SChloro_benchmark_res %>% 
    select(c(V1, V3)) %>% 
    group_by(V1) %>% 
    summarise(pred = paste0(V3, collapse = "; "))
  
  all_benchmark <- left_join(PlastoGram_evaluation[["Final_results"]],
                             schloro_df, by = c('seq_name' = 'V1')) %>% 
    mutate(dataset = gsub("_TAT|_SEC", "", gsub("N_|P_", "", dataset)),
           Localization = gsub("_TAT|_SEC", "", gsub("N_|P_", "", Localization))) %>% 
    mutate(dataset = ifelse(dataset == "IM", "E", dataset),
           Localization = ifelse(Localization == "IM", "E", Localization))
  
  datasets <- list("E" = "inner membrane|outer membrane", "S" = "stroma", 
                   "TM" = "thylakoid membrane", "TL" = "thylakoid lumen")
  
  lapply(names(datasets), function(ith_set) {
    dat <- filter(all_benchmark, dataset == ith_set) %>% 
      mutate(pred2 = ifelse(grepl(datasets[[ith_set]], pred), TRUE, FALSE))
    data.frame(
      Class = ith_set,
      PlastoGram = TPR(ifelse(dat[["dataset"]] == ith_set, TRUE, FALSE),
                       ifelse(dat[["Localization"]] == ith_set, TRUE, FALSE),
                       TRUE),
      SChloro = TPR(ifelse(dat[["dataset"]] == ith_set, TRUE, FALSE),
                    dat[["pred2"]],
                    TRUE)
    )
  }) %>% bind_rows()
}

get_benchmark_res_plot <- function(PlastoGram_evaluation, PlastoGram_evaluation_graphpart, schloro_res, schloro_res_graphpart, res_path) {
  h <- get_benchmark_res_table(PlastoGram_evaluation, schloro_res) %>%
    mutate(type = "Holdout")
  p <- get_benchmark_res_table(PlastoGram_evaluation_graphpart, schloro_res_graphpart) %>%
    mutate(type = "Partitioning")
  plot <- bind_rows(h, p) %>% 
    pivot_longer(2:3, names_to = "Software", values_to = "Class-specific accuracy") %>% 
    ggplot(aes(x = Class, y = `Class-specific accuracy`, group = Software, fill = Software)) +
    geom_col(position = "dodge") +
    facet_wrap(~type) +
    theme_bw() +
    scale_fill_manual("Software", values = c("PlastoGram" = "#76d872", "SChloro" = "#7281d8")) +
    theme(legend.position = "bottom")
  ggsave(paste0(res_path, "Benchmark_results.eps"), plot, width = 7, height = 4)
}

get_final_plastogram_model <- function(PlastoGram_ngram_models, PlastoGram_higher_level_model,
                                       PlastoGram_OM_IM_model, PlastoGram_informative_ngrams) {
  
  PlastoGram_model <- list("ngram_models" = PlastoGram_ngram_models,
                           "RF_model" = PlastoGram_higher_level_model,
                           "OM_IM_model" = PlastoGram_OM_IM_model,
                           "imp_ngrams" = PlastoGram_informative_ngrams)
  class(PlastoGram_model) <- "plastogram_model"
  PlastoGram_model
}

calculate_aa_comp_peptides <- function(seqs) {
  aa_comp_peptides <- lapply(names(seqs), function(i) {
    ds <- seqs[[i]]
    lapply(names(ds), function(ith_prot) {
      data.frame(table(factor(ds[[ith_prot]], levels = c("A", "C", "D", "E", "F", "G", "H",
                                                         "I", "K", "L", "M", "N", "P", "Q",
                                                         "R", "S", "T", "V", "W", "Y")))/length(ds[[ith_prot]])) %>%
        setNames(c("Amino acid", "Frequency"))  %>%
        mutate(dataset = i,
               prot_name = ith_prot)
    }) %>% bind_rows()
  }) %>% bind_rows() %>% 
    mutate(dataset = ifelse(dataset %in% c("N_TL_TAT", "N_TL_SEC"), "N_TL", dataset)) %>% 
    pivot_wider(names_from = "Amino acid", values_from = "Frequency") 
}

plot_aa_comp_pca <- function(aa_comp_peptides) {
  pca_res <- prcomp(select(aa_comp_peptides, -c("dataset", "prot_name")), center = TRUE)
  plot_dat <- mutate(aa_comp_peptides, 
                     PC1 = pca_res[["x"]][, "PC1"],
                     PC2 = pca_res[["x"]][, "PC2"]) 
  
  ggbiplot(pca_res, groups = aa_comp_peptides[["dataset"]], ellipse = TRUE, alpha = 0, var.axes = FALSE) +
    theme_bw() +
    geom_point(size = 0.25, aes(color = aa_comp_peptides[["dataset"]])) +
    scale_color_manual("Dataset", values =  c("N_OM" = "#8210c9", "N_IM" = "#344feb", "P_IM" = "#69cdff", 
                                              "N_S" = "#118a0c", "P_S" = "#9fff9c", "N_TM" = "#911a1a", 
                                              "P_TM" = "#ff7d7d", "N_TL" = "#ccbf10")) 
}

get_pca_and_props_plot <- function(seqs, traintest, traintest_data_df, res_path) {
  props_plot <- get_physicochemical_properties_plot(traintest, traintest_data_df, 
                                                    colors = c("N_OM" = "#b172d8", "N_IM" = "#7281d8", "N_TM" = "#d87272", "N_S" = "#76d872"))
  pca_plot <- plot_aa_comp_pca(calculate_aa_comp_peptides(seqs)) + 
    theme(legend.position = "bottom")
  res <- ggarrange(pca_plot + labs(tag = "A") + theme(plot.tag = element_text(size = '16')), props_plot,  widths = c(6, 3.5)) 
  ggsave(paste0(res_path, "PCA_props_plot.eps"), res, width = 11, height = 8)
}


# This is a modified version of ggbiplot function from ggbiplot package
# https://github.com/vqv/ggbiplot/blob/master/R/ggbiplot.r
# Copyright 2011 Vincent Q. Vu.
ggbiplot <- function(pcobj, choices = 1:2, scale = 1, pc.biplot = TRUE, 
                     obs.scale = 1 - scale, var.scale = scale, 
                     groups = NULL, ellipse = FALSE, ellipse.prob = 0.68, 
                     labels = NULL, labels.size = 3, alpha = 1, 
                     var.axes = TRUE, 
                     circle = FALSE, circle.prob = 0.69, 
                     varname.size = 3, varname.adjust = 1.5, 
                     varname.abbrev = FALSE, ...)
{
  library(ggplot2)
  library(scales)
  library(grid)
  
  stopifnot(length(choices) == 2)
  
  # Recover the SVD
  if(inherits(pcobj, 'prcomp')){
    nobs.factor <- sqrt(nrow(pcobj$x) - 1)
    d <- pcobj$sdev
    u <- sweep(pcobj$x, 2, 1 / (d * nobs.factor), FUN = '*')
    v <- pcobj$rotation
  } else if(inherits(pcobj, 'princomp')) {
    nobs.factor <- sqrt(pcobj$n.obs)
    d <- pcobj$sdev
    u <- sweep(pcobj$scores, 2, 1 / (d * nobs.factor), FUN = '*')
    v <- pcobj$loadings
  } else if(inherits(pcobj, 'PCA')) {
    nobs.factor <- sqrt(nrow(pcobj$call$X))
    d <- unlist(sqrt(pcobj$eig)[1])
    u <- sweep(pcobj$ind$coord, 2, 1 / (d * nobs.factor), FUN = '*')
    v <- sweep(pcobj$var$coord,2,sqrt(pcobj$eig[1:ncol(pcobj$var$coord),1]),FUN="/")
  } else if(inherits(pcobj, "lda")) {
    nobs.factor <- sqrt(pcobj$N)
    d <- pcobj$svd
    u <- predict(pcobj)$x/nobs.factor
    v <- pcobj$scaling
    d.total <- sum(d^2)
  } else {
    stop('Expected a object of class prcomp, princomp, PCA, or lda')
  }
  
  # Scores
  choices <- pmin(choices, ncol(u))
  df.u <- as.data.frame(sweep(u[,choices], 2, d[choices]^obs.scale, FUN='*'))
  
  # Directions
  v <- sweep(v, 2, d^var.scale, FUN='*')
  df.v <- as.data.frame(v[, choices])
  
  names(df.u) <- c('xvar', 'yvar')
  names(df.v) <- names(df.u)
  
  if(pc.biplot) {
    df.u <- df.u * nobs.factor
  }
  
  # Scale the radius of the correlation circle so that it corresponds to 
  # a data ellipse for the standardized PC scores
  r <- sqrt(qchisq(circle.prob, df = 2)) * prod(colMeans(df.u^2))^(1/4)
  
  # Scale directions
  v.scale <- rowSums(v^2)
  df.v <- r * df.v / sqrt(max(v.scale))
  
  # Change the labels for the axes
  if(obs.scale == 0) {
    u.axis.labs <- paste('standardized PC', choices, sep='')
  } else {
    u.axis.labs <- paste('PC', choices, sep='')
  }
  
  # Append the proportion of explained variance to the axis labels
  u.axis.labs <- paste(u.axis.labs, 
                       sprintf('(%0.1f%% explained var.)', 
                               100 * pcobj$sdev[choices]^2/sum(pcobj$sdev^2)))
  
  # Score Labels
  if(!is.null(labels)) {
    df.u$labels <- labels
  }
  
  # Grouping variable
  if(!is.null(groups)) {
    df.u$groups <- groups
  }
  
  # Variable Names
  if(varname.abbrev) {
    df.v$varname <- abbreviate(rownames(v))
  } else {
    df.v$varname <- rownames(v)
  }
  
  # Variables for text label placement
  df.v$angle <- with(df.v, (180/pi) * atan(yvar / xvar))
  df.v$hjust = with(df.v, (1 - varname.adjust * sign(xvar)) / 2)
  
  # Base plot
  g <- ggplot(data = df.u, aes(x = xvar, y = yvar)) + 
    xlab(u.axis.labs[1]) + ylab(u.axis.labs[2]) + coord_equal()
  
  if(var.axes) {
    # Draw circle
    if(circle) 
    {
      theta <- c(seq(-pi, pi, length = 50), seq(pi, -pi, length = 50))
      circle <- data.frame(xvar = r * cos(theta), yvar = r * sin(theta))
      g <- g + geom_path(data = circle, color = muted('white'), 
                         size = 1/2, alpha = 1/3)
    }
    
    # Draw directions
    g <- g +
      geom_segment(data = df.v,
                   aes(x = 0, y = 0, xend = xvar, yend = yvar),
                   arrow = arrow(length = unit(1/2, 'picas')), 
                   color = muted('red'))
  }
  
  # Draw either labels or points
  if(!is.null(df.u$labels)) {
    if(!is.null(df.u$groups)) {
      g <- g + geom_text(aes(label = labels, color = groups), 
                         size = labels.size)
    } else {
      g <- g + geom_text(aes(label = labels), size = labels.size)      
    }
  } else {
    if(!is.null(df.u$groups)) {
      g <- g + geom_point(aes(color = groups), alpha = alpha)
    } else {
      g <- g + geom_point(alpha = alpha)      
    }
  }
  
  # Overlay a concentration ellipse if there are groups
  if(!is.null(df.u$groups) && ellipse) {
    theta <- c(seq(-pi, pi, length = 50), seq(pi, -pi, length = 50))
    circle <- cbind(cos(theta), sin(theta))
    
    ell <- plyr::ddply(df.u, 'groups', function(x) {
      if(nrow(x) <= 2) {
        return(NULL)
      }
      sigma <- var(cbind(x$xvar, x$yvar))
      mu <- c(mean(x$xvar), mean(x$yvar))
      ed <- sqrt(qchisq(ellipse.prob, df = 2))
      data.frame(sweep(circle %*% chol(sigma) * ed, 2, mu, FUN = '+'), 
                 groups = x$groups[1])
    })
    names(ell)[1:2] <- c('xvar', 'yvar')
    g <- g + geom_path(data = ell, aes(color = groups, group = groups), size = 1)
  }
  
  # Label the variable axes
  if(var.axes) {
    g <- g + 
      geom_text(data = df.v, 
                aes(label = varname, x = xvar, y = yvar, 
                    angle = angle, hjust = hjust), 
                color = 'darkred', size = varname.size)
  }
  
  
  return(g)
}


get_mean_performance_of_om_im_models <- function(performance_results, data_path, outfile) {
  all_measures <- performance_results %>% 
    group_by(Replication) %>% 
    summarise(mean_rep_kappa = mean(Kappa),
              mean_rep_AUC = mean(AUC),
              mean_rep_N_OM_sens = mean(Sensitivity),
              mean_rep_N_IM_sens = mean(1-Sensitivity),
              sd_rep_kappa = sd(Kappa),
              sd_rep_AUC = sd(AUC),
              sd_rep_N_OM_sens = sd(Sensitivity),
              sd_rep_N_IM_sens = sd(1-Sensitivity)) %>% 
    ungroup() %>% 
    summarise(mean_kappa = mean(mean_rep_kappa),
              mean_AUC = mean(mean_rep_AUC),
              mean_N_OM_sens = mean(mean_rep_N_OM_sens),
              mean_N_IM_sens = mean(mean_rep_N_IM_sens),
              sd_kappa = mean(sd_rep_kappa),
              sd_AUC = mean(sd_rep_AUC),
              sd_N_OM_sens = mean(sd_rep_N_OM_sens),
              sd_N_IM_sens = mean(sd_rep_N_IM_sens)) %>% 
    pivot_longer(colnames(.), values_to = "value", names_to = "measure") %>% 
    # filter(measure %in% c("mean_kappa", "mean_AU1U") | (startsWith(measure, "mean") & endsWith(measure, "sens"))) %>% 
    mutate(measure = gsub("_sens", " accuracy", gsub("mean_", "", measure)))
  
  summ <- left_join(filter(all_measures, !grepl("sd_", measure)), 
                    mutate(filter(all_measures, grepl("sd_", measure)), measure = gsub("sd_", "", measure)), by = "measure") %>% 
    setNames(c("Measure", "Mean", "SD"))
    
    write.csv(summ, paste0(data_path, outfile), row.names = FALSE)
  summ
}

get_cv_data_for_statistical_test <- function(architectures_performance, best_architecture_name, version_name) {
  architectures_performance %>% 
    filter(model == best_architecture_name) %>% 
    select(-c(model, rep, fold, weighted_kappa)) %>% 
    pivot_longer(colnames(.), names_to = "Measure", values_to = "Value") %>% 
    mutate(Version = version_name)
}

get_cv_statistical_test <- function(holdout_architectures_performance, holdout_best_architecture_name, holdout_version_name,
                                    partitioning_architectures_performance, partitioning_best_architecture_name, partitioning_version_name,
                                    res_path) {
  test_dat <- bind_rows(
    get_cv_data_for_statistical_test(holdout_architectures_performance, 
                                     holdout_best_architecture_name, 
                                     "Holdout"),
    get_cv_data_for_statistical_test(partitioning_architectures_performance, 
                                     partitioning_best_architecture_name, 
                                     "Partitioning"))
  
  test_res <- lapply(unique(test_dat[["Measure"]]), function(ith_measure) {
    data.frame(Measure = ith_measure,
               pvalue = kruskal.test(filter(test_dat, Measure == ith_measure)[["Value"]],
                                     filter(test_dat, Measure == ith_measure)[["Version"]])[["p.value"]])
  }) %>% bind_rows() %>% 
    mutate(adjusted_pvalue = p.adjust(pvalue, method = "BH"))
  
  write.csv(test_res, paste0(res_path, "CV_results_statistical_test.csv"), row.names = FALSE)
  test_res
}

get_feature_size_table <- function(PlastoGram_ngram_models, PlastoGram_ngram_models_graphpart, res_path) {
  feature_size_df <- lapply(names(PlastoGram_ngram_models), function(ith_model) {
    data.frame(Model = ith_model,
               Holdout = PlastoGram_ngram_models[[ith_model]][["num.independent.variables"]],
               Partitioning = PlastoGram_ngram_models_graphpart[[ith_model]][["num.independent.variables"]])
  }) %>% bind_rows()
  write.csv(feature_size_df, paste0(res_path, "Feature_size.csv"), row.names = FALSE)
  feature_size_df
}

change_model_names <- function(architecture_df) {
  mutate(architecture_df, Model_name = case_when(Model_name == "Nuclear_membrane_model" ~ "N_E_vs_N_TM_model",
                                                 Model_name == "E_stroma_model" ~ "N_E_vs_N_S_model",
                                                 TRUE ~ Model_name))
}
