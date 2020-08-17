diff_abundance <- function(data, sample, grouping, condition, intensity, mean = NULL, sd = NULL, n_samples = NULL, ref_condition, filter_NA_missingness = TRUE, method = c("t-test", "t-test_mean_sd", "moderated_t-test", "proDA"), p_adj_method = "BH"){
  method <- match.arg(method)
  
  if(method == "t-test"){
    message("[1/3] Asign missingenss information ... ", appendLF = FALSE)
    
    t_test_missingness <- assign_missingness(
      data,
      sample = {{sample}},
      condition = {{condition}},
      grouping = {{grouping}},
      intensity = {{intensity}},
      ref_condition = ref_condition
    )
    
    message("DONE", appendLF = TRUE)
    message("[2/3] Create input for t-tests ... ", appendLF = FALSE)
    
    t_test_missingness_obs <- t_test_missingness %>%
      drop_na(missingness, {{intensity}}) %>%
      group_by(comparison, {{grouping}}) %>%
      mutate(n_obs = n()) %>%
      ungroup() %>%
      distinct({{grouping}}, comparison, missingness, n_obs)
    
    t_test_input <- t_test_missingness %>%
      dplyr::group_by(comparison, {{grouping}}, {{condition}}) %>%
      dplyr::summarize(intensity = list({{intensity}}), .groups = "drop") %>%
      dplyr::mutate(type = ifelse({{condition}} == "control_AE", "control", "treated")) %>%
      dplyr::select(-{{condition}}) %>%
      tidyr::pivot_wider(names_from = type, values_from = intensity)
    
    message("DONE", appendLF = TRUE)
    message("[3/3] Calculate t-tests ... ", appendLF = FALSE)
    
    t_test_result <- t_test_input %>%
      dplyr::mutate(pval = map2(
        .x = treated,
        .y = control,
        .f = function(.x, .y) {
          tryCatch({
            t.test(.x, .y)
          },
          error = function(error) {
            NA
          })
        }
      )) %>%
      dplyr::mutate(std_error = map_dbl(
        .x = pval,
        .f = ~ tryCatch({.x$stderr
        },
        error = function(error) {
          NA
        })
      )) %>%
      dplyr::mutate(pval = map_dbl(
        .x = pval,
        .f = ~ tryCatch({.x$p.value
        },
        error = function(error) {
          NA
        })
      )) %>%
      dplyr::mutate(diff = map2_dbl(
        .x = treated,
        .y = control,
        .f = function(.x, .y) {
            mean(.x, na.rm = TRUE) - mean(.y, na.rm = TRUE)
        }
      )) %>%
      mutate(diff = ifelse(diff == "NaN", NA, diff)) %>%
      dplyr::mutate(adj_pval = p.adjust(pval, method = p_adj_method)) %>%
      dplyr::select(-c(control, treated)) %>%
      left_join(t_test_missingness_obs, by = c(rlang::as_name(rlang::enquo(grouping)), "comparison"))
    
    message("DONE", appendLF = TRUE)
    
    return(t_test_result)
  }
  
  if(method == "t-test_mean_sd") {
    conditions_no_ref <- unique(pull(data, {{condition}}))[!unique(pull(data, {{condition}})) %in% ref_condition]
     t_test_mean_sd_result <- data %>%
       dplyr::mutate(comparison = ifelse(
         {{condition}} %in% conditions_no_ref,
         paste0({{condition}}, "_vs_", ref_condition),
         list(paste0(conditions_no_ref, "_vs_",  ref_condition))
       )) %>%
       unnest(comparison) %>%
       rename(mean = {{mean}}, sd = {{sd}}, n = {{n_samples}}) %>%
       split(.$comparison) %>%
       map_df( ~ .x %>%
                dplyr::mutate({{condition}} := ifelse({{condition}} == "control_AE", "control", "treated")) %>%
                pivot_wider(names_from = {{condition}}, values_from = c(.data$mean, .data$sd, .data$n)) %>%
                mutate(ttest_protti(.data$mean_control, .data$mean_treated, .data$sd_control, .data$sd_treated, .data$n_control, .data$n_treated))
               ) %>%
       drop_na(p_value) %>%
       mutate(diff = .data$mean_treated-.data$mean_control)

    return(t_test_mean_sd_result)
  }
  
  if(method == "moderated_t-test"){
  conditions_no_ref <- unique(pull(data, {{condition}}))[!unique(pull(data, {{condition}})) %in% ref_condition]
  
  message("[1/7] Creating moderated t-test input data ... ", appendLF = FALSE)
  moderated_t_test_input <- data %>%
    dplyr::distinct({{grouping}}, {{sample}}, {{intensity}}) %>%
    dplyr::arrange({{sample}}) %>% 
    tidyr::pivot_wider(names_from = {{sample}}, values_from = {{intensity}}) %>%
    tibble::column_to_rownames(var = rlang::as_name(rlang::enquo(grouping))) %>%
    as.matrix()   
    
  message("DONE", appendLF = TRUE)
  message("[2/7] Defining moderated t-test design ... ", appendLF = FALSE)
  
  moderated_t_test_map <- data %>%
    dplyr::distinct({{sample}}, {{condition}}) %>% 
    dplyr::mutate({{condition}} := paste0("x", {{condition}})) %>%
    dplyr::arrange({{sample}})
  
  moderated_t_test_design <- model.matrix(~ 0 + factor(dplyr::pull(moderated_t_test_map, {{condition}})))
  
  colnames(moderated_t_test_design) <- levels(factor(dplyr::pull(moderated_t_test_map, {{condition}})))
  
  message("DONE", appendLF = TRUE)
  message("[3/7] Fitting lmFit model ... ", appendLF = FALSE)
  
  moderated_t_test_fit <- suppressWarnings(lmFit(moderated_t_test_input, moderated_t_test_design))
  
  message("DONE", appendLF = TRUE)
  message("[4/7] Construc matrix of custom contrasts ... ", appendLF = FALSE)
  
  names <- map_chr(conditions_no_ref, ~paste0("x", .x, "_vs_x", ref_condition))
  comparisons <- map_chr(conditions_no_ref, ~paste0("x", .x, "-x", ref_condition))
  combinations <- map2(.x = names,
                      .y = comparisons,
                      .f = function(x, y){
                         deparse(expr(`=`(!!ensym(x), !!as_name(y))))
                       })
  flat_combinations <- purrr::reduce(combinations, ~ stringr::str_c(.x, .y, sep = ", ")) %>%
    paste0("makeContrasts(", ., ", levels = moderated_t_test_design)") %>%
    parse(text = .)
  
  contrast_matrix <- eval(flat_combinations)
  
  message("DONE", appendLF = TRUE)
  message("[5/7] Compute contrasts from linear model fit ... ", appendLF = FALSE)
  
  moderated_t_test_fit2 <- contrasts.fit(moderated_t_test_fit, contrast_matrix)
  
  message("DONE", appendLF = TRUE)
  message("[6/7] Compute empirical Bayes statistics ... ", appendLF = FALSE)
  
  moderated_t_test_fit3 <- eBayes(moderated_t_test_fit2)
  
  message("DONE", appendLF = TRUE)
  message("[7/7] Create result table ... ", appendLF = FALSE)
  
  moderated_t_test_missingness <- 
    assign_missingness(
      data,
      sample = {{sample}},
      condition = {{condition}},
      grouping = {{grouping}},
      intensity = {{intensity}},
      ref_condition = ref_condition
    ) %>% 
    drop_na(missingness, {{intensity}}) %>%
    group_by(comparison, {{grouping}}) %>%
    mutate(n_obs = n()) %>%
    ungroup() %>%
    distinct({{grouping}}, comparison, missingness, n_obs)
  
  moderated_t_test_result <- map_df(.x = names,
                .f = ~ topTable(moderated_t_test_fit3, coef = .x, number = Inf, confint = TRUE, sort.by = "p", adjust.method = p_adj_method) %>%
                  tibble::rownames_to_column(rlang::as_name(rlang::enquo(grouping))) %>%
                  mutate(comparison = .x)
  ) %>%
    mutate(comparison = stringr::str_replace_all(comparison, pattern = "^x|(?<=_vs_)x", replacement = "")) %>%
    rename(diff = logFC, CI_2.5 = CI.L, CI_97.5 = CI.R, t_statistic = t, avg_abundance = AveExpr, pval = P.Value, adj_pval = adj.P.Val) %>%
    left_join(moderated_t_test_missingness, by = c(rlang::as_name(rlang::enquo(grouping)), "comparison"))
  
  message("DONE", appendLF = TRUE)
  
  if(filter_NA_missingness == TRUE){
    moderated_t_test_result <- moderated_t_test_result %>%
      tidyr::drop_na(missingness)
    return(moderated_t_test_result)
  }
  if(filter_NA_missingness == FALSE) {
    return(moderated_t_test_result)
  }
  }
  
  if(method == "proDA"){
  message("[1/5] Creating proDA input data ... ", appendLF = FALSE)
    
  proDA_input <- data %>%
    dplyr::distinct({{grouping}}, {{sample}}, {{intensity}}) %>%
    dplyr::arrange({{sample}}) %>% 
    tidyr::pivot_wider(names_from = {{sample}}, values_from = {{intensity}}) %>%
    tibble::column_to_rownames(var = rlang::as_name(rlang::enquo(grouping))) %>%
    as.matrix()
  
  message("DONE", appendLF = TRUE)
  message("[2/5] Defining proDA design ... ", appendLF = FALSE)
  
  proDA_map <- data %>%
    dplyr::distinct({{sample}}, {{condition}}) %>% 
    dplyr::arrange({{sample}})
  
  proDA_design <- dplyr::pull(proDA_map, {{condition}})
  
  message("DONE", appendLF = TRUE)
  message("[3/5] Fitting proDA model (can take a few minutes) ... ", appendLF = FALSE)
  
  proDA_fit <- 
    proDA::proDA(proDA_input,
          design = proDA_design,
          reference_level = ref_condition)
  
  proDA_result_names <-
    proDA::result_names(proDA_fit)[!proDA::result_names(proDA_fit) %in% c("Intercept")]
  
  message("DONE", appendLF = TRUE)
  message("[4/5] Define missingness levels for filtering ... ", appendLF = FALSE)
  
  proDA_missingness <- 
    assign_missingness(
      data,
      sample = {{sample}},
      condition = {{condition}},
      grouping = {{grouping}},
      intensity = {{intensity}},
      ref_condition = ref_condition
    ) %>% 
    drop_na(missingness, {{intensity}}) %>%
    group_by(comparison, {{grouping}}) %>%
    mutate(n_obs = n()) %>%
    ungroup() %>%
    distinct({{grouping}}, comparison, missingness, n_obs)
  
  proDA_filter <- proDA_missingness %>%
    distinct({{grouping}}, .data$comparison) %>% 
    split(.$comparison) %>%
    map(select, -.data$comparison)
  
  message("DONE", appendLF = TRUE)
  message("[5/5] Extracting differential abundance from model and apply filter ... ", appendLF = FALSE)
  
  proDA_result <- proDA_result_names %>%
    map(~ proDA_fit) %>%
    set_names(nm = proDA_result_names) %>%
    map2(.y = names(.),
         .f = ~ proDA::test_diff(.x, contrast = .y, sort_by = "adj_pval")) 
  
  if(filter_NA_missingness == TRUE) {
  proDA_result <- proDA_result %>% 
    map2(.y = proDA_filter, 
         .f = ~ inner_join(.x, .y, by = c("name" = as_label(enquo(grouping))))) %>%
    map2(.y = names(.), 
         .f = ~ mutate(.x, comparison = str_replace_all(.y, pattern = "`", replacement = ""))) %>% 
    map_dfr(~ mutate(.x, adj_pval = p.adjust(pval, method = p_adj_method))) %>%
    select(-n_obs) %>%
    left_join(proDA_missingness, by = c("name" = rlang::as_name(rlang::enquo(grouping)), "comparison"))
  
  message("DONE", appendLF = TRUE)
  
  return(proDA_result)
  }
  if(filter_NA_missingness == FALSE) {
    proDA_result <- proDA_result %>% 
      map2(.y = names(.), 
           .f = ~ mutate(.x, comparison = str_replace_all(.y, pattern = "`", replacement = ""))) %>% 
      map_dfr(~ mutate(.x, adj_pval = p.adjust(pval, method = p_adj_method))) %>%
      select(-n_obs) %>%
      left_join(proDA_missingness, by = c("name" = rlang::as_name(rlang::enquo(grouping)), "comparison"))
    
    message("DONE", appendLF = TRUE)
    
    return(proDA_result)
  }
  }
}