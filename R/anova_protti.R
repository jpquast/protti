anova_protti <- function(data, grouping, condition, mean_ratio, sd, n){
  
  pb <- progress::progress_bar$new(total = length(unique(dplyr::pull(data, {{grouping}}))), format = "  Calculating ANOVA [:bar] :current/:total (:percent) :eta")  
  
  result <- data %>% 
    dplyr::distinct({{grouping}}, {{condition}}, {{mean_ratio}}, {{sd}}, {{n}}) %>% 
    dplyr::group_by({{grouping}}) %>% 
    dplyr::filter({{n}} != 0) %>% 
    dplyr::mutate(n_groups = dplyr::n_distinct(!!ensym(condition))) %>% 
    dplyr::mutate(grand_mean = mean({{mean_ratio}})) %>% 
    dplyr::mutate(total_n = sum({{n}})) %>% 
    dplyr::mutate(ms_group = sum(({{mean_ratio}} - .data$grand_mean)^2 * {{n}}) / (.data$n_groups - 1)) %>%
    dplyr::mutate(ms_error = sum({{sd}}^2 * ({{n}} - 1))/(.data$total_n - .data$n_groups)) %>%
    dplyr::mutate(f = ms_group/ms_error) %>%
    dplyr::mutate(pval = stats::pf(.data$f, .data$n_groups - 1, .data$total_n - .data$n_groups, lower.tail = FALSE))%>%
    dplyr::distinct({{grouping}}, .data$ms_group, .data$ms_error, .data$f, .data$pval) %>%
    dplyr::mutate(adj_pval = stats::p.adjust(pval, method = "BH"))
  
  result
}
