fit_drc_4p <- function(data, sample, grouping, response, dose){
  input <- data %>%
    dplyr::distinct({{sample}}, {{grouping}}, {{response}}, {{dose}}) %>%
    tidyr::drop_na({{response}}) %>%
    dplyr::group_by({{grouping}}) %>%
    dplyr::mutate(n_concentrations = n_distinct(!!ensym(dose))) %>%
    dplyr::ungroup() %>%
    split(., pull(., !!ensym(grouping))) %>%
    purrr::keep( ~ unique(.x$n_concentrations) > 4) 
  
  pb <- progress::progress_bar$new(total = length(input), format = "  Model fitting [:bar] :percent :eta")
  fit_objects <- input %>%
    purrr::map(~ drc_4p(data = ., {{response}}, {{dose}}, pb))
  
  pb <- progress::progress_bar$new(total = length(input), format = "  1/4 Calculating correlations [:bar] :percent :eta")
  correlation_output <- input%>%
    map2(fit_objects, function(x, y) {
      pb$tick()
      tryCatch({
        suppressWarnings(cor(pull(x, !!ensym(response)), stats::fitted(y), method = "pearson"))
      }, error = function(error){
        c("no_correlation")
      })
    })%>%
    keep(~ . != "no_correlation") %>%
    map( ~ tibble(correlation = .)) %>%
    map2_dfr(names(.), function(x, y)
      mutate(x, sequence = y))
  
  pb <- progress::progress_bar$new(total = length(input), format = "  2/4 Calculating p-values [:bar] :percent :eta")
  p_value_correlation <- input%>%
    map2(fit_objects, function(x, y) {
      pb$tick()
      tryCatch({
        suppressWarnings(cor.test(pull(x, !!ensym(response)), fitted(y), method = "pearson")$p.value)
      }, error = function(error) {
        c("no_pvalue")
      })
    })%>%
    keep(~ . != "no_pvalue") %>%
    map( ~ tibble(p_value = .)) %>%
    map2_dfr(names(.), function(x, y)
      mutate(x, sequence = y))
  
  predictions_range <- input %>%
    purrr::map(~ filter(., {{dose}} != 0)) %>%
    purrr::map( ~ tibble(min = min(pull(.,!!ensym(dose))), max = max(pull(.,!!ensym(dose))))) %>%
    purrr::map( ~ expand.grid(dose = exp(seq(log(.$max), log(.$min), length = 100))))
  
  pb <- progress::progress_bar$new(total = length(fit_objects), format = "  3/4 Predicting curves for plots [:bar] :percent :eta")
  predictions <- fit_objects %>%
    purrr::map2(predictions_range, function(x, y){
      pb$tick()
      tryCatch({suppressWarnings(predict(x, newdata = y, interval = "confidence"))},
               error = function(error){
                 tibble(Prediction = "no_fit")
               })
    }) %>%
    purrr::map(~ as_tibble(.x))
  
  line_fit <- predictions_range %>%
    purrr::map(~ as_tibble(.x)) %>%
    purrr::map2(predictions, function(x, y){
      bind_cols(x, y)
    }) %>%
    keep(~ !("no_fit" %in% unique(.x$Prediction)) )
  
  names_to_keep <- names(line_fit)
  
  pb <- progress::progress_bar$new(total = length(line_fit), format = "  4/4 Plotting curves [:bar] :percent :eta")
  plots <- input %>%
    .[names_to_keep] %>%
    purrr::map2(line_fit, function(x, y){ 
      pb$tick()
      ggplot(data = x, aes(x = {{dose}}, y = {{response}})) +
        geom_point() +
        suppressWarnings(geom_ribbon(data = y, aes(x = dose, y = Prediction, ymin = Lower, ymax = Upper), alpha = 0.2)) +
        geom_line(data = y, aes(x=dose, y = Prediction)) +
        labs(title = names(.)) +
        theme_bw()+
        theme(axis.text.x = element_text(angle = 45, hjust =1)) +
        scale_x_log10() 
      }
    )
  
  correlation_output <- correlation_output %>%
    left_join(p_value_correlation, by = "sequence")
  
  correlation_output <- fit_objects %>%
    map(function(x) {
      bind_rows(pluck(x, "coefficients"))}) %>%
    map2_dfr(names(.), function(x, y)
      mutate(x, sequence = y)) %>%
    left_join(correlation_output, by = "sequence")
  
  result <- list(fit_objects = fit_objects, correlations = correlation_output, plots = plots)
  result
}