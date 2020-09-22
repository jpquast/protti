#' Fitting four-parameter dose response curves
#'
#' Function for fitting four-parameter dose response curves for each group (precursor, peptide or protein). 
#'
#' @param data a data frame containing at least the input variables.
#' @param sample the name of the column containing the sample names.
#' @param grouping the name of the column containing precursor, peptide or protein identifiers.
#' @param response the name of the column containing response values, eg. log2 transformed intensities.
#' @param dose the name of the column containing dose values, eg. the treatment concentrations.
#' 
#' @return A list that contains: 
#' \itemize{
#' \item{\code{fit_objects}: }{The fit objects of type \code{drc} for each group. }
#' \item{\code{correlations}: }{A data frame that contains correlation of predicted to measured values as a measure of the goodness of the curve fit, an associated p-value and the four parameters of the model for each group.}
#' \item{\code{plots}: }{The ggplot objects containing the curve plot. Printing these will print the plot.}
#' } 
#' @import dplyr
#' @import tidyr
#' @import purrr
#' @import progress
#' @import progressor
#' @import ggplot2
#' @import furrr
#' @importFrom tibble tibble
#' @importFrom tibble as_tibble
#' @importFrom rlang .data
#' @importFrom magrittr %>%
#' @export
#'
#' @examples
#' \dontrun{
#' fit_drc_4p(
#' data,
#' sample = r_file_name,
#' grouping = eg_precursor_id,
#' response = intensity,
#' dose = concentration
#' )
#' }
fit_drc_4p <- function(data, sample, grouping, response, dose, parallel = FALSE){
  # to prevent no visible binding for global variable '.' note.
  . = NULL
  # normalised_intensity_log2 <- NULL
  # concentration <- NULL
  
  input <- data %>%
    dplyr::distinct({{sample}}, {{grouping}}, {{response}}, {{dose}}) %>%
    tidyr::drop_na({{response}}) %>%
    dplyr::group_by({{grouping}}) %>%
    dplyr::mutate(n_concentrations = dplyr::n_distinct(!!ensym(dose))) %>%
    dplyr::ungroup() %>%
    split(., dplyr::pull(., !!ensym(grouping))) %>%
    purrr::keep( ~ unique(.x$n_concentrations) > 4) 
  
  if(parallel == FALSE) { 
    pb <- progress::progress_bar$new(total = length(input), format = "  Model fitting [:bar] :percent :eta")
  fit_objects <- input %>%
    purrr::map(.f = ~ drc_4p(data = ., {{response}}, {{dose}}, pb)
               )
  } else {

  p <- progressor(steps = length(input))
  fit_objects <- furrr::future_map(.x = input,
                                   .f = ~{p(); drc_4p(data = .x, {{response}}, {{dose}})}, 
                                   .options = future_options(globals = FALSE)
      )
  }
  
  # if(parallel == FALSE){
  # pb <- progress::progress_bar$new(total = length(input), format = "  1/4 Calculating correlations [:bar] :percent :eta")
  # correlation_output <- input%>%
  #   purrr::map2(fit_objects, function(x, y) {
  #     pb$tick()
  #     tryCatch({
  #       suppressWarnings(stats::cor(pull(x, !!ensym(response)), stats::fitted(y), method = "pearson"))
  #     }, error = function(error){
  #       c("no_correlation")
  #     })
  #   })%>%
  #   purrr::keep(~ . != "no_correlation") %>%
  #   purrr::map( ~ tibble::tibble(correlation = .)) %>%
  #   purrr::map2_dfr(names(.), function(x, y)
  #     dplyr::mutate(x, sequence = y))
  # } else {
  #   correlation_output <- input%>%
  #     furrr::future_map2(.y = fit_objects, 
  #                        .f = function(x, y) {
  #       tryCatch({
  #         suppressWarnings(stats::cor(pull(x, !!ensym(response)), stats::fitted(y), method = "pearson"))
  #       }, error = function(error){
  #         c("no_correlation")
  #       })
  #     },
  #     .progress = TRUE) %>%
  #     purrr::keep(.f = ~ . != "no_correlation") %>%
  #     furrr::future_map(.f = ~ tibble::tibble(correlation = .)) %>%
  #     furrr::future_map2_dfr(.y = names(.), 
  #                            .f = ~ .x %>% 
  #                              dplyr::mutate(sequence = .y))
  # }
  # 
  # pb <- progress::progress_bar$new(total = length(input), format = "  2/4 Calculating p-values [:bar] :percent :eta")
  # p_value_correlation <- input%>%
  #   purrr::map2(fit_objects, function(x, y) {
  #     pb$tick()
  #     tryCatch({
  #       suppressWarnings(stats::cor.test(dplyr::pull(x, !!ensym(response)), stats::fitted(y), method = "pearson")$p.value)
  #     }, error = function(error) {
  #       c("no_pvalue")
  #     })
  #   })%>%
  #   purrr::keep(~ . != "no_pvalue") %>%
  #   purrr::map( ~ tibble::tibble(p_value = .)) %>%
  #   purrr::map2_dfr(names(.), function(x, y)
  #     dplyr::mutate(x, sequence = y))
  # 
  # predictions_range <- input %>%
  #   purrr::map(~ dplyr::filter(., {{dose}} != 0)) %>%
  #   purrr::map( ~ tibble::tibble(min = min(pull(.,!!ensym(dose))), max = max(pull(.,!!ensym(dose))))) %>%
  #   purrr::map( ~ expand.grid(dose = exp(seq(log(.$max), log(.$min), length = 100))))
  # 
  # pb <- progress::progress_bar$new(total = length(fit_objects), format = "  3/4 Predicting curves for plots [:bar] :percent :eta")
  # predictions <- fit_objects %>%
  #   purrr::map2(predictions_range, function(x, y){
  #     pb$tick()
  #     tryCatch({suppressWarnings(stats::predict(x, newdata = y, interval = "confidence"))},
  #              error = function(error){
  #                tibble::tibble(Prediction = "no_fit")
  #              })
  #   }) %>%
  #   purrr::map(~ tibble::as_tibble(.x))
  # 
  # line_fit <- predictions_range %>%
  #   purrr::map(~ tibble::as_tibble(.x)) %>%
  #   purrr::map2(predictions, function(x, y){
  #     bind_cols(x, y)
  #   }) %>%
  #   purrr::keep(~ !("no_fit" %in% unique(.x$Prediction)) )
  # 
  # names_to_keep <- names(line_fit)
  # 
  # pb <- progress::progress_bar$new(total = length(line_fit), format = "  4/4 Plotting curves [:bar] :percent :eta")
  # plots <- purrr::pmap(list(x = input[names_to_keep], y = line_fit, z = names(line_fit)), function(x, y, z){
  #   pb$tick()
  #   ggplot2::ggplot(data = x, ggplot2::aes(x = {{dose}}, y = {{response}})) +
  #     ggplot2::geom_point() +
  #     suppressWarnings(ggplot2::geom_ribbon(data = y, ggplot2::aes(x = dose, y = .data$Prediction, ymin = .data$Lower, ymax = .data$Upper), alpha = 0.2)) +
  #     ggplot2::geom_line(data = y, ggplot2::aes(x=dose, y = .data$Prediction)) +
  #     ggplot2::labs(title = z) +
  #     ggplot2::theme_bw()+
  #     ggplot2::theme(axis.text.x = element_text(angle = 45, hjust =1)) +
  #     ggplot2::scale_x_log10()
  # }
  # )
  # 
  # correlation_output <- correlation_output %>%
  #   dplyr::left_join(p_value_correlation, by = "sequence")
  # 
  # correlation_output <- fit_objects %>%
  #   purrr::map(function(x) {
  #     dplyr::bind_rows(purrr::pluck(x, "coefficients"))}) %>%
  #   purrr::map2_dfr(names(.), function(x, y)
  #     dplyr::mutate(x, sequence = y)) %>%
  #   dplyr::left_join(correlation_output, by = "sequence")
  # 
  # result <- list(fit_objects = fit_objects, correlations = correlation_output, plots = plots)
  # result
}