#' Fitting four-parameter dose response curves
#'
#' Function for fitting four-parameter dose response curves for each group (precursor, peptide or protein). 
#'
#' @param data A data frame containing at least the input variables.
#' @param sample The name of the column containing the sample names.
#' @param grouping The name of the column containing precursor, peptide or protein identifiers.
#' @param response The name of the column containing response values, eg. log2 transformed intensities.
#' @param dose The name of the column containing dose values, eg. the treatment concentrations.
#' @param include_models A logical indicating if model fit objects should be exported. These are usually very large and not necessary for further analysis.
#' 
#' @return If \code{include_models = FALSE} a data frame is returned that contains correlations of predicted to measured values as a measure of the goodness of the curve fit,
#' an associated p-value and the four parameters of the model for each group. Furthermore, input data for plots is returned in the columns \code{plot_curve} (curve and confidence interval) and \code{plot_points} (measured points).
#' If \ code{include_models = TURE}, a list is returned that contains: 
#' \itemize{
#' \item{\code{fit_objects}: }{The fit objects of type \code{drc} for each group.}
#' \item{\code{correlations}: }{The correlation data frame described above}
#' } 
#' @import dplyr
#' @import tidyr
#' @import progress
#' @importFrom purrr keep map map2 map2_df pluck
#' @importFrom tibble tibble as_tibble
#' @importFrom rlang .data ensym as_name enquo :=
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
fit_drc_4p <- function(data, sample, grouping, response, dose, include_models = FALSE){
  # to prevent no visible binding for global variable '.' note.
  . = NULL
  
  # prepare data
  
  input <- data %>%
    dplyr::distinct({{sample}}, {{grouping}}, {{response}}, {{dose}}) %>%
    tidyr::drop_na({{response}}) %>%
    dplyr::group_by({{grouping}}) %>%
    dplyr::mutate(n_concentrations = dplyr::n_distinct(!!ensym(dose))) %>%
    dplyr::ungroup() %>%
    split(dplyr::pull(., !!ensym(grouping))) %>%
    purrr::keep(.p = ~ unique(.x$n_concentrations) > 4) 
  
  # fit models
  
  pb <- progress::progress_bar$new(total = length(input), format = " 1/4 Model fitting [:bar] :current/:total (:percent) :eta")
  fit_objects <- purrr::map(.x = input,
                            .f = ~ drc_4p(data = .x, {{response}}, {{dose}}, pb))
  
  # extract information from fit objects and calculate correlation
  
  pb <- progress::progress_bar$new(total = length(input), format = "  2/4 Calculating correlations [:bar] :current/:total (:percent) :eta")
  correlation_output <- input %>%
    purrr::map2(.y = fit_objects, 
                .f = ~ {pb$tick(); 
                  tryCatch({
                    suppressWarnings(stats::cor(pull(.x, !!ensym(response)), stats::fitted(.y), method = "pearson"))
                  }, error = function(error){
                    c("no_correlation")
                  })
                }) %>%
    purrr::keep(.p = ~ .x != "no_correlation") %>%
    purrr::map(.f = ~ tibble::tibble(correlation = .x)) %>%
    purrr::map2_df(.y = names(.),
                   .f = ~ dplyr::mutate(.x, {{grouping}} := .y)
    )
  
  # calculate p-values
  
  pb <- progress::progress_bar$new(total = length(input), format = "  3/4 Calculating p-values [:bar] :current/:total (:percent) :eta")
  p_value_correlation <- input %>%
    purrr::map2(.y = fit_objects, 
                .f = ~ {pb$tick();
                  tryCatch({
                    suppressWarnings(stats::cor.test(dplyr::pull(.x, !!ensym(response)), stats::fitted(.y), method = "pearson")$p.value)
                  }, error = function(error) {
                    c("no_pvalue")
                  })
                }) %>%
    purrr::keep(.p = ~ .x != "no_pvalue") %>%
    purrr::map(.f = ~ tibble::tibble(p_value = .x)) %>%
    purrr::map2_df(.y = names(.),
                   .f = ~ dplyr::mutate(.x, {{grouping}} := .y)
    )
  
  # creating correlation output data frame
  
  correlation_output <- correlation_output %>%
    dplyr::left_join(p_value_correlation, by = rlang::as_name(rlang::enquo(grouping)))
  
  correlation_output <- fit_objects %>%
    purrr::map(.f = ~ dplyr::bind_rows(purrr::pluck(.x, "coefficients"))) %>%
    purrr::map2_df(.y = names(.), 
                   .f = ~ dplyr::mutate(.x, {{grouping}} := .y)) %>%
    dplyr::left_join(correlation_output, by = rlang::as_name(rlang::enquo(grouping)))
  
  # starting calculations for plot input data frames
  
  predictions_range <- input %>%
    purrr::map(.f = ~ dplyr::filter(.x, {{dose}} != 0)) %>%
    purrr::map(.f = ~ tibble::tibble(min = min(pull(.x,!!ensym(dose))), max = max(pull(.x,!!ensym(dose))))) %>%
    purrr::map(.f = ~ expand.grid(dose = exp(seq(log(.x$max), log(.x$min), length = 100))))
  
  pb <- progress::progress_bar$new(total = length(fit_objects), format = "  4/4 Predicting curves for plots [:bar] :current/:total (:percent) :eta")
  predictions <- fit_objects %>%
    purrr::map2(.y = predictions_range, 
                .f = ~{pb$tick();
                  tryCatch({suppressWarnings(stats::predict(.x, newdata = .y, interval = "confidence"))},
                           error = function(error){
                             tibble::tibble(Prediction = "no_fit")
                           })
                }) %>%
    purrr::map(.f = ~ tibble::as_tibble(.x))
  
  line_fit <- predictions_range %>%
    purrr::map(.f = ~ tibble::as_tibble(.x)) %>%
    purrr::map2(.y = predictions, 
                .f = function(x, y){bind_cols(x, y)}
    ) %>%
    purrr::keep(.p = ~ !("no_fit" %in% unique(.x$Prediction))) %>% 
    purrr::map2_df(.y = names(.),
                   .f = ~ dplyr::mutate(.x, {{grouping}} := .y))
  
  plot_points <- input[unique(dplyr::pull(line_fit, {{grouping}}))] %>% 
    purrr::map2_df(.y = names(.),
                   .f = ~ dplyr::mutate(.x, {{grouping}} := .y) %>% 
                     dplyr::select({{grouping}}, {{response}}, {{dose}}))
  
  # combining correlations with information for plot
  
  output <- correlation_output %>%
    dplyr::left_join(line_fit, by = rlang::as_name(rlang::enquo(grouping))) %>% 
    dplyr::group_by({{grouping}}) %>% 
    tidyr::nest(plot_curve = c(.data$dose, .data$Prediction, .data$Lower, .data$Upper)) %>% 
    dplyr::left_join(plot_points, by = rlang::as_name(rlang::enquo(grouping))) %>% 
    tidyr::nest(plot_points = c({{response}}, {{dose}})) 
  
  # return result
  
  if(include_models == TRUE){ 
    combined_output <- list(fit_objects = fit_objects, correlations = output)
    return(combined_output)
  }
  if(include_models == FALSE){ 
    rm(fit_objects)
    return(output)
  }
}