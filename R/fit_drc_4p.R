#' Fitting four-parameter dose response curves
#'
#' Function for fitting four-parameter dose response curves for each group (precursor, peptide or protein). In addition it can
#' filter data based on completeness, the completeness distribution and statistical testing using ANOVA. 
#'
#' @param data A data frame containing at least the input variables.
#' @param sample The name of the column containing the sample names.
#' @param grouping The name of the column containing precursor, peptide or protein identifiers.
#' @param response The name of the column containing response values, eg. log2 transformed intensities.
#' @param dose The name of the column containing dose values, eg. the treatment concentrations.
#' @param filter A character vector indicating if models should be filtered and if they should be filtered before or after the curve
#' fits. Filtering of models can be skipped with \code{filter = "none"}. The only criteria always used is that at least 5 data 
#' points are present for any peptide. Otherwise no proper curves can be fit. Data can be filtered prior to model fitting with 
#' \code{filter = "pre"}. In that case models will only be fitted for data that passed the filtering step. This will allow for faster
#' model fitting since only fewer models will be fit. If you plan on performing an enrichment analysis you have to choose 
#' \code{filter = "post"}. All models will be fit (even the ones that do not pass the filtering criteria) the subset of models that do 
#' not fulfill the filtering criteria. The models that are no hits should be the ones used in enrichment analysis in combination with 
#' the models that are hits. Therefore for post-filtering the full list is returned and it will only contain annotations that 
#' indicate if the filtering was passed or not. Default is "post". For ANOVA an adjusted p-value of 0.05 is used as a cutoff.
#' @param replicate_completeness Similar to \code{completenss_MAR} of the \code{assign_missingness} function this argument sets a 
#' threshold for the completeness of data. In contrast to \code{assign_missingness} it only determines the completeness for one 
#' condition and not the comparison of two conditions. The threshold is used to calculate a minimal degree of data completeness. 
#' The value provided to this argument has to be between 0 and 1, default is 0.7. It is multiplied with the number of replicates 
#' and then adjusted downward. The resulting number is the minimal number of observations that a condition needs to have to be considered
#' "complete enough" for the \code{condition_completeness} argument.
#' @param condition_completeness This argument determines how many conditions need to at least fulfill the "complete enough" criteria 
#' set with \code{replicate_completeness}. The value provided to this argument has to be between 0 and 1, default is 0.5. It is 
#' multiplied with the number of conditions and then adjusted downward. The resulting number is the minimal number of conditions that
#' need to fulfill the \code{replicate_completeness} argument for a peptide to pass the filtering.
#' @param log_logarithmic logical indicating if a logarithmic or log-logarithmic model is fitted. 
#' If response values form a symmetric curve for non-log transformed dose values, a logarithmic model instead
#' of a log-logarithmic model should be used. Usually biological dose response data has a log-logarithmic distribution, which is the 
#' reason this is the default. Log-logarithmic models are symmetric if dose values are log transformed. 
#' @param include_models A logical indicating if model fit objects should be exported. These are usually very large and not necessary for further analysis.
#' @param retain_columns A vector indicating if certain columns should be retained from the input data frame. Default is not retaining 
#' additional columns \code{retain_columns = NULL}. Specific columns can be retained by providing their names (not in quotations marks, 
#' just like other column names, but in a vector).
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
#' @importFrom stats p.adjust
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
fit_drc_4p <- function(data, sample, grouping, response, dose, filter = "post", replicate_completeness = 0.7, condition_completeness = 0.5, log_logarithmic = TRUE, include_models = FALSE, retain_columns = NULL){
  if (!requireNamespace("drc", quietly = TRUE)) {
    stop("Package \"drc\" is needed for this function to work. Please install it.", call. = FALSE)
  }
  # to prevent no visible binding for global variable '.' note.
  . = NULL
  
  # preprocessing of data 
  data_prep <- data %>% 
    dplyr::ungroup() %>% 
    dplyr::distinct({{sample}}, {{grouping}}, {{response}}, {{dose}}) %>%
    tidyr::complete(nesting(!!ensym(sample), !!ensym(dose)), !!ensym(grouping)) 
  
  if(filter != "none"){
    n_conditions = length(unique(dplyr::pull(data_prep, {{ dose }})))
    n_replicates = length(unique(dplyr::pull(data_prep, {{ sample }}))) / n_conditions
    n_replicates_completeness <- floor(replicate_completeness * n_replicates)
    n_conditions_completeness <- floor(condition_completeness * n_conditions)
    
    # perform anova on groups
    anova <- data_prep %>% 
      dplyr::group_by({{ grouping }}, {{ dose }}) %>% 
      dplyr::mutate(
        n = n_replicates,
        mean_ratio = mean({{ response }}, na.rm = TRUE),
        sd = sd({{ response }}, na.rm = TRUE)) %>% 
      dplyr::distinct({{ grouping }}, {{ dose }}, .data$mean_ratio, .data$sd, .data$n) %>%
      tidyr::drop_na(.data$mean_ratio, .data$sd) %>% 
      anova_protti({{ grouping }}, {{ dose }}, .data$mean_ratio, .data$sd, .data$n) %>% 
      dplyr::distinct({{ grouping }}, .data$pval) %>% 
      dplyr::mutate(anova_adj_pval = stats::p.adjust(.data$pval, method = "BH")) %>% 
      dplyr::rename(anova_pval = .data$pval)
    
    # extract elements that pass anova significant threshold
    anova_filtered <- anova %>% 
      dplyr::filter(.data$anova_adj_pval <= 0.05) %>% 
      dplyr::pull({{ grouping }})
    
    # filter for data completeness based on replicates and conditions and based on the distribution of condition missingness
    up <- ceiling(n_conditions * 0.5)
    low <- floor(n_conditions * 0.5)
    
    # create vectors that define which distribution of data missingness is accepted and should be retained
    vector <- c(up, low)
    vector_rev <- rev(vector)
    vector_add <- c(vector[1] + 1, vector[2] - 1)
    vector_add_rev <- rev(vector_add)
    
    concentrations <- sort(unique(dplyr::pull(data_prep, {{ dose }})))
    
    filter_completeness <- data_prep %>% 
      dplyr::group_by({{ grouping }}, {{ dose }}) %>% 
      dplyr::mutate(enough_replicates = sum(!is.na({{ response }})) >= n_replicates_completeness) %>% 
      dplyr::group_by({{ grouping }}) %>% 
      dplyr::mutate(enough_conditions = sum(.data$enough_replicates) /n_replicates >= n_conditions_completeness) %>% 
      dplyr::mutate(anova_significant = {{ grouping }} %in% anova_filtered) %>% 
      dplyr::mutate(lower_vector = {{ dose }} %in% concentrations[1:vector[1]],
                    lower_vector_rev = {{ dose }} %in% concentrations[1:vector_rev[1]],
                    lower_vector_add = {{ dose }} %in% concentrations[1:vector_add[1]],
                    lower_vector_add_rev = {{ dose }} %in% concentrations[1:vector_add_rev[1]]) %>% 
      dplyr::group_by({{ grouping }}, .data$lower_vector) %>% 
      dplyr::mutate(mean_vector = mean({{ response }}, na.rm = TRUE)) %>% 
      dplyr::group_by({{ grouping }}, .data$lower_vector_rev) %>% 
      dplyr::mutate(mean_vector_rev = mean({{ response }}, na.rm = TRUE)) %>% 
      dplyr::group_by({{ grouping }}, .data$lower_vector_add) %>% 
      dplyr::mutate(mean_vector_add = mean({{ response }}, na.rm = TRUE)) %>% 
      dplyr::group_by({{ grouping }}, .data$lower_vector_add_rev) %>% 
      dplyr::mutate(mean_vector_add_rev = mean({{ response }}, na.rm = TRUE)) %>% 
      dplyr::group_by({{ grouping }}) %>% 
      dplyr::mutate(mean_vector = min(.data$mean_vector) == .data$mean_vector) %>% 
      dplyr::mutate(mean_vector_rev = min(.data$mean_vector_rev) == .data$mean_vector_rev) %>% 
      dplyr::mutate(mean_vector_add = min(.data$mean_vector_add) == .data$mean_vector_add) %>% 
      dplyr::mutate(mean_vector_add_rev = min(.data$mean_vector_add_rev) == .data$mean_vector_add_rev) %>% 
      dplyr::mutate(dose_MNAR = ifelse(all((.data$lower_vector & .data$enough_replicates == FALSE & .data$mean_vector)| (!.data$lower_vector & .data$enough_replicates == TRUE & !.data$mean_vector)) |
                                      all((.data$lower_vector & .data$enough_replicates == TRUE & !.data$mean_vector)| (!.data$lower_vector & .data$enough_replicates == FALSE & .data$mean_vector)) |
                                      all((.data$lower_vector_rev & .data$enough_replicates == FALSE & .data$mean_vector_rev)| (!.data$lower_vector_rev & .data$enough_replicates == TRUE & !.data$mean_vector_rev)) |
                                      all((.data$lower_vector_rev & .data$enough_replicates == TRUE & !.data$mean_vector_rev)| (!.data$lower_vector_rev & .data$enough_replicates == FALSE & .data$mean_vector_rev)) |
                                      all((.data$lower_vector_add & .data$enough_replicates == FALSE & .data$mean_vector_add)| (!.data$lower_vector_add & .data$enough_replicates == TRUE & !.data$mean_vector_add)) |
                                      all((.data$lower_vector_add & .data$enough_replicates == TRUE & !.data$mean_vector_add)| (!.data$lower_vector_add & .data$enough_replicates == FALSE & .data$mean_vector_add)) |
                                      all((.data$lower_vector_add_rev & .data$enough_replicates == FALSE & .data$mean_vector_add_rev)| (!.data$lower_vector_add_rev & .data$enough_replicates == TRUE & !.data$mean_vector_add_rev)) |
                                      all((.data$lower_vector_add_rev & .data$enough_replicates == TRUE & !.data$mean_vector_add_rev)| (!.data$lower_vector_add_rev & .data$enough_replicates == FALSE & .data$mean_vector_add_rev)), TRUE, FALSE)) %>% 
      dplyr::distinct({{ grouping }}, .data$enough_conditions, .data$anova_significant, .data$dose_MNAR) %>% 
      dplyr::mutate(passed_filter = (.data$enough_conditions == TRUE & .data$anova_significant == TRUE) | .data$dose_MNAR == TRUE)

    if(filter == "pre"){
      filtered_vector <- filter_completeness %>% 
        dplyr::filter(.data$passed_filter == TRUE) %>% 
        dplyr::pull({{ grouping }})
      
      data_prep <- data_prep %>% 
        filter({{ grouping }} %in% filtered_vector)
    }
  }
  # prepare data
  
  input <- data_prep %>%
    tidyr::drop_na({{response}}) %>%
    dplyr::group_by({{grouping}}) %>%
    dplyr::mutate(n_concentrations = dplyr::n_distinct(!!ensym(dose))) %>%
    dplyr::ungroup() %>%
    split(dplyr::pull(., !!ensym(grouping))) %>%
    purrr::keep(.p = ~ unique(.x$n_concentrations) > 4) # make sure that there are enough data points to even fit a curve. This is not really filtering.
  
  # fit models
  
  pb <- progress::progress_bar$new(total = length(input), format = " 1/4 Model fitting [:bar] :current/:total (:percent) :eta")
  fit_objects <- purrr::map(.x = input,
                            .f = ~ drc_4p(data = .x, response = {{response}}, dose = {{dose}}, log_logarithmic = log_logarithmic, pb = pb))
  
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
    purrr::keep(.p = ~ .x != "no_correlation" & !is.na(.x)) %>%
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
    purrr::keep(.p = ~ .x != "no_pvalue" & !is.na(.x)) %>%
    purrr::map(.f = ~ tibble::tibble(pval = .x)) %>%
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
  
  if(log_logarithmic == TRUE){
  predictions_range <- input %>%
    purrr::map(.f = ~ dplyr::filter(.x, {{dose}} != 0)) %>%
    purrr::map(.f = ~ tibble::tibble(min = min(pull(.x,!!ensym(dose))), max = max(pull(.x,!!ensym(dose))))) %>%
    purrr::map(.f = ~ expand.grid(dose = exp(seq(log(.x$max), log(.x$min), length = 100))))
  }
  if(log_logarithmic == FALSE){ # dose values are ideally not log transformed for logarithmic models
  predictions_range <- input %>%
    purrr::map(.f = ~ dplyr::filter(.x, {{dose}} != 0)) %>%
    purrr::map(.f = ~ tibble::tibble(min = min(pull(.x,!!ensym(dose))), max = max(pull(.x,!!ensym(dose))))) %>%
    purrr::map(.f = ~ expand.grid(dose = seq(.x$max, .x$min, length = 100)))
  }
  
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
    tidyr::nest(plot_points = c({{response}}, {{dose}})) %>% 
    dplyr::rename(hill_coefficient = .data$`hill:(Intercept)`,
                  min_model = .data$`min_value:(Intercept)`,
                  max_model = .data$`max_value:(Intercept)`,
                  ec_50 = .data$`ec_50:(Intercept)`) %>% 
    dplyr::ungroup()

  # post filter and addition of columns if prefilter was performed
  if(filter != "none"){
    output <- output %>% 
      dplyr::left_join(filter_completeness, by = rlang::as_name(rlang::enquo(grouping))) %>% 
      dplyr::left_join(anova, by = rlang::as_name(rlang::enquo(grouping))) %>% 
      dplyr::mutate(passed_filter = .data$passed_filter & .data$correlation >= 0.7 & .data$min_model < .data$max_model)
  }
  
  if(filter == "pre"){
    output <- output %>% 
      filter(passed_filter)
  }
  # return result

  if(!missing(retain_columns)){
  output <- data %>% 
    dplyr::select(!!enquo(retain_columns), colnames(output)[!colnames(output) %in% c("pval", "hill_coefficient", "min_model", "max_model", "ec_50", "correlation", "plot_curve", "plot_points", "enough_conditions", "anova_significant", "dose_MNAR", "anova_pval", "anova_adj_pval", "passed_filter")]) %>% 
    dplyr::distinct() %>% 
    dplyr::right_join(output, by = colnames(output)[!colnames(output) %in% c("pval", "hill_coefficient", "min_model", "max_model", "ec_50", "correlation", "plot_curve", "plot_points", "enough_conditions", "anova_significant", "dose_MNAR", "anova_pval", "anova_adj_pval", "passed_filter")]) %>% 
    dplyr::arrange(dplyr::desc(.data$correlation))
  }
  
  if(include_models == TRUE){ 
    combined_output <- list(fit_objects = fit_objects, correlations = output)
    return(combined_output)
  }
  if(include_models == FALSE){ 
    rm(fit_objects)
    return(output)
  }
}