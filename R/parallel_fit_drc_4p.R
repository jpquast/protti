#' Fitting four-parameter dose response curves (using parallel processing)
#'
#' This function is a wrapper around \code{fit_drc_4p} that allows the use of all system cores for model fitting. It should only be used on systems that have enough memory available. 
#' Workers can either be set up manually before running the function with \code{future::plan(multiprocess)} or automatically by the function (maximum number of workers is 12 in this case). If workers are set up manually the 
#' number of cores should be provided to \code{n_cores}. Worker can be terminated after completion with \code{future::plan(sequential)}. It is not possible to export the 
#' individual fit objects when using this function as compared to the non parallel function as they are too large for efficient export from the workers. 
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
#' indicate if the filtering was passed or not. Default is "post".
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
#' @param n_cores Optional, the number of cores used if workers are set up manually.
#' 
#' @return A data frame is returned that contains correlations of predicted to measured values as a measure of the goodness of the curve fit, 
#' an associated p-value and the four parameters of the model for each group. Furthermore, input data for plots is returned in the columns \code{plot_curve} 
#' (curve and confidence interval) and \code{plot_points} (measured points).
#' 
#' @import dplyr
#' @importFrom rlang .data as_name enquo
#' @importFrom magrittr %>%
#' @export
#'
#' @examples
#' \dontrun{
#' parallel_fit_drc_4p(
#' data,
#' sample = r_file_name,
#' grouping = eg_precursor_id,
#' response = intensity,
#' dose = concentration
#' )
#' }
parallel_fit_drc_4p <- function(data, sample, grouping, response, dose, filter = "post", replicate_completeness = 0.7, condition_completeness = 0.5, log_logarithmic = TRUE, retain_columns = NULL, n_cores = NULL){
  dependency_test <- c(furrr = !requireNamespace("furrr", quietly = TRUE), future = !requireNamespace("future", quietly = TRUE), parallel = !requireNamespace("parallel", quietly = TRUE))
  if (any(dependency_test)) {
    dependency_name <- names(dependency_test[dependency_test == TRUE])
    if(length(dependency_name) == 1){
      stop("Package \"", paste(dependency_name), "\" is needed for this function to work. Please install it.", call. = FALSE)
    } else{
      stop("Packages \"", paste(dependency_name, collapse = "\" and \""), "\" are needed for this function to work. Please install them.", call. = FALSE)
    }
  }
  . = NULL
  terminate = FALSE
  
  if(missing(n_cores)){
    message("Setting up workers ... ", appendLF = FALSE)
    terminate = TRUE
    n_cores <- parallel::detectCores()
    n_cores <- ifelse(n_cores > 12, 12, n_cores)
    future::plan(future::multiprocess, workers = n_cores)
    message("DONE", appendLF = TRUE)
  }
  
  pieces <- rep(1:n_cores, round(length(unique(dplyr::pull(data, {{grouping}})))/n_cores) +1)[1:length(unique(dplyr::pull(data, {{grouping}})))]
  
  pieces_mapping <- tibble::tibble({{grouping}} := unique(dplyr::pull(data, {{grouping}})), piece = pieces)
  
  input <- data %>%
    dplyr::left_join(pieces_mapping, by = rlang::as_name(rlang::enquo(grouping))) %>%
    split(dplyr::pull(., .data$piece))

  message("Performing model fit (this may take a while) ... ", appendLF = FALSE)
  
  result <- furrr::future_map_dfr(.x = input,
                                  .f = ~ protti::fit_drc_4p(.x, sample = {{sample}}, grouping = {{grouping}}, response = {{response}}, dose = {{dose}}, filter = filter, replicate_completeness = replicate_completeness, condition_completeness = condition_completeness, log_logarithmic = log_logarithmic, retain_columns = retain_columns, include_models = FALSE),
                                  .options = furrr::future_options(globals = FALSE)
  )
  
  message("DONE", appendLF = TRUE)
  
  if(terminate == TRUE){
    future::plan(future::sequential)
  }
  
  result %>% 
    dplyr::arrange(desc(.data$correlation))
}