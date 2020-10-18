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
#' @param n_cores Optional, the number of cores used if workers are set up manually.
#' 
#' @return A data frame is returned that contains correlations of predicted to measured values as a measure of the goodness of the curve fit, 
#' an associated p-value and the four parameters of the model for each group. Furthermore, input data for plots is returned in the columns \code{plot_curve} 
#' (curve and confidence interval) and \code{plot_points} (measured points).
#' 
#' @import dplyr
#' @importFrom furrr future_map_dfr
#' @importFrom rlang .data as_name enquo
#' @importFrom magrittr %>%
#' @importFrom parallel detectCores
#' @importFrom future plan
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
parallel_fit_drc_4p <- function(data, sample, grouping, response, dose, n_cores = NULL){
  . = NULL
  multiprocess = NULL
  sequential = NULL
  terminate = FALSE
  
  if(missing(n_cores)){
    message("Setting up workers ... ", appendLF = FALSE)
    terminate = TRUE
    n_cores <- parallel::detectCores()
    n_cores <- ifelse(n_cores > 12, 12, n_cores)
    future::plan(multiprocess)
    message("DONE", appendLF = TRUE)
  }
  
  pieces <- rep(1:n_cores, round(length(unique(dplyr::pull(data, {{grouping}})))/n_cores) +1)[1:length(unique(dplyr::pull(data, {{grouping}})))]
  
  pieces_mapping <- tibble::tibble({{grouping}} := unique(dplyr::pull(data, {{grouping}})), piece = pieces)
  
  input <- data %>%
    dplyr::left_join(pieces_mapping, by = rlang::as_name(rlang::enquo(grouping))) %>%
    split(dplyr::pull(., .data$piece))
  
  message("Performing model fit (this may take a while) ... ", appendLF = FALSE)
  
  result <- furrr::future_map_dfr(.x = input,
                                  .f = ~ protti::fit_drc_4p(.x, {{sample}}, {{grouping}}, {{response}}, {{dose}}, include_models = FALSE),
                                  .options = furrr::future_options(globals = FALSE)
  )
  
  message("DONE", appendLF = TRUE)
  
  if(terminate == TRUE){
    future::plan(sequential)
  }
  
  result
}