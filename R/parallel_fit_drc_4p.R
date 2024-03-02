#' Fitting four-parameter dose response curves (using parallel processing)
#'
#' This function is a wrapper around \code{fit_drc_4p} that allows the use of all system cores for
#' model fitting. It should only be used on systems that have enough memory available. Workers can
#' either be set up manually before running the function with \code{future::plan(multisession)} or
#' automatically by the function (maximum number of workers is 12 in this case). If workers are set
#' up manually the number of cores should be provided to \code{n_cores}. Worker can be terminated
#' after completion with \code{future::plan(sequential)}. It is not possible to export the
#' individual fit objects when using this function as compared to the non parallel function as
#' they are too large for efficient export from the workers.
#'
#' @details If data filtering options are selected, data is annotated based on multiple criteria.
#' If `"post"` is selected the data is annotated based on completeness, the completeness distribution, the
#' adjusted ANOVA p-value cutoff and a correlation cutoff. Completeness of features is determined based on
#' the `n_replicate_completeness` and `n_condition_completeness` arguments. The completeness distribution determines
#' if there is a distribution of not random missingness of data along the dose. For this it is checked if half of a
#' features values (+/-1 value) pass the replicate completeness criteria and half do not pass it. In order to fall into
#' this category, the values that fulfill the completeness cutoff and the ones that do not fulfill it
#' need to be consecutive, meaning located next to each other based on their concentration values. Furthermore,
#' the values that do not pass the completeness cutoff need to be lower in intensity. Lastly, the difference
#' between the two groups is tested for statistical significance using a Welch's t-test and a
#' cutoff of p <= 0.1 (we want to mainly discard curves that falsely fit the other criteria but that
#' have clearly non-significant differences in mean). This allows curves to be considered that have
#' missing values in half of their observations due to a decrease in intensity. It can be thought
#' of as conditions that are missing not at random (MNAR). It is often the case that those entities
#' do not have a significant p-value since half of their conditions are not considered due to data
#' missingness. The ANOVA test is performed on the features by concentration. If it is significant it is
#' likely that there is some response. However, this test would also be significant even if there is one
#' outlier concentration so it should only be used only in combination with other cutoffs to determine
#' if a feature is significant. The `passed_filter` column is `TRUE` for all the
#' features that pass the above mentioned criteria and that have a correlation greater than the cutoff
#' (default is 0.8) and the adjusted ANOVA p-value below the cutoff (default is 0.05).
#'
#' The final list is ranked based on a score calculated on entities that pass the filter.
#' The score is the negative log10 of the adjusted ANOVA p-value scaled between 0 and 1 and the
#' correlation scaled between 0 and 1 summed up and divided by 2. Thus, the highest score an
#' entity can have is 1 with both the highest correlation and adjusted p-value. The rank is
#' corresponding to this score. Please note, that entities with MNAR conditions might have a
#' lower score due to the missing or non-significant ANOVA p-value. If no score could be calculated
#' the usual way these cases receive a score of 0. You should have a look at curves that are TRUE
#' for \code{dose_MNAR} in more detail.
#'
#' If the `"pre"` option is selected for the `filter` argument then the data is filtered for completeness
#' prior to curve fitting and the ANOVA test. Otherwise annotation is performed exactly as mentioned above.
#' We recommend the `"pre"` option because it leaves you with not only the likely hits of your treatment, but
#' also with rather high confidence true negative results. This is because the filtered data has a high
#' degree of completeness making it unlikely that a real dose-response curve is missed due to data missingness.
#'
#' Please note that in general, curves are only fitted if there are at least 5 conditions with data points present
#' to ensure that there is potential for a good curve fit. This is done independent of the selected filtering option.
#'
#' @param data a data frame that contains at least the input variables.
#' @param sample  a character column in the \code{data} data frame that contains the sample names.
#' @param grouping  a character column in the \code{data} data frame that contains the precursor,
#' peptide or protein identifiers.
#' @param response  a numeric column in the \code{data} data frame that contains the response
#' values, e.g. log2 transformed intensities.
#' @param dose  a numeric column in the \code{data} data frame that contains the dose values, e.g.
#' the treatment concentrations.
#' @param filter a character value that can either be `"pre"`, `"post"` or `"none"`. The data is
#' annotated for completeness, ANOVA significance and the completeness distribution along
#' the doses (`"pre"` and `"post"`). The combined output of this filtering step can be found in
#' the `passed_filter` column and depends on the cutoffs provided to the function. Note that this
#' is only an annotation and nothing is removed from the output. If `"pre"` is selected then, in
#' addition to the annotation, the data is filtered for completeness based on the condition completeness
#' prior to the curve fitting and ANOVA calculation and p-value adjustment. This has the benefit that less
#' curves need to be fitted and that the ANOVA p-value adjustment is done only on the relevant set of tests.
#' If `"none"` is selected the data will be neither annotated nor filtered.
#' @param replicate_completeness `r lifecycle::badge("deprecated")` please use `n_replicate_completeness` instead. 
#' A numeric value which similar to \code{completenss_MAR} of the
#' \code{assign_missingness} function sets a threshold for the completeness of data. In contrast
#' to \code{assign_missingness} it only determines the completeness for one condition and not the
#' comparison of two conditions. The threshold is used to calculate a minimal degree of data
#' completeness. The value provided to this argument has to be between 0 and 1, default is 0.7.
#' It is multiplied with the number of replicates and then adjusted downward. The resulting number
#' is the minimal number of observations that a condition needs to have to be considered "complete
#' enough" for the \code{condition_completeness} argument.
#' @param condition_completeness `r lifecycle::badge("deprecated")` please use `n_condition_completeness` instead. 
#' A numeric value which determines how many conditions need to at
#' least fulfill the "complete enough" criteria set with \code{replicate_completeness}. The
#' value provided to this argument has to be between 0 and 1, default is 0.5. It is multiplied with
#' the number of conditions and then adjusted downward. The resulting number is the minimal number
#' of conditions that need to fulfill the \code{replicate_completeness} argument for a peptide to
#' pass the filtering.
#' @param n_replicate_completeness a numeric value that defines the minimal number of observations that a 
#' condition (concentration) needs to have to be considered "complete enough" for the `n_condition_completeness` 
#' argument. E.g. if each concentration has 4 replicates this argument could be set to 3 to allow for one 
#' replicate to be missing for the completeness criteria.
#' @param n_condition_completeness a numeric value that defines the minimal number
#' of conditions that need to fulfill the `n_replicate_completeness` argument for a feature to
#' pass the filtering. E.g. if an experiment has 12 concentrations, this argument could be set
#' to 6 to define that at least 6 of 12 concentrations need to make the replicate completeness cutoff.
#' @param anova_cutoff a numeric value that specifies the ANOVA adjusted p-value cutoff used for
#' data filtering. Any fits with an adjusted ANOVA p-value bellow the cutoff will be considered
#' for scoring.
#' @param correlation_cutoff a numeric value that specifies the correlation cutoff used for data
#' filtering. Any fits with a correlation above the cutoff will be considered for scoring.
#' @param log_logarithmic a logical value that indicates if a logarithmic or log-logarithmic model
#' is fitted. If response values form a symmetric curve for non-log transformed dose values, a
#' logarithmic model instead of a log-logarithmic model should be used. Usually biological dose
#' response data has a log-logarithmic distribution, which is the reason this is the default.
#' Log-logarithmic models are symmetric if dose values are log transformed.
#' @param retain_columns a vector that specifies columns that should be retained from the input
#' data frame. Default is not retaining additional columns \code{retain_columns = NULL}. Specific
#' columns can be retained by providing their names (not in quotations marks, just like other
#' column names, but in a vector).
#' @param n_cores optional, a numeric value that specifies the number of cores used if workers
#' are set up manually.
#'
#' @return A data frame is returned that contains correlations of predicted to measured values as
#' a measure of the goodness of the curve fit, an associated p-value and the four parameters of
#' the model for each group. Furthermore, input data for plots is returned in the columns \code{plot_curve}
#' (curve and confidence interval) and \code{plot_points} (measured points).
#'
#' @import dplyr
#' @importFrom stats p.adjust
#' @importFrom rlang .data as_name enquo
#' @importFrom magrittr %>%
#' @export
#'
#' @examples
#' \dontrun{
#' # Load libraries
#' library(dplyr)
#'
#' set.seed(123) # Makes example reproducible
#'
#' # Create example data
#' data <- create_synthetic_data(
#'   n_proteins = 2,
#'   frac_change = 1,
#'   n_replicates = 3,
#'   n_conditions = 8,
#'   method = "dose_response",
#'   concentrations = c(0, 1, 10, 50, 100, 500, 1000, 5000),
#'   additional_metadata = FALSE
#' )
#'
#' # Perform dose response curve fit
#' drc_fit <- parallel_fit_drc_4p(
#'   data = data,
#'   sample = sample,
#'   grouping = peptide,
#'   response = peptide_intensity_missing,
#'   dose = concentration,
#'   n_replicate_completeness = 2,
#'   n_condition_completeness = 5,
#'   retain_columns = c(protein, change_peptide)
#' )
#'
#' glimpse(drc_fit)
#'
#' head(drc_fit, n = 10)
#' }
parallel_fit_drc_4p <- function(data,
                                sample,
                                grouping,
                                response,
                                dose,
                                filter = "post",
                                replicate_completeness = 0.7,
                                condition_completeness = 0.5,
                                n_replicate_completeness = NULL,
                                n_condition_completeness = NULL,
                                anova_cutoff = 0.05,
                                correlation_cutoff = 0.8,
                                log_logarithmic = TRUE,
                                retain_columns = NULL,
                                n_cores = NULL) {
  dependency_test <- c(
    furrr = !requireNamespace("furrr", quietly = TRUE),
    future = !requireNamespace("future", quietly = TRUE),
    parallel = !requireNamespace("parallel", quietly = TRUE),
    drc = !requireNamespace("drc", quietly = TRUE)
  )
  if (any(dependency_test)) {
    dependency_name <- names(dependency_test[dependency_test == TRUE])
    if (length(dependency_name) == 1) {
      message("Package \"",
        paste(dependency_name),
        "\" is needed for this function to work. Please install it.",
        call. = FALSE
      )
      return(invisible(NULL))
    } else {
      message("Packages \"",
        paste(dependency_name, collapse = "\" and \""),
        "\" are needed for this function to work. Please install them.",
        call. = FALSE
      )
      return(invisible(NULL))
    }
  }
  . <- NULL
  terminate <- FALSE

  if (missing(n_cores)) {
    message("Setting up workers ... ", appendLF = FALSE)
    terminate <- TRUE
    n_cores <- parallel::detectCores()
    n_cores <- ifelse(n_cores > 12, 12, n_cores)
    oplan <- future::plan(future::multisession, workers = n_cores)
    on.exit(future::plan(oplan), add = TRUE)
    message("DONE", appendLF = TRUE)
  }

  pieces <- rep(
    1:n_cores,
    round(
      length(
        unique(dplyr::pull(data, {{ grouping }}))
      ) / n_cores
    ) + 1
  )[1:length(unique(dplyr::pull(data, {{ grouping }})))]

  pieces_mapping <- tibble::tibble({{ grouping }} := unique(dplyr::pull(data, {{ grouping }})), piece = pieces)

  input <- data %>%
    dplyr::left_join(pieces_mapping, by = rlang::as_name(rlang::enquo(grouping))) %>%
    split(dplyr::pull(., .data$piece))

  message("Performing model fit (this may take a while) ... ", appendLF = FALSE)

  result <- furrr::future_map_dfr(
    .x = input,
    .f = ~ protti::fit_drc_4p(.x,
      sample = {{ sample }},
      grouping = {{ grouping }},
      response = {{ response }},
      dose = {{ dose }},
      filter = filter,
      replicate_completeness = replicate_completeness,
      condition_completeness = condition_completeness,
      n_replicate_completeness = n_replicate_completeness,
      n_condition_completeness = n_condition_completeness,
      anova_cutoff = anova_cutoff,
      correlation_cutoff = correlation_cutoff,
      log_logarithmic = log_logarithmic,
      retain_columns = {{ retain_columns }},
      include_models = FALSE
    ),
    .options = furrr::furrr_options(globals = FALSE)
  )

  message("DONE", appendLF = TRUE)

  if (terminate == TRUE) {
    future::plan(future::sequential)
  }

  result <- result %>%
    dplyr::arrange(desc(.data$correlation))

  if (filter != "none") {
    result <- result %>%
      dplyr::mutate(anova_adj_pval = stats::p.adjust(.data$anova_pval, method = "BH")) %>%
      dplyr::mutate(anova_significant = ifelse(.data$anova_adj_pval > anova_cutoff | is.na(.data$anova_adj_pval),
        FALSE,
        TRUE
      )) %>%
      dplyr::mutate(passed_filter = (((.data$enough_conditions == TRUE &
        .data$anova_significant == TRUE) |
        .data$dose_MNAR == TRUE) &
        .data$correlation >= correlation_cutoff &
        .data$min_model < .data$max_model) |
          .data$dose_MNAR) %>%
      dplyr::group_by(.data$passed_filter) %>%
      dplyr::mutate(score = ifelse(.data$passed_filter & .data$anova_significant & .data$correlation >= correlation_cutoff,
        (scale_protti(-log10(.data$anova_pval), method = "01") + scale_protti(.data$correlation, method = "01")) / 2,
        NA
      )) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(score = ifelse(is.na(.data$score) & .data$passed_filter == TRUE, 0, .data$score)) %>%
      dplyr::arrange(dplyr::desc(.data$correlation)) %>%
      dplyr::arrange(dplyr::desc(.data$score)) %>%
      dplyr::select(-.data$rank) %>%
      tibble::rownames_to_column(var = "rank") %>%
      dplyr::mutate(rank = ifelse(!is.na(.data$score), as.numeric(.data$rank), NA))
  }

  return(result)
}
