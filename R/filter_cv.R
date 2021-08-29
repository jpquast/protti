#' Data filtering based on coefficients of variation (CV)
#'
#' Filters the input data based on precursor, peptide or protein intensity coefficients of variation.
#' The function should be used to ensure that only robust measurements and quantifications are used for
#' data analysis. It is advised to use the function after inspection of raw values (quality control)
#' and median normalisation. Generally, the function calculates CVs of each peptide, precursor or
#' protein for each condition and removes peptides, precursors or proteins that have a CV above
#' the cutoff in less than the (user-defined) required number of conditions. Since the user-defined
#' cutoff is fixed and does not depend on the number of conditions that have detected values, the
#' function might bias for data completeness.
#'
#' @param data a data frame that contains at least the input variables.
#' @param grouping a character column in the \code{data} data frame that contains the grouping
#' variable that can be either precursors, peptides or proteins.
#' @param condition a character or numeric column in the \code{data} data frame that contains
#' information on the sample condition.
#' @param log2_intensity a numeric column in the \code{data} data frame that contains log2
#' transformed intensities.
#' @param cv_limit optional, a numeric value that specifies the CV cutoff that will be applied.
#' Default is 0.25.
#' @param min_conditions a numeric value that specifies the minimum number of conditions for
#' which grouping CVs should be below the cutoff.
#' @param silent a logical value that specifies if a message with the number of filtered out
#' conditions should be returned. Default is FALSE.
#'
#' @return The CV filtered data frame.
#' @importFrom magrittr %>%
#' @importFrom dplyr group_by mutate filter ungroup pull
#' @export
#'
#' @examples
#' set.seed(123) # Makes example reproducible
#'
#' # Create synthetic data
#' data <- create_synthetic_data(
#'   n_proteins = 50,
#'   frac_change = 0.05,
#'   n_replicates = 3,
#'   n_conditions = 2,
#'   method = "effect_random",
#'   additional_metadata = FALSE
#' )
#'
#' # Filter coefficients of variation
#' data_filtered <- filter_cv(
#'   data = data,
#'   grouping = peptide,
#'   condition = condition,
#'   log2_intensity = peptide_intensity_missing,
#'   cv_limit = 0.25,
#'   min_conditions = 2
#' )
filter_cv <- function(data,
                      grouping,
                      condition,
                      log2_intensity,
                      cv_limit = 0.25,
                      min_conditions,
                      silent = FALSE) {
  if (max(dplyr::pull(data, {{ log2_intensity }}), na.rm = TRUE) > 50) {
    stop(strwrap("Please transform your data. The function requires your data to be log2 transformed.",
      prefix = "\n", initial = ""
    ))
  }
  n_groups_start <- length(unique(dplyr::pull(data, {{ grouping }})))

  peptide_list <- data %>%
    dplyr::group_by({{ grouping }}, {{ condition }}) %>%
    dplyr::mutate(cv_small = (sd(2^{{ log2_intensity }}, na.rm = TRUE) /
      mean(2^{{ log2_intensity }}, na.rm = TRUE)
    < {{ cv_limit }}) / dplyr::n()) %>%
    dplyr::group_by({{ grouping }}) %>%
    dplyr::mutate(cv_count = sum(.data$cv_small, na.rm = TRUE)) %>%
    dplyr::filter(.data$cv_count >= {{ min_conditions }}) %>%
    dplyr::select(-c(.data$cv_small, .data$cv_count)) %>%
    dplyr::ungroup()

  n_groups_end <- length(unique(dplyr::pull(peptide_list, {{ grouping }})))

  if (silent == FALSE) {
    message(
      n_groups_start - n_groups_end,
      " groups of ",
      n_groups_start,
      " were filtered out. ",
      round(n_groups_end / n_groups_start, digits = 2) * 100,
      "% of data remains."
    )
  }

  peptide_list
}
