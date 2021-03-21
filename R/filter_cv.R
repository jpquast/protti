#' Data filtering based on coefficients of variation (CV)
#'
#' Filters the input data based on precursor, peptide or protein intensity coefficients of variation.
#' The function should be used to ensure that only robust measurements and quantifications are used for
#' data analysis. It is advised to use the function after inspection of raw values (quality control) and
#' median normalisation. Generally, the function calculates CVs of each peptide, precursor or protein for
#' each condition and removes peptides, precursors or proteins that have a CV above the cutoff in less
#' than the (user-defined) required number of conditions. Since the user-defined cutoff is fixed and does not
#' depend on the number of conditions that have detected values, the function might bias for data completeness.
#'
#' @param data Data frame containing at least the input variables.
#' @param grouping Column in the data frame containing the grouping variable that can be either precursors, peptides or proteins.
#' @param condition Column in the data frame containing information on the sample condition.
#' @param log2_intensity Column in the data frame containing log2 transformed intensities.
#' @param cv_limit Optional argument specifying the CV cutoff that will be applied. Default is 0.25.
#' @param min_conditions The minimum number of conditions for which grouping CVs should be below the cutoff.
#' @param silent Logical argument specifiying if a message with the number of filtered out conditions should be returned.
#'
#' @return The CV filtered data frame.
#' @importFrom magrittr %>%
#' @importFrom dplyr group_by mutate filter ungroup pull
#' @export
#'
#' @examples
#' \dontrun{
#' filter_cv(
#'   data,
#'   grouping = eg_precursor_id,
#'   condition = r_condition,
#'   log2_intensity = normalised_intensity_log2,
#'   cv_limit = 0.25,
#'   min_conditions = 5
#' )
#' }
filter_cv <- function(data, grouping, condition, log2_intensity, cv_limit = 0.25, min_conditions, silent = FALSE) {
  if (max(dplyr::pull(data, {{ log2_intensity }}), na.rm = TRUE) > 50) {
    stop("Please transform your data. The function requires your data to be log2 transformed.")
  }
  n_groups_start <- length(unique(dplyr::pull(data, {{ grouping }})))

  peptide_list <- data %>%
    dplyr::group_by({{ grouping }}, {{ condition }}) %>%
    dplyr::mutate(cv_small = (sd(2^{{ log2_intensity }}, na.rm = TRUE) / mean(2^{{ log2_intensity }}, na.rm = TRUE) < {{ cv_limit }}) / dplyr::n()) %>%
    dplyr::group_by({{ grouping }}) %>%
    dplyr::mutate(cv_count = sum(.data$cv_small, na.rm = TRUE)) %>%
    dplyr::filter(.data$cv_count >= {{ min_conditions }}) %>%
    dplyr::select(-c(.data$cv_small, .data$cv_count)) %>%
    dplyr::ungroup()

  n_groups_end <- length(unique(dplyr::pull(peptide_list, {{ grouping }})))

  if (silent == FALSE) {
    message(n_groups_start - n_groups_end, " groups of ", n_groups_start, " were filtered out. ", round(n_groups_end / n_groups_start, digits = 2) * 100, "% of data remains.")
  }

  peptide_list
}
