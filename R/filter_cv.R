#' Data filtering based on coefficients of variation (CV)
#'
#' Filters the input data based on precursor, peptide or protein intensity coefficients of variation.
#'
#' @param data Data frame containing at least the input variables.
#' @param grouping Column in the data frame containing the grouping variable that can be either precursors, peptides or proteins.
#' @param condition Column in the data frame containing information on the sample condition.
#' @param log2_intensity Column in the data frame containing log2 transformed intensities.
#' @param cv_limit Optional argument specifying the CV cutoff that will be applied. Default is 0.25.
#' @param min_conditions The minimum number of conditions for which grouping CVs should be below the cutoff.
#'
#' @return The CV filtered data frame.
#' @import ggplot2
#' @importFrom magrittr %>%
#' @importFrom dplyr group_by summarise pull filter
#' @importFrom utils data
#' @export
#'
#' @examples
#' \dontrun{
#' filter_cv(
#'   data,
#'   grouping = eg_precursor_id,
#'   intensity = log2_intensity,
#'   condition = r_condition,
#'   log2_intensity = normalised_intensity_log2,
#'   cv_limit = 0.25,
#'   min_conditions = 5
#' )
#' }
filter_cv <- function (data, grouping, condition, log2_intensity, cv_limit = 0.25, min_conditions)
{

  if(max(dplyr::pull(data, {{log2_intensity}}), na.rm = TRUE) > 50)
  {
    stop("Please transform your data. The function requires your data to be log2 transformed.")
  }

  peptide_list <- data %>%
    dplyr::group_by({{grouping}}, {{condition}}) %>%
    dplyr::summarise(cv_count = sum( sd(2^{{log2_intensity}}, na.rm = TRUE) / mean(2^{{log2_intensity}}, na.rm = TRUE) < {{cv_limit}} ), .groups = 'drop') %>%
    dplyr::filter(.data$cv_count >= 1) %>%
    dplyr::group_by({{grouping}}) %>%
    dplyr::summarise(cv_count2 = sum(.data$cv_count), .groups = 'drop') %>%
    dplyr::filter(.data$cv_count2 >= {{min_conditions}}) %>%
    dplyr::pull({{grouping}})

  dplyr::filter(data, {{grouping}} %in% peptide_list)
}
