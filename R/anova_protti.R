#' Perform ANOVA
#'
#' Performs an ANOVA statistical test
#'
#' @param data A data frame containing at least the input variables.
#' @param grouping The column in the data frame containing precursor or peptide identifiers.
#' @param condition The column in the data frame containing the conditions.
#' @param mean_ratio The column in the data frame containing mean intensities or mean intensity ratios.
#' @param sd The column in the data frame containing the standard deviation corresponding to the mean.
#' @param n The column in the data frame containing the number of replicates for which the corresponding mean was calculated.
#'
#' @return A data frame that contains the within group error (\code{ms_group}) and the between group error (\code{ms_error}), f statistic and p-values.
#' @import dplyr
#' @export
#'
#' @examples
#' \dontrun{
#' anova_protti(
#'   data,
#'   grouping = eg_precursor_id,
#'   condition = r_condition,
#'   mean = mean,
#'   sd = sd,
#'   n = n
#' )
#' }
anova_protti <- function(data, grouping, condition, mean_ratio, sd, n) {
  result <- data %>%
    dplyr::distinct({{ grouping }}, {{ condition }}, {{ mean_ratio }}, {{ sd }}, {{ n }}) %>%
    dplyr::group_by({{ grouping }}) %>%
    dplyr::filter({{ n }} != 0) %>%
    dplyr::mutate(n_groups = dplyr::n_distinct(!!ensym(condition))) %>%
    dplyr::mutate(grand_mean = mean({{ mean_ratio }})) %>%
    dplyr::mutate(total_n = sum({{ n }})) %>%
    dplyr::mutate(ms_group = sum(({{ mean_ratio }} - .data$grand_mean)^2 * {{ n }}) / (.data$n_groups - 1)) %>%
    dplyr::mutate(ms_error = sum({{ sd }}^2 * ({{ n }} - 1)) / (.data$total_n - .data$n_groups)) %>%
    dplyr::mutate(f = .data$ms_group / .data$ms_error) %>%
    dplyr::mutate(pval = stats::pf(.data$f, .data$n_groups - 1, .data$total_n - .data$n_groups, lower.tail = FALSE)) %>%
    dplyr::distinct({{ grouping }}, .data$ms_group, .data$ms_error, .data$f, .data$pval) %>%
    dplyr::ungroup()

  result
}
