#' Perform ANOVA
#'
#' Performs an ANOVA statistical test
#'
#' @param data a data frame containing at least the input variables.
#' @param grouping a character column in the \code{data} data frame that contains precursor or
#' peptide identifiers.
#' @param condition a character or numeric column in the \code{data} data frame that contains the
#' conditions.
#' @param mean_ratio a numeric column in the \code{data} data frame that contains mean intensities
#' or mean intensity ratios.
#' @param sd a numeric column in the \code{data} data frame that contains the standard deviation
#' corresponding to the mean.
#' @param n a numeric column in the \code{data} data frame that contains the number of replicates
#' for which the corresponding mean was calculated.
#'
#' @return a data frame that contains the within group error (\code{ms_group}) and the between
#' group error (\code{ms_error}), f statistic and p-values.
#' @import dplyr
#' @export
#'
#' @examples
#' data <- data.frame(
#'   precursor = c("A", "A", "A", "B", "B", "B"),
#'   condition = c("C1", "C2", "C3", "C1", "C2", "C3"),
#'   mean = c(10, 12, 20, 11, 12, 8),
#'   sd = c(2, 1, 1.5, 1, 2, 4),
#'   n = c(4, 4, 4, 4, 4, 4)
#' )
#'
#' anova_protti(
#'   data,
#'   grouping = precursor,
#'   condition = condition,
#'   mean = mean,
#'   sd = sd,
#'   n = n
#' )
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
