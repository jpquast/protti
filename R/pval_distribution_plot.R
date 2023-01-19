#' Plot histogram of p-value distribution
#'
#' `r lifecycle::badge('deprecated')`
#' This function was deprecated due to its name changing to `pval_distribution_plot()`.
#'
#' @return A histogram plot that shows the p-value distribution.
#' @keywords internal
#' @export
plot_pval_distribution <- function(...) {
  # This function has been renamed and is therefore deprecated.
  lifecycle::deprecate_warn("0.2.0",
    "plot_pval_distribution()",
    "pval_distribution_plot()",
    details = "This function has been renamed."
  )

  pval_distribution_plot(...)
}
#' Plot histogram of p-value distribution
#'
#' Plots the distribution of p-values derived from any statistical test as a histogram.
#'
#' @param data a data frame that contains at least grouping identifiers (precursor, peptide or
#' protein) and p-values derived from any statistical test.
#' @param grouping a character column in the \code{data} data frame that contains either precursor,
#' peptide or protein identifiers. For each entry in this column there should be one unique p-value.
#' That means the statistical test that created the p-value should have been performed on the
#' level of the content of this column.
#' @param pval a numeric column in the \code{data} data frame that contains p-values.
#' @param facet_by optional, a character column that contains information by which the data should
#' be faceted into multiple plots.
#'
#' @return A histogram plot that shows the p-value distribution.
#' @import ggplot2
#' @importFrom magrittr %>%
#' @importFrom dplyr distinct
#' @importFrom tidyr drop_na
#' @export
#'
#' @examples
#' set.seed(123) # Makes example reproducible
#'
#' # Create example data
#' data <- data.frame(
#'   peptide = paste0("peptide", 1:1000),
#'   pval = runif(n = 1000)
#' )
#'
#' # Plot p-values
#' pval_distribution_plot(
#'   data = data,
#'   grouping = peptide,
#'   pval = pval
#' )
pval_distribution_plot <- function(data, grouping, pval, facet_by = NULL) {
  input <- data %>%
    dplyr::distinct({{ grouping }}, {{ pval }}, {{ facet_by }}) %>%
    tidyr::drop_na()

  plot <- input %>%
    ggplot2::ggplot(ggplot2::aes(x = {{ pval }})) +
    ggplot2::geom_histogram(
      binwidth = 0.05,
      boundary = 0,
      color = "black",
      fill = "#5680C1",
      size = 1
    ) +
    ggplot2::labs(title = "P-Value Distribution", x = "P-Value", y = "Frequency") +
    {
      if (!missing(facet_by)) {
        ggplot2::facet_wrap(rlang::new_formula(NULL, rlang::enquo(facet_by)),
          scales = "fixed"
        )
      }
    } +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 20),
      axis.title.x = ggplot2::element_text(size = 15),
      axis.text.y = ggplot2::element_text(size = 15),
      axis.text.x = ggplot2::element_text(size = 15),
      axis.title.y = ggplot2::element_text(size = 15),
      strip.text = ggplot2::element_text(size = 15),
      strip.background = element_blank()
    )
  plot
}
