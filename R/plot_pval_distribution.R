#' Plot histogram of p-value distribution
#'
#' Plots the distribution of p-values derived from any statistical test as a histogram.
#'
#' @param data a data frame containing at least grouping identifiers (precursor, peptide or protein) and p-values derived from any
#' statistical test.
#' @param grouping the column in the data frame containing either precursor, peptide or protein identifiers.
#' @param pval the column in the data frame containing p-values.
#'
#' @return A histogram or boxplot that shows the intensity distribution over all samples or by sample.
#' @import ggplot2
#' @importFrom magrittr %>%
#' @importFrom dplyr distinct
#' @importFrom tidyr drop_na
#' @export
#'
#' @examples
#' \dontrun{
#' plot_pval_distribution(
#'   data,
#'   grouping = eg_precursor_id,
#'   pval = pval
#' )
#' }
plot_pval_distribution <- function(data, grouping, pval) {
  input <- data %>%
    dplyr::distinct({{ grouping }}, {{ pval }}) %>%
    tidyr::drop_na()

  plot <- input %>%
    ggplot2::ggplot(ggplot2::aes(x = {{ pval }})) +
    ggplot2::geom_histogram(binwidth = 0.05, color = "black", fill = "#5680C1", size = 1) +
    ggplot2::labs(title = "P-Value Distribution", x = "P-Value", y = "Frequency") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 20),
      axis.title.x = ggplot2::element_text(size = 15),
      axis.text.y = ggplot2::element_text(size = 15),
      axis.text.x = ggplot2::element_text(size = 15),
      axis.title.y = ggplot2::element_text(size = 15)
    )
  plot
}
