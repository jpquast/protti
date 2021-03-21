#' Protein coverage distribution
#'
#' Plots the distribution of protein coverages in a histogram.
#'
#' @param data A data frame containing at least the input variables.
#' @param protein_identifier Column in the data frame containing protein identifiers.
#' @param coverage Column in the data frame containing protein coverage in percent. This information can be obtained using
#' the \code{\link{sequence_coverage}} function.
#' @param sample Column in the data frame containing sample names. Please only provide this argument if you want to
#' facet the distribution plot by sample otherwise do not provide this argument.
#' @param interactive A logical, if TRUE the plot is interactive (default is TRUE).
#'
#' @return A protein coverage histogram with 5 percent binwidth. The vertical dotted line indicates the median.
#' @import dplyr
#' @import ggplot2
#' @importFrom tidyr drop_na
#' @importFrom magrittr %>%
#' @importFrom plotly ggplotly
#' @importFrom stringr str_sort
#' @export
#'
#' @examples
#' \dontrun{
#' qc_sequence_coverage(
#'   data,
#'   protein_identifier = pg_protein_accessions,
#'   coverage = coverage
#' )
#' }
#'
#' @seealso \code{\link{sequence_coverage}}
qc_sequence_coverage <- function(data, protein_identifier, coverage, sample = NULL, interactive = FALSE) {
  result <- data %>%
    dplyr::distinct({{ protein_identifier }}, {{ coverage }}, {{ sample }}) %>%
    tidyr::drop_na({{ coverage }})

  if (!missing(sample)) {
    result <- result %>%
      dplyr::mutate({{ sample }} := factor({{ sample }}, levels = unique(stringr::str_sort({{ sample }}, numeric = TRUE))))
  }

  plot <- result %>%
    ggplot2::ggplot(ggplot2::aes({{ coverage }})) +
    ggplot2::geom_histogram(binwidth = 5, col = "black", fill = "#5680C1", boundary = 0, size = 1) +
    ggplot2::geom_vline(xintercept = stats::median(dplyr::pull(result, {{ coverage }}), na.rm = TRUE), linetype = "dashed") +
    ggplot2::labs(title = "Protein coverage distribution", x = "Coverage [%]", y = "Number of proteins") +
    ggplot2::scale_x_continuous(breaks = seq(from = 0, to = 100, by = 10)) +
    {
      if (!missing(sample)) ggplot2::facet_wrap(rlang::new_formula(NULL, rlang::enquo(sample)), scales = "free", ncol = 4)
    } +
    ggplot2::theme_bw() +
    theme(
      plot.title = ggplot2::element_text(size = 20),
      axis.title.x = ggplot2::element_text(size = 15),
      axis.text.y = ggplot2::element_text(size = 15),
      axis.text.x = ggplot2::element_text(size = 12),
      axis.title.y = ggplot2::element_text(size = 15),
      strip.text = ggplot2::element_text(size = 15),
      panel.border = ggplot2::element_rect(fill = NA),
      strip.background = element_rect(fill = "white")
    )

  if (interactive == FALSE) {
    return(plot)
  }
  plotly::ggplotly(plot)
}
