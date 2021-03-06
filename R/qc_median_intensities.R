#' Median run intensities
#'
#' Median intensities per run are returned either as a plot or a table.
#'
#' @param data a data frame containing at least the input variables.
#' @param sample the name of the column containing the sample names.
#' @param grouping the name of the column containing precursor or peptide identifiers.
#' @param intensity the name of the column containing intensity values. The intensity should be ideally log2 transformed, but also non-transformed values can be used.
#' @param plot logical, if TRUE a plot is returned. If FALSE a table is returned.
#' @param interactive logical, if TRUE the plot is interactive using plotly.
#'
#' @return A plot that displays median intensity over all samples. If \code{plot = FALSE} a data frame containing median intensities is returned.
#' @import dplyr
#' @import ggplot2
#' @importFrom plotly ggplotly
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @importFrom stringr str_sort
#' @export
#'
#' @examples
#' \dontrun{
#' qc_median_intensities(
#'   data,
#'   sample = r_file_name,
#'   grouping = eg_precursor_id,
#'   intensity = log2_intensity
#' )
#' }
qc_median_intensities <- function(data, sample, grouping, intensity, plot = TRUE, interactive = FALSE) {
  table <- data %>%
    dplyr::distinct({{ sample }}, {{ grouping }}, {{ intensity }}) %>%
    dplyr::group_by({{ sample }}) %>%
    dplyr::summarize(median_intensity = stats::median({{ intensity }}, na.rm = TRUE), .groups = "drop")

  if (plot == FALSE) {
    return(table)
  }

  plot <- table %>%
    mutate({{ sample }} := factor({{ sample }}, levels = unique(stringr::str_sort({{ sample }}, numeric = TRUE)))) %>%
    ggplot2::ggplot(ggplot2::aes({{ sample }}, .data$median_intensity, group = 1)) +
    ggplot2::geom_line(size = 1) +
    ggplot2::labs(title = "Medians of run intensities", x = "", y = "Intensity") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 20),
      axis.title.x = ggplot2::element_text(size = 15),
      axis.text.y = ggplot2::element_text(size = 15),
      axis.text.x = ggplot2::element_text(size = 12, angle = 75, hjust = 1),
      axis.title.y = ggplot2::element_text(size = 15)
    )

  if (interactive == FALSE) {
    return(plot)
  }

  suppressWarnings(plotly::ggplotly(plot))
}
