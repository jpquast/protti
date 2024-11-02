#' Median run intensities
#'
#' Median intensities per run are returned either as a plot or a table.
#'
#' @param data a data frame that contains at least the input variables.
#' @param sample a character or factor column in the \code{data} data frame that contains the sample name.
#' @param grouping a character column in the \code{data} data frame that contains either precursor or
#' peptide identifiers.
#' @param intensity a numeric column in the \code{data} data frame that contains intensity values.
#' The intensity should be ideally log2 transformed, but also non-transformed values can be used.
#' @param plot a logical value that indicates whether the result should be plotted.
#' @param interactive a logical value that specifies whether the plot should be interactive
#' (default is FALSE).
#'
#' @return A plot that displays median intensity over all samples. If \code{plot = FALSE} a data
#' frame containing median intensities is returned.
#' @import dplyr
#' @import ggplot2
#' @importFrom plotly ggplotly
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @importFrom stringr str_sort
#' @export
#'
#' @examples
#' set.seed(123) # Makes example reproducible
#'
#' # Create example data
#' data <- create_synthetic_data(
#'   n_proteins = 100,
#'   frac_change = 0.05,
#'   n_replicates = 3,
#'   n_conditions = 2,
#'   method = "effect_random"
#' )
#'
#' # Calculate median intensities
#' qc_median_intensities(
#'   data = data,
#'   sample = sample,
#'   grouping = peptide,
#'   intensity = peptide_intensity_missing,
#'   plot = FALSE
#' )
#'
#' # Plot median intensities
#' qc_median_intensities(
#'   data = data,
#'   sample = sample,
#'   grouping = peptide,
#'   intensity = peptide_intensity_missing,
#'   plot = TRUE
#' )
qc_median_intensities <- function(data,
                                  sample,
                                  grouping,
                                  intensity,
                                  plot = TRUE,
                                  interactive = FALSE) {
  table <- data %>%
    dplyr::distinct({{ sample }}, {{ grouping }}, {{ intensity }}) %>%
    dplyr::group_by({{ sample }}) %>%
    dplyr::summarize(
      median_intensity = stats::median({{ intensity }}, na.rm = TRUE),
      .groups = "drop"
    )

  if (plot == FALSE) {
    return(table)
  }

  if (is(dplyr::pull(table, {{ sample }}), "character")) {
    table <- table %>%
      mutate({{ sample }} := factor({{ sample }},
        levels = unique(stringr::str_sort({{ sample }}, numeric = TRUE))
      ))
  }

  plot <- table %>%
    ggplot2::ggplot(ggplot2::aes({{ sample }}, .data$median_intensity, group = 1)) +
    ggplot2::geom_line(linewidth = 1) +
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
