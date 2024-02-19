#' Protein coverage distribution
#'
#' Plots the distribution of protein coverages in a histogram.
#'
#' @param data a data frame that contains at least the input variables.
#' @param protein_identifier a character column in the \code{data} data frame that contains protein
#' identifiers.
#' @param coverage a numeric column in the \code{data} data frame that contains protein coverage
#' in percent. This information can be obtained using the \code{\link{sequence_coverage}} function.
#' @param sample optional, a character or factor column in the \code{data} data frame that contains sample names.
#' Please only provide this argument if you want to facet the distribution plot by sample
#' otherwise do not provide this argument.
#' @param interactive a logical value that specifies whether the plot should be interactive
#' (default is FALSE).
#'
#' @return A protein coverage histogram with 5 percent binwidth. The vertical dotted line
#' indicates the median.
#' @import dplyr
#' @import ggplot2
#' @importFrom tidyr drop_na
#' @importFrom magrittr %>%
#' @importFrom plotly ggplotly
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
#' # Plot sequence coverage
#' qc_sequence_coverage(
#'   data = data,
#'   protein_identifier = protein,
#'   coverage = coverage
#' )
#' @seealso \code{\link{sequence_coverage}}
qc_sequence_coverage <- function(data,
                                 protein_identifier,
                                 coverage,
                                 sample = NULL,
                                 interactive = FALSE) {

  # Validate inputs
  if (!all(c(protein_identifier, coverage) %in% colnames(data))) {
    stop("Column names for protein_identifier and coverage must exist in the dataset.")
  }
  
  result <- data %>%
    dplyr::distinct({{ protein_identifier }}, {{ coverage }}, {{ sample }}) %>%
    tidyr::drop_na({{ coverage }}) %>%
    dplyr::group_by({{ sample }}) %>%
    dplyr::mutate(median_coverage = median({{ coverage }})) %>%
    dplyr::ungroup()

  if (!missing(sample) && is(dplyr::pull(result, {{ sample }}), "character")) {
    result <- result %>%
      dplyr::mutate({{ sample }} := factor({{ sample }},
        levels = unique(stringr::str_sort({{ sample }}, numeric = TRUE))
      ))
  }

  plot <- result %>%
    ggplot2::ggplot(ggplot2::aes({{ coverage }})) +
    ggplot2::geom_histogram(
      binwidth = 5,
      col = "black",
      fill = "#5680C1",
      boundary = 0,
      size = 1
    ) +
    ggplot2::geom_vline(data = result %>% dplyr::distinct(.data$median_coverage, {{ sample }}),
                        mapping = aes(xintercept = median_coverage),
                        linewidth = 1,
                        linetype = "dashed") +
    ggplot2::labs(
      title = "Protein coverage distribution",
      x = "Coverage [%]",
      y = "Number of proteins"
    ) +
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
