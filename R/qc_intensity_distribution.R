#' Check intensity distribution per sample and overall
#'
#' Plots the overall or sample-wise distribution of all peptide intensities as a boxplot or
#' histogram.
#'
#' @param data a data frame that contains at least sample names, grouping identifiers (precursor,
#' peptide or protein) and log2 transformed intensities for each grouping identifier.
#' @param sample an optional character or factor column in the \code{data} data frame that contains the
#' sample name. If the sample column is of type factor, the ordering is based on the factor
#' levels. NOTE: If the overall distribution should be returned please do not provide the name of the
#' sample column.
#' @param grouping a character column in the \code{data} data frame that contains the grouping
#' variables (e.g. peptides, precursors or proteins).
#' @param intensity_log2 a numeric column in the \code{data} data frame that contains the log2
#' transformed intensities of each grouping identifier sample combination.
#' @param plot_style a character value that indicates the plot type. This can be either
#' "histogram", "boxplot" or "violin". Plot style "boxplot" and "violin" can only be used if a
#' sample column is provided.
#'
#' @return A histogram or boxplot that shows the intensity distribution over all samples or by
#' sample.
#' @import ggplot2
#' @importFrom dplyr distinct mutate
#' @importFrom magrittr %>%
#' @importFrom tidyr drop_na
#' @importFrom stringr str_sort
#' @importFrom rlang new_formula enquo
#' @importFrom methods is
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
#' # Plot intensity distribution
#' # The plot style can be changed
#' qc_intensity_distribution(
#'   data = data,
#'   sample = sample,
#'   grouping = peptide,
#'   intensity_log2 = peptide_intensity_missing,
#'   plot_style = "boxplot"
#' )
qc_intensity_distribution <- function(data,
                                      sample = NULL,
                                      grouping,
                                      intensity_log2,
                                      plot_style) {
  if (missing(plot_style)) stop("Please provide a plot type in the plot_style argument!")
  input <- data %>%
    dplyr::distinct({{ sample }}, {{ grouping }}, {{ intensity_log2 }}) %>%
    tidyr::drop_na({{ intensity_log2 }})

  if (!missing(sample) && is(dplyr::pull(input, {{ sample }}), "character")) {
    input <- input %>%
      dplyr::mutate({{ sample }} := factor({{ sample }},
        levels = unique(stringr::str_sort({{ sample }}, numeric = TRUE))
      ))
  }

  if (plot_style == "histogram") {
    plot <- input %>%
      ggplot2::ggplot(ggplot2::aes(x = {{ intensity_log2 }})) +
      ggplot2::geom_histogram(
        binwidth = 0.5,
        color = "black",
        fill = "#5680C1"
      ) +
      ggplot2::labs(
        title = "Overall log2 Intensity Distribution",
        x = "Log2 Intensity",
        y = "Frequency"
      ) +
      {
        if (!missing(sample)) ggplot2::facet_wrap(rlang::new_formula(NULL, rlang::enquo(sample)))
      } +
      ggplot2::theme_bw() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(size = 20),
        axis.title.x = ggplot2::element_text(size = 15),
        axis.text.y = ggplot2::element_text(size = 15),
        axis.text.x = ggplot2::element_text(size = 12),
        axis.title.y = ggplot2::element_text(size = 15),
        strip.text = ggplot2::element_text(size = 15),
        strip.background = element_blank()
      )
    return(plot)
  }

  if (plot_style == "boxplot") {
    if (missing(sample)) {
      stop(strwrap("Please provide a column with sample
information when choosing boxplot as plot style!",
        prefix = "\n", initial = ""
      ))
    }
    plot <- input %>%
      ggplot2::ggplot(aes(x = {{ sample }}, y = {{ intensity_log2 }})) +
      geom_boxplot(fill = "#5680C1", outlier.color = "#B96DAD") +
      labs(title = "Run intensities", x = "", y = "Intensity") +
      theme_bw() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(size = 20),
        axis.title.x = ggplot2::element_text(size = 15),
        axis.text.y = ggplot2::element_text(size = 15),
        axis.text.x = ggplot2::element_text(size = 12, angle = 75, hjust = 1),
        axis.title.y = ggplot2::element_text(size = 15)
      )
    return(plot)
  }

  if (plot_style == "violin") {
    if (missing(sample)) {
      stop(strwrap("Please provide a column with sample
information when choosing violin as plot style!",
        prefix = "\n", initial = ""
      ))
    }
    plot <- input %>%
      ggplot2::ggplot(aes(x = {{ sample }}, y = {{ intensity_log2 }})) +
      ggplot2::geom_violin(fill = "#5680C1", na.rm = TRUE) +
      geom_boxplot(width = 0.15, fill = "white", na.rm = TRUE, alpha = 0.6) +
      labs(title = "Run intensities", x = "", y = "Intensity") +
      theme_bw() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(size = 20),
        axis.title.x = ggplot2::element_text(size = 15),
        axis.text.y = ggplot2::element_text(size = 15),
        axis.text.x = ggplot2::element_text(size = 12, angle = 75, hjust = 1),
        axis.title.y = ggplot2::element_text(size = 15)
      )
    return(plot)
  }
}
