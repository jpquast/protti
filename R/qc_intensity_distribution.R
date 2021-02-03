#' Check intensity distribution per sample and overall
#'
#' Plots the overall or sample-wise distribution of all peptide intensities as a boxplot or histogram.
#'
#' @param data A data frame containing at least sample names, grouping identifiers (precursor, peptide or protein) and
#' log2 transformed intensities for each grouping identifier.
#' @param sample The column in the data frame containing the sample name. NOTE: If the overall distribution should be returned please
#' do not provide the name of the sample column.
#' @param grouping The column in the data frame containing either precursor, peptide or protein identifiers.
#' @param intensity The column in the data frame containing the log2 transformed intensities of each grouping identifier sample combination.
#' @param method A character vector indicating the plot type. This can be either "histogram" or "boxplot". Method "boxplot" can
#' only be used if a sample column is provided.
#'
#' @return A histogram or boxplot that shows the intensity distribution over all samples or by sample.
#' @import ggplot2
#' @importFrom dplyr distinct mutate
#' @importFrom magrittr %>%
#' @importFrom tidyr drop_na
#' @importFrom stringr str_sort
#' @importFrom rlang new_formula enquo
#' @export
#'
#' @examples
#' \dontrun{
#' qc_intensity_distribution(
#'   data,
#'   sample = r_file_name,
#'   grouping = eg_precursor_id,
#'   intensity = normalised_intensity_log2,
#'   method = "boxplot"
#' )
#' }
qc_intensity_distribution <- function(data, sample = NULL, grouping, intensity, method) {
  if (missing(method)) stop("Please provide a plot type in the method argument!")
  input <- data %>%
    dplyr::distinct({{ sample }}, {{ grouping }}, {{ intensity }}) %>%
    tidyr::drop_na({{ intensity }})

  if (!missing(sample)) {
    input <- input %>%
      dplyr::mutate({{ sample }} := factor({{ sample }}, levels = unique(stringr::str_sort({{ sample }}, numeric = TRUE))))
  }

  if (method == "histogram") {
    plot <- input %>%
      ggplot2::ggplot(ggplot2::aes(x = {{ intensity }})) +
      ggplot2::geom_histogram(binwidth = 0.5, color = "black", fill = "#5680C1") +
      ggplot2::labs(title = "Overall log2 Intensity Distribution", x = "Log2 Intensity", y = "Frequency") +
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

  if (method == "boxplot") {
    if (missing(sample)) stop("Please provide a column with sample information when choosing boxplot as method!")
    plot <- input %>%
      ggplot2::ggplot(aes(x = {{ sample }}, y = {{ intensity }})) +
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
}
