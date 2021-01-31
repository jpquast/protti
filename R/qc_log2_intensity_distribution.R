#' Check Log2 intensity distribution
#'
#' Plots the overall distribution of all peptide intensities
#'
#' @param data A dataframe containing at least sample names, peptide or precursor identifiers and log2 transformed intensities for each peptide or precursor.
#' @param sample The column in the data dataframe containing the sample name.
#' @param grouping The column in the data dataframe containing either precursor or peptide identifiers.
#' @param log2_intensity The column in the data dataframe containing the log2 transformed intensities of each precursor or peptide.
#'
#' @return A plot that shows the intensity distribution over all samples
#' @import dplyr
#' @import ggplot2
#' @importFrom magrittr %>%
#' @export
#'
#' @examples
#' \dontrun{
#' qc_log2_intensity_distribution(
#' data,
#' sample = r_file_name,
#' grouping = pep_stripped_sequence,
#' log2_intensity = normalised_intensity_log2)
#' }
qc_log2_intensity_distribution <-
  function(data, sample, grouping, log2_intensity)
  {
    result <- data %>%
      dplyr::distinct({{sample}}, {{grouping}}, {{log2_intensity}}) %>%
      dplyr::filter(!is.na({{log2_intensity}}))

    plot <- result %>%
      ggplot2::ggplot(aes(x = {{log2_intensity}})) +
      geom_histogram(binwidth = 0.5, color = "black", fill = "#5680C1") +
      labs(title = "Log2(intensity) distribution over all samples", x = "log2(intensity)", y = "frequency") +
      ggplot2::theme_bw() +
      theme(plot.title = ggplot2::element_text(size = 20),
            axis.title.x = ggplot2::element_text(size = 15),
            axis.text.y = ggplot2::element_text(size = 15),
            axis.text.x = ggplot2::element_text(size = 12),
            axis.title.y = ggplot2::element_text(size = 15))
    return(plot)
  }
