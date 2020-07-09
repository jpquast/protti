#' Check run intensities of all samples
#'
#' Plots the intensities in each sample as a boxplot
#'
#' @param data A dataframe containing at least sample names, peptide or precursor identifiers and (log2 transformed) intensities for each peptide or precursor.
#' @param sample The column in the data dataframe containing the sample name.
#' @param grouping The column in the data dataframe containing either precursor or peptide identifiers.
#' @param intensity The column in the data dataframe containing the (log2 transformed) intensities of each precursor or peptide.
#'
#' @return A plot with one boxplot for each sample, depicting the peptide or precurosor intensities.
#' @import dplyr
#' @import ggplot2
#' @importFrom magrittr %>%
#' @export
#'
#' @examples
#' \dontrun{
#' qc_run_intensity(
#' data,
#' sample = r_file_name,
#' grouping = pep_stripped_sequence,
#' intensity = normalised_intensity_log2)
#' }
qc_run_intensity <-
  function(data, sample, grouping, intensity)
  {
    result <- data %>%
      dplyr::distinct({{sample}}, {{grouping}}, {{intensity}}) %>%
      dplyr::filter(!is.na({{intensity}}))

    plot <- result %>%
      ggplot2::ggplot(aes(x = {{sample}}, y = {{intensity}})) +
      geom_boxplot(fill = "cornflowerblue", outlier.color = "orchid3") +
      theme_bw() +
      labs(title = "Run intensities", x = "sample", y = "intensity")
    return(plot)
  }

