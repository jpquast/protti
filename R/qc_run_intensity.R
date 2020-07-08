#' Check run intensities of all samples
#'
#' Plots the intensities in each sample as a boxplot
#'
#' @param data A dataframe containing at least sample names, peptide or precursor identifiers and log2 transformed intensities for each peptide or precursor.
#' @param sample The column in the data dataframe containing the sample name.
#' @param grouping The column in the data dataframe containing either precursor or peptide identifiers.
#' @param log2intensity The column in the data dataframe containing the log2 transformed intensities of each precursor or peptide.
#'
#' @return A plot with one boxplot for each sample, depicting the peptide or precurosor intensities.
#' @import dplyr
#' @import ggplot2
#' @importFrom magrittr %>%
#' @export
#'
#' @examples
#' \dontrun{
#' qc_run_intensity(data,
#' sample,
#' grouping,
#' log2intensity)
#' }
qc_run_intensity <-
  function(data, sample, grouping, log2intensity)
  {
    result <- data %>%
      dplyr::distinct({{sample}}, {{grouping}}, {{log2intensity}}) %>%
      dplyr::group_by({{sample}}) %>%
      dplyr::filter(!is.na({{log2intensity}}))

    plot <- result %>%
      ggplot2::ggplot(aes(x = {{sample}}, y = {{log2intensity}})) +
      geom_boxplot(fill = "cornflowerblue", outlier.color = "orchid3") +
      theme_bw() +
      labs(title = "Boxplots of log2 transformed peptide intensities in each sample", x = "sample", y = "log2(intensity)")
    return(plot)
  }

