#' Check Log2 intensity distribution
#'
#' Plots the overall distribution of all peptide intensities
#'
#' @param data A dataframe containing at least sample names, peptide or precursor identifiers and log2 transformed intensities for each peptide or precursor.
#' @param sample The column in the data dataframe containing the sample name.
#' @param grouping The column in the data dataframe containing either precursor or peptide identifiers.
#' @param log2intensity The column in the data dataframe containing the log2 transformed intensities of each precursor or peptide.
#'
#' @return A plot that shows the intensity distribution over all samples
#' @import dplyr
#' @import ggplot2
#' @importFrom magrittr %>%
#' @export
#'
#' @examples
#' \dontrun{
#' qc_log2_distribution(data,
#' sample,
#' grouping,
#' log2intensity)
#' }
qc_log2_intensity_distribution <-
  function(data, sample, grouping, log2intensity)
  {
    result <- data %>%
    dplyr::distinct({{sample}}, {{grouping}}, {{log2intensity}}) %>%
    dplyr::filter(!is.na({{log2intensity}}))

        plot <- result %>%
          ggplot2::ggplot(aes(x = {{log2intensity}})) +
          geom_histogram(binwidth = 0.5, color = "black", fill = "cornflowerblue") +
          theme_bw() +
          labs(title = "Log2(intensity) distribution over all samples", x = "log2(intensity)", y = "frequency")
        return(plot)
    }
