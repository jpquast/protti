#' Check number of precursor/preptide/protein IDs
#'
#' Plots the number of IDs for each sample
#'
#' @param data Dataframe containing at least sample names and precursor IDs.
#' @param sample Column in the dataframe specifying the sample name.
#' @param grouping Column in the data dataframe containing either precursor or peptide identifiers.
#' @param condition Optional column in the dataframe specifying the condition of the sample (e.g. LiP_reated, LiP_untreated).
#' @param title Optional argument specifying the plot title (default is "ID count per sample").
#'
#'
#' @return A bar plot with the height corresponding to the number of IDs, each bar represents one sample.
#' @import dplyr
#' @import ggplot2
#' @importFrom magrittr %>%
#' @importFrom plotly ggplotly
#' @export
#'
#' @examples
#' \dontrun{
#' qc_ids(
#' data,
#' sample = r_file_name
#' grouping = eg_precursor_id,
#' condition = r_condition,
#' title = "Number of peptide IDs per sample"
#' )
#' }
qc_ids <-
  function(data, sample, grouping, condition = NULL, title = "ID count per sample")
  {
    result <- data %>%
      dplyr::select({{sample}}, {{grouping}}, {{condition}}) %>%
      dplyr::distinct() %>%
      dplyr::group_by({{sample}}) %>%
      dplyr::mutate(count = length({{grouping}})) %>%
      dplyr::distinct()

    plot <- result %>%
      ggplot2::ggplot(aes(x = {{sample}}, y = .data$count, col = {{condition}})) +
      geom_bar(stat = "identity", fill = "cornflowerblue", alpha = 1) +
      labs(title = title,
           x = "sample",
           y = "count") +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 75, hjust = 1),
            legend.position = "none")
    return(plotly::ggplotly(plot))
  }
