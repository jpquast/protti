#' Check number of precursor, peptide or protein IDs
#'
#' Returns a plot or table of the number of IDs for each sample.
#'
#' @param data Dataframe containing at least sample names and precursor/peptide/protein IDs.
#' @param sample Column in the data frame specifying the sample name.
#' @param grouping Column in the data frame containing either precursor, peptide or protein identifiers.
#' @param intensity Column in the data frame containing intensities. If \code{remove_na_intensities = FALSE}, this argument is not required.
#' @param remove_na_intensities Logical specifying if sample/grouping combinations with intensities that are NA (not quantified IDs) should be droped from the data frame.
#' Default is TRUE since we are usually interested in the number of quantifiable IDs.
#' @param condition Optional column in the data frame specifying the condition of the sample (e.g. LiP_treated, LiP_untreated),
#' if column is provided, the bars in the plot will be coloured according to the condition.
#' @param title Optional argument specifying the plot title (default is "ID count per sample").
#' @param plot Argument specifying whether the output of the function should be plotted (default is TRUE).
#' @param interactive Argument specifying whether the plot should be interactive (default is TRUE).
#'
#' @return A bar plot with the height corresponding to the number of IDs, each bar represents one sample
#' (if \code{plot = TRUE}). If \code{plot = FALSE} a table with ID counts is returned.
#' @import dplyr
#' @import ggplot2
#' @importFrom magrittr %>%
#' @importFrom tidyr drop_na
#' @importFrom plotly ggplotly
#' @importFrom utils data
#' @export
#'
#' @examples
#' \dontrun{
#' qc_ids(
#' data,
#' sample = r_file_name,
#' grouping = eg_precursor_id,
#' condition = r_condition,
#' title = "Number of peptide IDs per sample"
#' )
#' }
qc_ids <-
  function(data, sample, grouping, intensity = NULL, remove_na_intensities = TRUE, condition = NULL, title = "ID count per sample", plot = TRUE, interactive = TRUE) {
    protti_colors <- "placeholder" # assign a placeholder to prevent a missing global variable warning
    utils::data("protti_colors", envir=environment()) # then overwrite it with real data
    if(remove_na_intensities == TRUE){
      data <- data %>%
        dplyr::distinct({{sample}}, {{grouping}}, {{condition}}, {{intensity}}) %>%
        tidyr::drop_na({{intensity}})
    }
   result <- data %>%
      dplyr::distinct({{sample}}, {{grouping}}, {{condition}}) %>%
      dplyr::group_by({{sample}}) %>%
      dplyr::mutate(count = n()) %>%
      dplyr::select(-c({{grouping}})) %>%
      dplyr::distinct()

 if (plot == TRUE)
    {
    plot <- result %>%
      ggplot2::ggplot(aes(x = {{sample}}, y = .data$count, fill = {{condition}})) +
      geom_col(col = "black") +
      {if(length(result %>% ungroup() %>% select({{condition}})) == 0 ) geom_col(fill = "#5680C1", col = "black")}  +
      labs(title = title,
           x = "sample",
           y = "count") +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 75, hjust = 1)) +
      scale_fill_manual(values = protti_colors)
    if (interactive == TRUE)
    {
      return(plotly::ggplotly(plot))
    } else {
      return(plot)
    }
    } else {
      return(result)
    }
  }
