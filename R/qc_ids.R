#' Check number of precursor, peptide or protein IDs
#'
#' Returns a plot or table of the number of IDs for each sample. The default settings remove grouping variables without quantitative information (intensity is NA).
#' These will not be counted as IDs.
#'
#' @param data Data frame containing at least sample names and precursor/peptide/protein IDs.
#' @param sample the column in the data frame specifying the sample name.
#' @param grouping the column in the data frame containing either precursor, peptide or protein identifiers.
#' @param intensity the column in the data frame containing raw or log2 intensities. If \code{remove_na_intensities = FALSE}, this argument is not required.
#' @param remove_na_intensities Logical specifying if sample/grouping combinations with intensities that are NA (not quantified IDs) should be dropped from the data frame.
#' Default is TRUE since we are usually interested in the number of quantifiable IDs.
#' @param condition Optional column in the data frame specifying the condition of the sample (e.g. LiP_treated, LiP_untreated),
#' if column is provided, the bars in the plot will be coloured according to the condition.
#' @param title Optional argument specifying the plot title (default is "ID count per sample").
#' @param plot Argument specifying whether the output of the function should be plotted (default is TRUE).
#' @param interactive Argument specifying whether the plot should be interactive (default is FALSE).
#'
#' @return A bar plot with the height corresponding to the number of IDs, each bar represents one sample
#' (if \code{plot = TRUE}). If \code{plot = FALSE} a table with ID counts is returned.
#' @import dplyr
#' @import ggplot2
#' @importFrom magrittr %>%
#' @importFrom tidyr drop_na
#' @importFrom plotly ggplotly
#' @importFrom utils data
#' @importFrom stringr str_sort
#' @export
#'
#' @examples
#' \dontrun{
#' qc_ids(
#' data,
#' sample = r_file_name,
#' grouping = eg_precursor_id,
#' intensity = fg_quantity,
#' condition = r_condition,
#' title = "Number of peptide IDs per sample"
#' )
#' }
qc_ids <-
  function(data, sample, grouping, intensity, remove_na_intensities = TRUE, condition = NULL, title = "ID count per sample", plot = TRUE, interactive = FALSE) {
    protti_colours <- "placeholder" # assign a placeholder to prevent a missing global variable warning
    utils::data("protti_colours", envir=environment()) # then overwrite it with real data
    if(remove_na_intensities == TRUE){

      if(missing(intensity)) stop("Please provide a column containing intensities or set remove_na_intensities to FALSE")

      data <- data %>%
        tidyr::drop_na({{intensity}}) %>% 
        dplyr::mutate({{condition}} := as.character({{condition}}))
    }

   result <- data %>%
      dplyr::distinct({{sample}}, {{grouping}}, {{condition}}) %>%
      dplyr::group_by({{sample}}) %>%
      dplyr::mutate(count = dplyr::n()) %>%
      dplyr::select(-c({{grouping}})) %>%
      dplyr::distinct() %>%
      dplyr::ungroup()

 if (plot == TRUE)
    {
    plot <- result %>%
      dplyr::mutate({{sample}} := factor({{sample}}, levels = unique(stringr::str_sort({{sample}}, numeric = TRUE)))) %>%
      ggplot2::ggplot(aes(x = {{sample}}, y = .data$count, fill = {{condition}})) +
      ggplot2::geom_col(col = "black", size = 1) +
      {if(missing(condition)) ggplot2::geom_col(fill = "#5680C1", col = "black")}  +
      ggplot2::labs(title = title,
           x = "",
           y = "Count",
           fill = "Condition") +
      ggplot2::theme_bw() +
      ggplot2::scale_fill_manual(values = protti_colours) +
      ggplot2::theme(plot.title = ggplot2::element_text(size = 20),
                     axis.title.x = ggplot2::element_text(size = 15),
                     axis.text.y = ggplot2::element_text(size = 15),
                     axis.text.x = ggplot2::element_text(size = 12, angle = 75, hjust =1),
                     axis.title.y = ggplot2::element_text(size = 15),
                     legend.title = ggplot2::element_text(size = 15),
                     legend.text = ggplot2::element_text(size = 15))
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
