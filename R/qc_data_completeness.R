#' Data completeness
#'
#' Calculates the percentage of data completeness. That means, what percentage of all detected precursors is present in each sample. 
#' 
#' @param data A data frame containing at least the input variables.
#' @param sample The name of the column containing the sample names.
#' @param grouping The name of the column containing either precursor or peptide identifiers. 
#' @param intensity The name of the column containing raw intensity values.
#' @param digestion Optional column indicating the mode of digestion (limited proteolysis or tryptic digest). Alternatively, any other variable
#' by which the data should be split can be provided.
#' @param plot Logical, if TRUE a plot is returned. If FALSE a table is returned.
#' @param interactive Logical, if TRUE the plot is interactive using plotly.
#'
#' @return A bar plot that displays the percentage of data completeness over all samples. If \code{plot = FALSE} a data frame is returned. 
#' If \code{interactive = TRUE}, the plot is interactive.
#' @import dplyr
#' @import ggplot2
#' @importFrom purrr map2_df
#' @importFrom tidyr complete
#' @importFrom plotly ggplotly
#' @importFrom magrittr %>%
#' @importFrom rlang .data := new_formula enquo
#' @importFrom stringr str_sort
#' @export
#'
#' @examples
#' \dontrun{
#' qc_data_completeness(
#' data, 
#' sample = r_file_name, 
#' grouping = eg_precursor_id,
#' intensity = fg_quantity,
#' digestion = digestion
#' )
#' }
qc_data_completeness <- function(data, sample, grouping, intensity, digestion = NULL, plot = TRUE, interactive = FALSE){
  .=NULL
  
  if(!missing(digestion)){
    result <- data %>% 
      split(dplyr::pull(data, {{digestion}})) %>% 
      purrr::map2_df(.y = names(.),
                     .f = ~ .x %>% 
                       dplyr::distinct({{sample}}, {{grouping}}, {{intensity}}) %>% 
                       tidyr::complete({{sample}}, {{grouping}}) %>% 
                       dplyr::group_by({{sample}}) %>% 
                       dplyr::summarise(completeness = sum(!is.na({{intensity}})) / dplyr::n() * 100, .groups = "drop") %>% 
                       dplyr::mutate({{digestion}} := .y)
      )
  } else {
    result <- data %>% 
      dplyr::distinct({{sample}}, {{grouping}}, {{intensity}}) %>% 
      tidyr::complete({{sample}}, {{grouping}}) %>% 
      dplyr::group_by({{sample}}) %>% 
      dplyr::summarise(completeness = sum(!is.na({{intensity}})) / dplyr::n() * 100, .groups = "drop")
  }
  
  if(plot == FALSE) return(result)
  
  completeness_plot <- result %>% 
    dplyr::mutate({{sample}} := factor({{sample}}, levels = unique(stringr::str_sort({{sample}}, numeric = TRUE)))) %>% 
    ggplot2::ggplot(ggplot2::aes({{sample}}, .data$completeness)) +
    ggplot2::geom_col(fill = "#5680C1", col = "black", size = 1.5) +
    {if(interactive == FALSE) geom_text(
      data = result,
      aes(label = round(.data$completeness, digits = 1)),
      position = position_stack(vjust = 0.5)
    )} +
    ggplot2::scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) +
    {if(!missing(digestion)) ggplot2::facet_wrap(rlang::new_formula(NULL, rlang::enquo(digestion)), scales = "free", ncol = 2)} +
    ggplot2::labs(title = "Data completeness per .raw file", x = "", y = "Data Completeness [%]") +
    ggplot2::theme_bw() +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 20),
                   axis.title.x = ggplot2::element_text(size = 15),
                   axis.text.y = ggplot2::element_text(size = 15),
                   axis.text.x = ggplot2::element_text(size = 12, angle = 75, hjust =1),
                   axis.title.y = ggplot2::element_text(size = 15),
                   panel.border = ggplot2::element_rect(fill = NA),
                   strip.text = ggplot2::element_text(size = 15),
                   strip.background = element_rect(fill = "white")
    )
  
  if(interactive == FALSE) return(completeness_plot)
  
  suppressWarnings(plotly::ggplotly(completeness_plot))
}