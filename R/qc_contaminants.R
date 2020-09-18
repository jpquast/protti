#' Percentage of contaminants per sample
#'
#' Calculates the percentage of contaminating proteins as the share of total intensity
#' 
#' @param data A data frame containing at least the input variables.
#' @param sample The name of the column containing the sample names.
#' @param protein The name of the column containing protein IDs or protein names.
#' @param is_contaminant The name of the column containing a logical indicating if the protein is a contaminant.
#' @param intensity The name of the column containing raw intensity values.
#' @param n_contaminants Numeric, indicating how many contaminants should be displayed individually. The rest is combined to a group called "other". The 
#' default is 5.
#' @param plot Logical, if TRUE a plot is returned. If FALSE a table is returned.
#' @param interactive logical, if TRUE the plot is interactive using plotly.
#'
#' @return A bar plot that displays the percentage of contaminating proteins over all samples. If \code{plot = FALSE} a data frame is returned.
#' @import dplyr
#' @import ggplot2
#' @importFrom utils head
#' @importFrom forcats fct_inorder
#' @importFrom plotly ggplotly
#' @importFrom magrittr %>%
#' @importFrom rlang .data :=
#' @export
#'
#' @examples
#' \dontrun{
#' qc_contaminants(
#' data, 
#' sample = sample, 
#' protein = leading_razor_protein,
#' is_contaminant = potential_contaminant
#' intensity = intensity
#' )
#' }
qc_contaminants <- function(data, sample, protein, is_contaminant, intensity, n_contaminants = 5, plot = TRUE, interactive = FALSE){
  result <- data %>% 
    dplyr::distinct({{sample}}, {{protein}}, {{is_contaminant}}, {{intensity}}) %>%
    dplyr::group_by({{sample}}) %>%
    dplyr::mutate(total_intensity = sum({{intensity}}, na.rm = TRUE)) %>%
    dplyr::filter({{is_contaminant}} == TRUE) %>%
    dplyr::group_by({{sample}}, {{protein}}, .data$total_intensity) %>%
    dplyr::summarise(sum_protein_intensity = sum({{intensity}}, na.rm = TRUE), .groups = "drop") %>%
    dplyr::mutate(contaminant_percentage = .data$sum_protein_intensity / .data$total_intensity * 100) %>%
    dplyr::group_by({{sample}}) %>%
    dplyr::mutate({{protein}} := ifelse(.data$contaminant_percentage %in% utils::head(sort(.data$contaminant_percentage, decreasing = TRUE), n_contaminants),
      {{protein}},
      "other"
    )) %>%
    dplyr::arrange(.data$contaminant_percentage) %>%
    dplyr::ungroup() %>%
    dplyr::mutate({{protein}} := forcats::fct_inorder(factor({{protein}}))) %>%
    dplyr::group_by({{sample}}, {{protein}}) %>%
    dplyr::summarise(contaminant_percentage = sum(.data$contaminant_percentage), .groups = "drop")

  if(plot == FALSE) return(result)

  plot_result <- result %>%
    ggplot2::ggplot(ggplot2::aes({{sample}}, .data$contaminant_percentage, fill = {{protein}})) +
    ggplot2::geom_col(col = "black", size = 1.2) +
    # ggplot2::scale_fill_manual(values = c("#C1AE7C", "#F1DABF", "#F3CA40", "#F2A541","#F0B67F", "#F5E6E8", "#B9FAF8", "#AAA1C8", "#2F4858", "#33658A")) +
    # https://coolors.co/33658a-2f4858-aaa1c8-b9faf8-f5e6e8-f0b67f-f2a541-f3ca40-f1dabf-c1ae7c There should be more colors added to also support n_contaminants higher than 5
    ggplot2::labs(title = "Contaminants per .raw file", x = "Sample", y = "% of total intensity", fill = "Contaminant Protein") +
    ggplot2::theme_bw() +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 20),
                   axis.title.x = ggplot2::element_text(size = 15),
                   axis.text.y = ggplot2::element_text(size = 15),
                   axis.text.x = ggplot2::element_text(size = 12, angle = 45, hjust =1),
                   axis.title.y = ggplot2::element_text(size = 15),
                   legend.title = ggplot2::element_text(size = 15),
                   legend.text = ggplot2::element_text(size = 15),
                   panel.border = ggplot2::element_rect(fill = NA)
    )

  if(interactive == FALSE) return(plot_result)
  
  suppressWarnings(plotly::ggplotly(plot_result))
  
}