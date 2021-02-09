#' Check peptide type percentage share
#'
#' Calculates the percentage share of each peptide types (fully-tryptic, semi-tryptic, non-tryptic) for each sample.
#'
#' @param data A dataframe containing at least the input columns.
#' @param sample the name of the column containing the sample names.
#' @param peptide the name of the column containing the peptide sequence.
#' @param pep_type the name of the column containing the peptide type. Can be obtained using the \code{find_peptide} and \code{peptide_type} function together.
#' @param intensity column containing peptide intensity values (not log2 transformed).
#' @param remove_na_intensities Logical specifying if sample/peptide combinations with intensities that are NA (not quantified IDs) should 
#' be dropped from the data frame for analysis of peptide type distributions. Default is TRUE since we are usually 
#' interested in the peptide type distribution of quantifiable IDs.
#' @param method method used for calculation. \code{method = "intensity"} calculates the peptide type percentage by intensity, whereas \code{method = "count"} calculates the percentage by peptide ID count. Default is \code{method = count}.
#' @param plot a logical indicating whether the result should be plotted.
#' @param interactive a logical indicating whether the plot should be interactive.
#'
#' @return A data frame that contains the calculated percentage shares of each peptide type per sample. The \code{count} column contains the number of peptides with a specific type. The \code{peptide_type_percent} column contains the percentage share of a specific peptide type.
#' @import dplyr
#' @import tidyr
#' @import ggplot2
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @importFrom plotly ggplotly
#' @importFrom tidyr drop_na
#' @importFrom utils data
#' @importFrom stringr str_sort
#' @export
#'
#' @examples
#' \dontrun{
#' qc_peptide_type(
#' data,
#' sample = r_file_name,
#' peptide = pep_stripped_sequence,
#' pep_type = pep_type,
#' intensity = fg_quantity,
#' method = "intensity",
#' plot = TRUE)
#' }
qc_peptide_type <- function(data, sample, peptide, pep_type, intensity, remove_na_intensities = TRUE, method = "count", plot = FALSE, interactive = FALSE) {
  protti_colours <- "placeholder" # assign a placeholder to prevent a missing global variable warning
  utils::data("protti_colours", envir=environment()) # then overwrite it with real data
  if(remove_na_intensities == TRUE){
    data <- data %>% 
      tidyr::drop_na({{intensity}})
  }
  if(method == "count")
  {
    result <- data %>%
      dplyr::distinct({{sample}}, {{peptide}}, {{pep_type}}) %>%
      tidyr::drop_na({{pep_type}}) %>%
      dplyr::count({{sample}}, {{pep_type}}, name = "count") %>%
      dplyr::group_by({{sample}}) %>%
      dplyr::mutate(peptide_type_percent = .data$count/sum(.data$count)*100) %>%
      dplyr::ungroup() %>% 
      dplyr::distinct({{sample}}, {{pep_type}}, .data$peptide_type_percent) %>%
      dplyr::mutate(pep_type = factor({{pep_type}}, levels = c("fully-tryptic", "semi-tryptic", "non-tryptic")))

    if (plot == TRUE & interactive == FALSE){
      plot <- result %>%
        dplyr::mutate({{sample}} := factor({{sample}}, levels = unique(stringr::str_sort({{sample}}, numeric = TRUE)))) %>% 
        ggplot2::ggplot(ggplot2::aes({{sample}}, .data$peptide_type_percent, fill = .data$pep_type)) +
        ggplot2::geom_col(col = "black", size = 1) +
        geom_text(
          data = result %>% dplyr::filter(.data$peptide_type_percent > 10),
          aes(label = round(.data$peptide_type_percent, digits = 1)),
          position = position_stack(vjust = 0.9)
        ) +
        ggplot2::labs(title = "Peptide types per .raw file", x = "", y = "Percentage of peptides", fill = "Type") +
        ggplot2::theme_bw() +
        ggplot2::theme(plot.title = ggplot2::element_text(size = 20),
                       axis.title.x = ggplot2::element_text(size = 15),
                       axis.text.y = ggplot2::element_text(size = 15),
                       axis.text.x = ggplot2::element_text(size = 12, angle = 75, hjust =1),
                       axis.title.y = ggplot2::element_text(size = 15),
                       legend.title = ggplot2::element_text(size = 15),
                       legend.text = ggplot2::element_text(size = 15)) +
        ggplot2::scale_fill_manual(values = protti_colours)

      return(plot)
    }
    if (plot == TRUE & interactive == TRUE) {
      plot <- result %>%
        dplyr::mutate({{sample}} := factor({{sample}}, levels = unique(stringr::str_sort({{sample}}, numeric = TRUE)))) %>% 
        ggplot2::ggplot(ggplot2::aes({{sample}}, .data$peptide_type_percent, fill = .data$pep_type)) +
        ggplot2::geom_col(col = "black", size = 1) +
        ggplot2::labs(title = "Peptide types per .raw file", x = "Sample", y = "Percentage of peptides", fill = "Type") +
        ggplot2::theme_bw() +
        ggplot2::theme(plot.title = ggplot2::element_text(size = 20),
                       axis.title.x = ggplot2::element_text(size = 15),
                       axis.text.y = ggplot2::element_text(size = 15),
                       axis.text.x = ggplot2::element_text(size = 12, angle = 75, hjust =1),
                       axis.title.y = ggplot2::element_text(size = 15),
                       legend.title = ggplot2::element_text(size = 15),
                       legend.text = ggplot2::element_text(size = 15)) +
        ggplot2::scale_fill_manual(values = protti_colours)

      interactive_plot <- plotly::ggplotly(plot)

      return(interactive_plot)
    }
    if (plot == FALSE){
      return(result)
    }
  }

  if(method == "intensity")
  {
    result <- data %>%
      tidyr::drop_na({{intensity}}) %>% 
      dplyr::distinct({{sample}}, {{peptide}}, {{pep_type}}, {{intensity}}) %>%
      tidyr::drop_na({{pep_type}}) %>%
      dplyr::group_by({{sample}}) %>%
      dplyr::mutate(total_int = sum({{intensity}})) %>%
      dplyr::group_by({{sample}}, {{pep_type}}) %>%
      dplyr::mutate(pep_type_int = sum({{intensity}})) %>%
      dplyr::group_by({{sample}}) %>%
      dplyr::mutate(peptide_type_percent = (.data$pep_type_int / .data$total_int) * 100) %>%
      dplyr::ungroup() %>% 
      dplyr::distinct({{sample}}, {{pep_type}}, .data$peptide_type_percent) %>%
      dplyr::mutate(pep_type = factor({{pep_type}}, levels = c("fully-tryptic", "semi-tryptic", "non-tryptic")))

    if (plot == TRUE & interactive == FALSE){
      plot <- result %>%
        dplyr::mutate({{sample}} := factor({{sample}}, levels = unique(stringr::str_sort({{sample}}, numeric = TRUE)))) %>% 
        ggplot2::ggplot(ggplot2::aes({{sample}}, .data$peptide_type_percent, fill = .data$pep_type)) +
        ggplot2::geom_col(col = "black", size = 1) +
        geom_text(
          data = result %>% dplyr::filter(.data$peptide_type_percent > 10),
          aes(label = round(.data$peptide_type_percent, digits = 1)),
          position = position_stack(vjust = 0.9)
        ) +
        ggplot2::labs(title = "Peptide type intensity per .raw file", x = "Sample", y = "Percentage of total peptide intensity", fill = "Type") +
        ggplot2::theme_bw() +
        ggplot2::theme(plot.title = ggplot2::element_text(size = 20),
                       axis.title.x = ggplot2::element_text(size = 15),
                       axis.text.y = ggplot2::element_text(size = 15),
                       axis.text.x = ggplot2::element_text(size = 12, angle = 75, hjust =1),
                       axis.title.y = ggplot2::element_text(size = 15),
                       legend.title = ggplot2::element_text(size = 15),
                       legend.text = ggplot2::element_text(size = 15)) +
        ggplot2::scale_fill_manual(values = protti_colours)
      return(plot)
    }
    if (plot == TRUE & interactive == TRUE) {
      plot <- result %>%
        dplyr::mutate({{sample}} := factor({{sample}}, levels = unique(stringr::str_sort({{sample}}, numeric = TRUE)))) %>% 
        ggplot2::ggplot(ggplot2::aes({{sample}}, .data$peptide_type_percent, fill = .data$pep_type)) +
        ggplot2::geom_col(col = "black", size = 1) +
        ggplot2::labs(title = "Peptide type intensity per .raw file", x = "Sample", y = "Percentage of total peptide intensity", fill = "Type") +
        ggplot2::theme_bw() +
        ggplot2::theme(plot.title = ggplot2::element_text(size = 20),
                       axis.title.x = ggplot2::element_text(size = 15),
                       axis.text.y = ggplot2::element_text(size = 15),
                       axis.text.x = ggplot2::element_text(size = 12, angle = 75, hjust =1),
                       axis.title.y = ggplot2::element_text(size = 15),
                       legend.title = ggplot2::element_text(size = 15),
                       legend.text = ggplot2::element_text(size = 15)) +
        ggplot2::scale_fill_manual(values = protti_colours)

      interactive_plot <- plotly::ggplotly(plot)

      return(interactive_plot)
    }
    if (plot == FALSE){
      return(result)
    }

  }
}
