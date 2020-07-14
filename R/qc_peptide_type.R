#' Check peptide type percentage share
#'
#' Calculates the percentage share of each peptide types (fully-tryptic, semi-tryptic, non-tryptic) for each sample.
#'
#' @param data A dataframe containing at least the input columns.
#' @param sample the name of the column containing the sample names.
#' @param peptide the name of the column containing the peptide sequence.
#' @param pep_type the name of the column containing the peptide type. Can be obtained using the \code{find_peptide} and \code{peptide_type} function together.
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
#' @export
#'
#' @examples
#' \dontrun{
#' qc_peptide_type(
#' data,
#' sample = r_file_name,
#' peptide = pep_stripped_sequence,
#' pep_type = pep_type,
#' plot = TRUE)
#' }
qc_peptide_type <- function(data, sample, peptide, pep_type, plot = FALSE, interactive = FALSE){
  result <- data %>%
    dplyr::distinct({{sample}}, {{peptide}}, {{pep_type}}) %>%
    tidyr::drop_na({{pep_type}}) %>%
    dplyr::count({{sample}}, {{pep_type}}, name = "count") %>%
    dplyr::group_by({{sample}}) %>%
    dplyr::mutate(peptide_type_percent = .data$count/sum(.data$count)*100)
  
  if (plot == TRUE & interactive == FALSE){
    plot <- result %>%
      ggplot2::ggplot(ggplot2::aes({{sample}}, .data$peptide_type_percent, fill = {{pep_type}})) +
      ggplot2::geom_bar(stat = "identity", position = "dodge", col = "black") + 
      ggplot2::labs(title = "Peptide types per .raw file", x = "Sample", y = "Percentage of peptides", fill = "Type") + 
      ggplot2::geom_text(ggplot2::aes(label = round(.data$peptide_type_percent, digits = 1), y = .data$peptide_type_percent + max(.data$peptide_type_percent) * 0.03), position = position_dodge(1), size = 4) +
      ggplot2::theme_bw() + 
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 75, hjust = 1)) + 
      ggplot2::coord_cartesian(ylim = c(0, (max(result$peptide_type_percent))))
    
    return(plot)
  } 
  if (plot == TRUE & interactive == TRUE) {
    plot <- result %>%
      ggplot2::ggplot(ggplot2::aes({{sample}}, .data$peptide_type_percent, fill = {{pep_type}})) +
      ggplot2::geom_bar(stat = "identity", position = "dodge", col = "black") + 
      ggplot2::labs(title = "Peptide types per .raw file", x = "Sample", y = "Percentage of peptides", fill = "Type") + 
      ggplot2::theme_bw() + 
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 75, hjust = 1)) + 
      ggplot2::coord_cartesian(ylim = c(0, (max(result$peptide_type_percent))))
    
    interactive_plot <- plotly::ggplotly(plot)
    
    return(interactive_plot)
  }
  if (plot == FALSE){
    return(result)
  }
}