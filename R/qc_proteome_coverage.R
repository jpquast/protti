#' Proteome coverage per sample and total
#'
#' Calculates the proteome coverage for each samples and for all samples combined. In other words the fraction of detected proteins to all proteins in the proteome is calculated.
#'
#' @param data A data frame containing at least sample names and protein ID's.
#' @param sample The column in the data data frame containing the sample name.
#' @param protein_id The column in the data data frame containing protein identifiers such as UniProt accessions.
#' @param organism_id The NCBI taxonomy identifier (TaxId) of the organism used. Human: 9606, S. cerevisiae: 559292, E. coli: 83333.
#' @param plot A logical indicating whether the result should be plotted (default is TRUE).
#' @param interactive A logical indicating whether the plot should be interactive (default is TRUE).
#'
#' @return A bar plot showing the percentage of of the proteome detected and undetected in total and for each sample. If \code{plot = FALSE} a data frame containing the numbers is returned.
#' @import dplyr
#' @import ggplot2
#' @importFrom magrittr %>%
#' @importFrom forcats fct_rev
#' @importFrom rlang .data := ensym
#' @importFrom tidyr pivot_longer
#' @importFrom plotly ggplotly
#' @export
#'
#' @examples
#' \dontrun{
#' qc_proteome_coverage(
#' data,
#' sample = r_file_name,
#' protein_id = pg_protein_accession,
#' organism_id = 9606
#' )
#' }
qc_proteome_coverage <- function(data, sample, protein_id, organism_id, plot = TRUE, interactive = TRUE) {
  proteins_total <- data %>%
    dplyr::summarize(proteins_detected = dplyr::n_distinct(!!ensym(protein_id)), .groups = "drop") %>%
    dplyr::mutate({{sample}} := "Total")

  proteome <- fetch_uniprot_proteome({{organism_id}}) %>%
    dplyr::summarize(proteins_proteome = dplyr::n_distinct(.data$id), .groups = "drop")
  
  proteome_coverage <- data %>%
    dplyr::group_by({{sample}}) %>%
    dplyr::summarize(proteins_detected = dplyr::n_distinct(!!ensym(protein_id)), .groups = "drop") %>%
    dplyr::bind_rows(proteins_total) %>%
    dplyr::mutate(proteins_undetected = proteome$proteins_proteome - .data$proteins_detected) %>%
    dplyr::mutate(proteins_undetected = .data$proteins_undetected/proteome$proteins_proteome*100, 
                  proteins_detected = .data$proteins_detected/proteome$proteins_proteome*100) %>%
    tidyr::pivot_longer(cols = c(.data$proteins_detected, .data$proteins_undetected), names_to = "type", values_to = "percentage") %>%
    dplyr::mutate({{sample}} := stats::relevel(factor({{sample}}), ref = "Total"),
           type = forcats::fct_rev(factor(.data$type)))

  proteome_coverage_plot <- proteome_coverage %>%
    ggplot2::ggplot(ggplot2::aes({{sample}}, .data$percentage, fill = .data$type)) +
    ggplot2::geom_col(col = "black") +
    ggplot2::labs(title = "Proteome coverage", x = "", y = "Proteome [%]")+
    ggplot2::scale_fill_manual(values = c("proteins_detected" = "cornflowerblue", "proteins_undetected" = "brown1"), name = "Proteins", labels = c("Not detected", "Detected"))+
    ggplot2::geom_text(ggplot2::aes(label = round(.data$percentage, digits = 1)), position = ggplot2::position_stack(vjust = 0.5), size=4)+
    ggplot2::theme_bw()+
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust =1))

  if(plot == FALSE) return(proteome_coverage)
  if(plot == TRUE){
    if(interactive == FALSE) return(proteome_coverage_plot)
    if(interactive == TRUE) {
      proteome_coverage_plot <- proteome_coverage %>%
        ggplot2::ggplot(ggplot2::aes({{sample}}, .data$percentage, fill = .data$type)) +
        ggplot2::geom_col(col = "black") +
        ggplot2::labs(title = "Proteome coverage per .raw file", x = "", y = "Proteome [%]")+
        ggplot2::scale_fill_manual(values = c("proteins_detected" = "cornflowerblue", "proteins_undetected" = "brown1"), name = "Proteins")+
        ggplot2::theme_bw()+
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust =1))
      plotly::ggplotly(proteome_coverage_plot)
    }
  }
}