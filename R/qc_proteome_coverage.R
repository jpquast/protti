#' Proteome coverage per sample and total
#'
#' Calculates the proteome coverage for each samples and for all samples combined. In other words t
#' he fraction of detected proteins to all proteins in the proteome is calculated.
#'
#' @param data a data frame that contains at least sample names and protein ID's.
#' @param sample a character column in the \code{data} data frame that contains the sample name.
#' @param protein_id a character or numeric column in the \code{data} data frame that contains
#' protein identifiers such as UniProt accessions.
#' @param organism_id a numeric value that specifies a NCBI taxonomy identifier (TaxId) of the
#' organism used. Human: 9606, S. cerevisiae: 559292, E. coli: 83333.
#' @param reviewed a logical value that determines if only reviewed protein entries will be considered
#' as the full proteome. Default is TRUE.
#' @param plot a logical value that specifies whether the result should be plotted.
#' @param interactive a logical value that indicates whether the plot should be interactive
#' (default is FALSE).
#'
#' @return A bar plot showing the percentage of of the proteome detected and undetected in total
#' and for each sample. If \code{plot = FALSE} a data frame containing the numbers is returned.
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
#' \donttest{
#' # Create example data
#' proteome <- data.frame(id = 1:4518)
#' data <- data.frame(
#'   sample = c(rep("A", 101), rep("B", 1000), rep("C", 1000)),
#'   protein_id = c(proteome$id[1:100], proteome$id[1:1000], proteome$id[1000:2000])
#' )
#'
#' # Calculate proteome coverage
#' qc_proteome_coverage(
#'   data = data,
#'   sample = sample,
#'   protein_id = protein_id,
#'   organism_id = 83333,
#'   plot = FALSE
#' )
#'
#' # Plot proteome coverage
#' qc_proteome_coverage(
#'   data = data,
#'   sample = sample,
#'   protein_id = protein_id,
#'   organism_id = 83333,
#'   plot = TRUE
#' )
#' }
qc_proteome_coverage <- function(data,
                                 sample,
                                 protein_id,
                                 organism_id,
                                 reviewed = TRUE,
                                 plot = TRUE,
                                 interactive = FALSE) {
  proteins_total <- data %>%
    dplyr::summarize(proteins_detected = dplyr::n_distinct(!!ensym(protein_id)), .groups = "drop") %>%
    dplyr::mutate({{ sample }} := "Total")

  proteome <- fetch_uniprot_proteome(organism_id, reviewed = reviewed)

  if (is(proteome, "character")) {
    # UniProt information could not be fetched.
    message("UniProt information could not be fetched")
    return(NULL)
  }

  proteome <- proteome %>%
    dplyr::summarize(proteins_proteome = dplyr::n_distinct(.data$accession), .groups = "drop")

  proteome_coverage <- data %>%
    dplyr::group_by({{ sample }}) %>%
    dplyr::summarize(proteins_detected = dplyr::n_distinct(!!ensym(protein_id)), .groups = "drop") %>%
    dplyr::bind_rows(proteins_total) %>%
    dplyr::mutate(proteins_undetected = proteome$proteins_proteome - .data$proteins_detected) %>%
    dplyr::mutate(
      proteins_undetected = .data$proteins_undetected / proteome$proteins_proteome * 100,
      proteins_detected = .data$proteins_detected / proteome$proteins_proteome * 100
    ) %>%
    tidyr::pivot_longer(
      cols = c("proteins_detected", "proteins_undetected"),
      names_to = "type",
      values_to = "percentage"
    ) %>%
    dplyr::mutate({{ sample }} := stats::relevel(factor({{ sample }}), ref = "Total"),
      type = forcats::fct_rev(factor(.data$type))
    )

  proteome_coverage_plot <- proteome_coverage %>%
    ggplot2::ggplot(ggplot2::aes({{ sample }}, .data$percentage, fill = .data$type)) +
    ggplot2::geom_col(col = "black", size = 1) +
    ggplot2::labs(
      title = "Proteome coverage per .raw file",
      x = "",
      y = "Proteome [%]"
    ) +
    ggplot2::scale_fill_manual(
      values = c("proteins_detected" = "#5680C1", "proteins_undetected" = "#B96DAD"),
      name = "Proteins",
      labels = c("Not detected", "Detected")
    ) +
    ggplot2::geom_text(
      data = proteome_coverage %>% dplyr::filter(.data$percentage > 5),
      ggplot2::aes(label = round(.data$percentage, digits = 1)),
      position = ggplot2::position_stack(vjust = 0.5),
      size = 4
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 20),
      axis.title.x = ggplot2::element_text(size = 15),
      axis.text.y = ggplot2::element_text(size = 15),
      axis.text.x = ggplot2::element_text(size = 12, angle = 75, hjust = 1),
      axis.title.y = ggplot2::element_text(size = 15),
      legend.title = ggplot2::element_text(size = 15),
      legend.text = ggplot2::element_text(size = 15)
    )

  if (plot == FALSE) {
    return(proteome_coverage)
  }
  if (plot == TRUE) {
    if (interactive == FALSE) {
      return(proteome_coverage_plot)
    }
    if (interactive == TRUE) {
      proteome_coverage_plot <- proteome_coverage %>%
        ggplot2::ggplot(ggplot2::aes({{ sample }}, .data$percentage, fill = .data$type)) +
        ggplot2::geom_col(col = "black", size = 1) +
        ggplot2::labs(
          title = "Proteome coverage per .raw file",
          x = "",
          y = "Proteome [%]"
        ) +
        ggplot2::scale_fill_manual(
          values = c("proteins_detected" = "#5680C1", "proteins_undetected" = "#B96DAD"),
          name = "Proteins"
        ) +
        ggplot2::theme_bw() +
        ggplot2::theme(
          plot.title = ggplot2::element_text(size = 20),
          axis.title.x = ggplot2::element_text(size = 15),
          axis.text.y = ggplot2::element_text(size = 15),
          axis.text.x = ggplot2::element_text(size = 12, angle = 75, hjust = 1),
          axis.title.y = ggplot2::element_text(size = 15),
          legend.title = ggplot2::element_text(size = 15),
          legend.text = ggplot2::element_text(size = 15)
        )
      plotly::ggplotly(proteome_coverage_plot)
    }
  }
}
