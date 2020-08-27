#' Wood's plot
#'
#' Creates a Wood's plot that plots log2 fold change of peptides or precursors along the protein sequence.
#'
#' @param data Data frame containing differential abundance, start and end peptide or precursor positions, protein 
#' length and optionally a variable based on which peptides or precursors should be colored. 
#' @param fold_change Column in the data data frame containing log2 fold changes.
#' @param start_position Column in the data frame containing the start positions for each peptide or precursor.
#' @param end_position Column in the data frame containing the end positions for each peptide or precursor.
#' @param protein_length Column in the data frame containing the length of the protein.
#' @param protein_id Optional argument, column in the data frame containing protein identifiers. Required if only one protein 
#' should be plotted and the data frame contains only information for this protein.
#' @param facet Optional argument, column in the data frame containing information by which data should be faceted. This can be 
#' protein identifiers. 
#' @param coloring Optional argument, column in the data frame containing information by which peptide or precursors should
#' be colored.
#' @param fold_change_cutoff Optional argument specifying the fold change cutoff used for assessing whether changes are significant. The default value is 2.
#'
#' @return A Wood's plot is returned.
#' @import ggplot2
#' @importFrom magrittr %>%
#' @importFrom dplyr pull
#' @importFrom rlang new_formula enquo as_name
#' @export
#'
#' @examples
#' \dontrun{
#' woods_plot(test,
#' fold_change = diff,
#' start_position = start,
#' end_position = end,
#' protein_length = length,
#' coloring = pep_type,
#' facet = pg_protein_accessions
#')
#' }
woods_plot <- function(data, fold_change, start_position, end_position, protein_length, protein_id = NULL, facet = NULL, coloring = NULL, fold_change_cutoff = 2){
  if(!missing(protein_id)){
    if(length(unique(dplyr::pull(data, {{protein_id}}))) > 1){
      stop("If data contains information of multiple proteins use the facet argument, not the protein_id argument")
    }
  }
  if(!missing(facet)){
    if(length(unique(dplyr::pull(data, {{facet}}))) > 20){
      n_proteins <- length(unique(dplyr::pull(data, {{facet}})))
      twenty_proteins <- unique(dplyr::pull(data, {{facet}}))[1:20]
      data <- data %>%
        filter({{facet}} %in% twenty_proteins)
      warning(paste("Only the first 20 proteins from", rlang::as_name(enquo(facet)), 
                    "have been used for plotting since there are", n_proteins, 
                    "proteins. Consider mapping over subsetted datasets." ))
    }
  }
  data %>%
    ggplot2::ggplot() +
    ggplot2::geom_rect(ggplot2::aes(xmin = {{fold_change}} - 0.5, xmax = {{fold_change}} + 0.5, ymin = {{start_position}}, ymax= {{end_position}}, fill = {{coloring}}), col = "black", size = 0.7) +
    ggplot2::geom_linerange(ggplot2::aes(x = 0, ymax = {{protein_length}}, ymin = 0)) +
    ggplot2::geom_vline(xintercept = -{{fold_change_cutoff}}, col = "blue", alpha = .8, size = 0.7) +
    ggplot2::geom_vline(xintercept = {{fold_change_cutoff}}, col = "blue", alpha = .8, size = 0.7) +
    ggplot2::coord_flip() +
    ggplot2::xlim(min(-2.5, dplyr::pull(data, {{fold_change}}))-0.5, max(2.5, dplyr::pull(data, {{fold_change}}))+0.5) +
    ggplot2::scale_y_continuous(limits = NULL, expand = c(0, 0)) +
    ggplot2::labs(x = "log2 Fold change", y = "Protein Sequence", title = {if(!missing(protein_id)) unique(dplyr::pull(data, {{protein_id}}))}) +
    {if(!missing(facet)) ggplot2::facet_wrap(rlang::new_formula(NULL, rlang::enquo(facet)), scales = "free", ncol = 4)} +
    ggplot2::guides(size = FALSE) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        size = 20
      ),
      axis.text.x = ggplot2::element_text(
        size = 15
      ),
      axis.text.y = ggplot2::element_text(
        size = 15
      ),
      axis.title.y = ggplot2::element_text(
        size = 15
      ),
      axis.title.x = ggplot2::element_text(
        size = 15
      ),
      legend.title = ggplot2::element_text(
        size = 15
      ),
      legend.text = ggplot2::element_text(
        size = 15
      ),
      strip.text.x = ggplot2::element_text(
        size = 15
      )
    )
}