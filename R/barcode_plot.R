#' Barcode plot
#'
#' Plots a "barcode plot" - a line for each identified peptide at the identified position with peptides changing in abundance highlighted.
#'
#' @param data Dataframe containing at least the input variables.
#' @param protein_id Column in the data dataframe containing the protein IDs.
#' @param target Protein ID for which a barcode plot should be plotted.
#' @param peptides Column containing the sequences.
#' @param start Start position of the peptide in the protein sequence.
#' @param end End position of the peptide in the protein sequence.
#' @param protein_length Length of the protein.
#' @param foldchange Fold change values for each peptide.
#' @param foldchange_cutoff Optional argument specifying the cutoff for highlighted peptides in the plot, default = 1.
#' 
#' @return A barcode plot for the target protein.
#' @import dplyr
#' @import ggplot2
#' @importFrom magrittr %>%
#' @importFrom rlang as_name
#' @importFrom rlang :=
#' @importFrom rlang enquo
#' @export
#'
#' @examples
#' \dontrun{
#' barcode_plot(
#' data,
#' protein_id = uniprot_id,
#' target = "Q9Y6K9",
#' peptides = pep_stripped_sequence,
#' start = start, 
#' end = end,
#' protein_length = length_protein,
#' foldchange = Log2FC,
#' foldchange_cutoff = 0.5
#' )
#' }
barcode_plot <- function(data, protein_id, target, peptides, start, end, protein_length, foldchange, foldchange_cutoff = 1)
{
  result <- data %>%
    filter({{protein_id}} == target)
  
  plot <- result %>%
    ggplot2::ggplot() +
    ggplot2::geom_rect(aes(
      ymin = -2.5,
      ymax = 2.5,
      xmax = (({{start}} + nchar({{peptides}})) / {{protein_length}}) * 100,
      xmin = ({{start}} / {{protein_length}}) * 100
    )) +
    ggplot2::geom_rect(
      data = dplyr::filter(result, {{foldchange}} < -1 * foldchange_cutoff | {{foldchange}} > foldchange_cutoff),
      aes(
        ymin = -2.5,
        ymax = 2.5,
        xmax = (({{start}} + nchar({{peptides}})) / {{protein_length}}) * 100,
        xmin = ({{start}} / {{protein_length}}) * 100),
      fill = "cornflowerblue"
    ) +
    ggplot2::theme(axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          legend.position = 0, 
          panel.background = element_blank(), 
          panel.border = element_rect(fill = NA)
    ) +
    ggplot2::xlab("Protein Sequence") +
    ggplot2::labs(title = rlang::as_name(rlang::enquo(target)))
 return(plot)
}

