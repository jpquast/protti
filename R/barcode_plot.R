#' Barcode plot
#'
#' Plots a "barcode plot" - a vertical line for each identified peptide. Peptides can be colored based on an additional variable. Also differential
#' abundance can be displayed.
#'
#' @param data Data frame containing differential abundance, start and end peptide or precursor positions and protein length.
#' @param start_position Column in the data frame containing the start positions for each peptide or precursor.
#' @param end_position Column in the data frame containing the end positions for each peptide or precursor.
#' @param protein_length Column in the data frame containing the length of the protein.
#' @param coverage Optional, column in the data frame containing coverage in percent. Will appear in the title of the barcode if provided.
#' @param colouring Optional argument, column in the data frame containing information by which peptide or precursors should
#' be colored.
#' @param protein_id Optional argument, column in the data frame containing protein identifiers. Required if only one protein
#' should be plotted and the data frame contains only information for this protein.
#' @param facet Optional argument, column in the data frame containing information by which data should be faceted. This can be
#' protein identifiers. Only 20 proteins are plotted at a time, the rest is ignored. If more should be plotted, a mapper over a
#' subsetted data frame should be created.
#' @param cutoffs Optional argument specifying the log2 fold change and significance cutoffs used for highlighting peptides.
#' If this argument is provided colouring information will be overwritten with peptides that fulfill this condition.
#' The cutoff should be provided in a vector of the form c(diff = 2, pval = 0.05). The name of the cutoff should reflect the 
#' column name that contains this information (log2 fold changes, p-values or adjusted p-values). 
#'
#' @return A barcode plot is returned.
#' @import dplyr
#' @import ggplot2
#' @importFrom magrittr %>%
#' @importFrom rlang as_name := enquo new_formula !! ensym
#' @importFrom forcats fct_rev
#' @export
#'
#' @examples
#' \dontrun{
#' barcode_plot(
#' data,
#' start_position = start,
#' end_position = end,
#' protein_length = length,
#' facet = pg_protein_accessions,
#' cutoffs = c(diff = 2, pval = 0.05)
#' )
#' }
barcode_plot <- function(data, start_position, end_position, protein_length, coverage = NULL, colouring = NULL, protein_id = NULL, facet = NULL, cutoffs = NULL)
{
  # Check if there is more than one protein even though protein_id was specified.
  if(!missing(protein_id)){
    if(length(unique(dplyr::pull(data, {{protein_id}}))) > 1){
      stop("If data contains information of multiple proteins use the facet argument, not the protein_id argument")
    }
  }
  # Check if there are more than 20 proteins for faceting.
  if(!missing(facet)){
    if(length(unique(dplyr::pull(data, {{facet}}))) > 20){
      n_proteins <- length(unique(dplyr::pull(data, {{facet}})))
      twenty_proteins <- unique(dplyr::pull(data, {{facet}}))[1:20]
      data <- data %>%
        dplyr::filter({{facet}} %in% twenty_proteins)
      warning(paste("Only the first 20 proteins from", rlang::as_name(enquo(facet)),
                    "have been used for plotting since there are", n_proteins,
                    "proteins. Consider mapping over subsetted datasets." ))
    }
  }
  # Apply fold change  and significance cutoff if fold change is provided
  if(!missing(cutoffs)){
    fc_name <- names(cutoffs)[1]
    sig_name <- names(cutoffs)[2]
    fc <- cutoffs[1]
    sig <- cutoffs[2]
    
    colouring = sym("change")
    
    data <- data %>%
      dplyr::mutate({{colouring}} := ifelse(((!!ensym(fc_name) >= fc | !!ensym(fc_name) <= -fc) & !!ensym(sig_name) <= sig), "Structural change", "Unchanged")) %>%
      dplyr::mutate({{colouring}} := forcats::fct_rev(ifelse(is.na({{colouring}}), "Unchanged", {{colouring}}))) %>%
      dplyr::arrange({{colouring}})
  }
  # Add coverage to protein ID name if present.
  if(!missing(coverage) & !missing(facet)){
    data <- data %>%
      mutate({{facet}} := paste0({{facet}}, " (", round({{coverage}}, digits = 1), "%)"))
  }
  if(!missing(coverage) & !missing(protein_id)){
    data <- data %>%
      mutate({{protein_id}} := paste0({{protein_id}}, " (", round({{coverage}}, digits = 1), "%)"))
  }
  # Create plot
  data %>%
    ggplot2::ggplot() +
    ggplot2::geom_rect(ggplot2::aes(
      ymin = -2.5,
      ymax = 2.5,
      xmax = {{end_position}} / {{protein_length}} * 100,
      xmin = {{start_position}} / {{protein_length}} * 100,
      fill = {{colouring}}),
      size = 0.7
    ) +
    ggplot2::scale_fill_manual(values = c("#999999", "#5680C1", "#B96DAD", "#64CACA", "#81ABE9", "#F6B8D1", "#99F1E4", "#9AD1FF", "#548BDF", "#A55098", "#3EB6B6",
                                          "#87AEE8", "#CA91C1", "#A4E0E0", "#1D4F9A", "#D7ACD2","#49C1C1")) +
    ggplot2::scale_x_continuous(limits = c(0, 100), expand = c(0, 0)) +
    ggplot2::scale_y_continuous(limits = NULL, expand = c(0, 0)) +
    ggplot2::labs(x = "Protein Sequence", title = {if(!missing(protein_id)) unique(dplyr::pull(data, {{protein_id}}))}) +
    {if(!missing(facet)) ggplot2::facet_wrap(rlang::new_formula(NULL, rlang::enquo(facet)))} +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 20),
                   axis.title.y = element_blank(),
                   axis.text.y = element_blank(),
                   axis.ticks.y = element_blank(),
                   axis.text.x = element_blank(),
                   axis.title.x = ggplot2::element_text(size = 15),
                   axis.ticks.x = element_blank(),
                   legend.title = ggplot2::element_text(size = 15),
                   legend.text = ggplot2::element_text(size = 15),
                   strip.text = ggplot2::element_text(size = 15),
                   strip.background = element_blank(),
                   panel.background = element_blank(),
                   panel.border = element_rect(fill = NA)
    )
}
