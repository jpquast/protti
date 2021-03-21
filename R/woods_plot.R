#' Wood's plot
#'
#' Creates a Wood's plot that plots log2 fold change of peptides or precursors along the protein sequence.
#'
#' @param data Data frame containing differential abundance, start and end peptide or precursor positions, protein
#' length and optionally a variable based on which peptides or precursors should be coloured.
#' @param fold_change Column in the data frame containing log2 fold changes.
#' @param start_position Column in the data frame containing the start positions for each peptide or precursor.
#' @param end_position Column in the data frame containing the end positions for each peptide or precursor.
#' @param protein_length Column in the data frame containing the length of the protein.
#' @param coverage Optional, column in the data frame containing coverage in percent. Will appear in the title of the barcode if provided.
#' @param protein_id Optional argument, column in the data frame containing protein identifiers. Required if only one protein
#' should be plotted and the data frame contains only information for this protein.
#' @param facet Optional argument, column in the data frame containing information by which data should be faceted. This can be
#' protein identifiers.
#' @param colouring Optional argument, column in the data frame containing information by which peptide or precursors should
#' be coloured.
#' @param fold_change_cutoff Optional argument specifying the log2 fold change cutoff used for assessing whether changes are significant. The default value is 2.
#'
#' @return A Wood's plot is returned. Plotting peptide or precursor fold changes accross protein sequence.
#' @import ggplot2
#' @importFrom magrittr %>%
#' @importFrom dplyr pull mutate filter
#' @importFrom rlang new_formula enquo as_name :=
#' @importFrom utils data
#' @export
#'
#' @examples
#' \dontrun{
#' woods_plot(test,
#'   fold_change = diff,
#'   start_position = start,
#'   end_position = end,
#'   protein_length = length,
#'   colouring = pep_type,
#'   facet = pg_protein_accessions
#' )
#' }
woods_plot <- function(data, fold_change, start_position, end_position, protein_length, coverage = NULL, protein_id = NULL, facet = NULL, colouring = NULL, fold_change_cutoff = 1) {
  protti_colours <- "placeholder" # assign a placeholder to prevent a missing global variable warning
  utils::data("protti_colours", envir = environment()) # then overwrite it with real data
  # Check if there are more than one protein even though protein_id was specified.
  if (!missing(protein_id)) {
    if (length(unique(dplyr::pull(data, {{ protein_id }}))) > 1) {
      stop("If data contains information of multiple proteins use the facet argument, not the protein_id argument")
    }
  }
  # Check if there are more than 20 proteins for faceting.
  if (!missing(facet)) {
    if (length(unique(dplyr::pull(data, {{ facet }}))) > 20) {
      n_proteins <- length(unique(dplyr::pull(data, {{ facet }})))
      twenty_proteins <- unique(dplyr::pull(data, {{ facet }}))[1:20]
      data <- data %>%
        dplyr::filter({{ facet }} %in% twenty_proteins)
      warning(paste(
        "Only the first 20 proteins from", rlang::as_name(enquo(facet)),
        "have been used for plotting since there are", n_proteins,
        "proteins. Consider mapping over subsetted datasets."
      ))
    }
  }
  # Add coverage to protein ID name if present.
  if (!missing(coverage) & !missing(facet)) {
    data <- data %>%
      dplyr::mutate({{ facet }} := paste0({{ facet }}, " (", round({{ coverage }}, digits = 1), "%)"))
  }
  if (!missing(coverage) & !missing(protein_id)) {
    data <- data %>%
      dplyr::mutate({{ protein_id }} := paste0({{ protein_id }}, " (", round({{ coverage }}, digits = 1), "%)"))
  }
  # Create plot
  data %>%
    ggplot2::ggplot() +
    ggplot2::geom_rect(ggplot2::aes(
      xmin = 0,
      xmax = {{ protein_length }},
      ymin = -0.01,
      ymax = 0.01
    ),
    fill = "black"
    ) +
    ggplot2::geom_rect(ggplot2::aes(
      xmin = {{ start_position }},
      xmax = {{ end_position }},
      ymin = {{ fold_change }} - 0.5,
      ymax = {{ fold_change }} + 0.5,
      fill = {{ colouring }}
    ),
    col = "black",
    size = 0.7
    ) +
    ggplot2::geom_hline(
      yintercept = -{{ fold_change_cutoff }},
      col = "blue",
      alpha = .8,
      size = 0.7
    ) +
    ggplot2::geom_hline(
      yintercept = {{ fold_change_cutoff }},
      col = "blue",
      alpha = .8,
      size = 0.7
    ) +
    ggplot2::ylim(min(-2.5, dplyr::pull(data, {{ fold_change }})) - 0.5, max(2.5, dplyr::pull(data, {{ fold_change }})) + 0.5) +
    ggplot2::scale_x_continuous(limits = NULL, expand = c(0, 0)) +
    ggplot2::labs(
      x = "Protein Sequence",
      y = "log2(fold change)",
      title = {
        if (!missing(protein_id)) unique(dplyr::pull(data, {{ protein_id }}))
      }
    ) +
    {
      if (!missing(facet)) ggplot2::facet_wrap(rlang::new_formula(NULL, rlang::enquo(facet)), scales = "free", ncol = 4)
    } +
    ggplot2::guides(size = FALSE) +
    {
      if (!missing(colouring) && !is.numeric(dplyr::pull(data, {{ colouring }}))) ggplot2::scale_fill_manual(values = protti_colours)
    } +
    {
      if (!missing(colouring) && is.numeric(dplyr::pull(data, {{ colouring }}))) ggplot2::scale_fill_gradient(low = protti_colours[1], high = protti_colours[2])
    } +
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
      ),
      strip.text = ggplot2::element_text(
        size = 15
      ),
      strip.background = element_blank()
    )
}
