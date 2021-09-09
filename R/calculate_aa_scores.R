#' Calculate scores for each amino acid position in a protein sequence
#'
#' `r lifecycle::badge("experimental")`
#' Calculate a score for each amino acid position in a protein sequence based on the product of the
#' -log10(adjusted p-value) and the absolute log2(fold change) per peptide covering this amino acid. In detail, all the
#' peptides are aligned along the sequence of the corresponding protein, and the average score per
#' amino acid position is computed. In a limited proteolysis coupled to mass spectrometry (LiP-MS)
#' experiment, the score allows to prioritize and narrow down structurally affected regions.
#'
#' @param data a data frame containing at least the input columns.
#' @param adj_pval a numeric column in the \code{data} data frame containing the adjusted p-value.
#' @param diff a numeric column in the \code{data} data frame containing the log2 fold change.
#' @param start_position a numeric column \code{data} in the data frame containing the start position
#' of a peptide or precursor.
#' @param end_position a numeric column in the data frame containing the end position of a peptide or
#' precursor.
#' @param protein a character column in the data frame containing the protein identifier or name.
#' @param retain_columns a vector indicating if certain columns should be retained from the input
#' data frame. Default is not retaining additional columns \code{retain_columns = NULL}. Specific
#' columns can be retained by providing their names (not in quotations marks, just like other
#' column names, but in a vector).
#'
#' @return A data frame that contains the aggregated scores per amino acid position, enabling to
#' draw fingerprints for each individual protein.
#'
#' @author Patrick Stalder
#' @import dplyr
#' @import tidyr
#' @export
#'
#' @examples
#' \dontrun{
#' calculate_aa_scores(
#'   data,
#'   protein = pg_protein_accessions,
#'   diff = diff,
#'   adj_pval = adj_pval,
#'   start_position = start,
#'   end_position = end
#' )
#' }
#'
calculate_aa_scores <- function(data,
                                protein,
                                diff = diff,
                                adj_pval = adj_pval,
                                start_position,
                                end_position,
                                retain_columns = NULL) {
  output <- data %>%
    dplyr::ungroup() %>%
    dplyr::distinct({{ protein }}, {{ diff }}, {{ adj_pval }}, {{ start_position }}, {{ end_position }}) %>%
    tidyr::drop_na({{ diff }}, {{ adj_pval }}) %>% 
    dplyr::mutate(score = -log10({{ adj_pval }}) * abs({{ diff }})) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(residue = list(seq({{ start_position }}, {{ end_position }}))) %>%
    tidyr::unnest(.data$residue) %>%
    dplyr::group_by({{ protein }}, .data$residue) %>%
    dplyr::mutate(amino_acid_score = mean(.data$score)) %>%
    dplyr::distinct({{ protein }}, .data$residue, .data$amino_acid_score) 
  
  
  if (!missing(retain_columns)) {
    output <- data %>%
      dplyr::select(!!enquo(retain_columns), colnames(output)[!colnames(output) %in% c(
        "residue",
        "amino_acid_score"
      )]) %>%
      dplyr::distinct() %>%
      dplyr::right_join(output, by = colnames(output)[!colnames(output) %in% c(
        "residue",
        "amino_acid_score"
      )])
  }
  
  output
}
