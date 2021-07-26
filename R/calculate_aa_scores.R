#' Calculate amino acid scores
#' 
#' Calculate a score for each amino acid of a protein based on the product of the
#' -log10(adjusted p-value) and the absolute log2 fold change per peptide. In detail, all the peptides are 
#' aligned along the sequence of the corresponding protein, and the average score per amino acid is computed. 
#' In a limited proteolysis coupled to mass spectrometry (LiP-MS) experiment, the score allows to prioritize and 
#' narrow down structurally affected regions. 
#'
#' @param data a data frame containing at least the input columns.
#' @param adj_pval a numeric column in the data frame containing the adjusted p-value.
#' @param diff a numeric column in the data frame containing the log2 fold change.
#' @param grouping a character column in the data frame containing peptide or precursor identifiers.
#' @param start a numeric column in the data frame containing the start position of a peptide or precursor.
#' @param end a numeric column in the data frame containing the end position of a peptide or precursor.
#' @param protein a character column in the data frame containing the protein identifier or name.
#' 
#' @return a data frame that contains the aggregated scores per amino acid, enabling to draw amino acid fingerprints 
#' for each individual protein.
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
#'   grouping = pep_stripped_sequence,
#'   diff = diff,
#'   adj_pval = adj_pval,
#'   start = start,
#'   end = end
#' )
#' }


calculate_aa_scores <- function(data, protein, grouping, diff, adj_pval, start, end) {
  
data %>%
  dplyr::ungroup() %>% 
  dplyr::distinct({{ protein }}, {{ grouping }}, {{ diff }}, {{ adj_pval }}, {{ start }}, {{ end }}) %>% 
  dplyr::mutate(score = -log10({{adj_pval}}) * abs({{diff}})) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(residue = list(seq({{start}}, {{end}}))) %>% 
  tidyr::unnest(.data$residue) %>%
  dplyr::group_by({{ protein }}, .data$residue) %>% 
  dplyr::mutate(amino_acid_score = mean(.data$score)) %>%
  dplyr::distinct({{ protein }}, .data$residue, .data$amino_acid_score)
}