#' Protein sequence coverage
#'
#' Calculate sequence coverage for each identified protein.
#'
#' @param data A dataframe containing at least the protein sequence and the identified peptides as columns.
#' @param protein_sequence A column containing protein sequences, can be obtained by using the function \code{fetch_uniprot()}
#' @param peptides A column containing the identified peptides.
#'
#' @return A new column containing the calculated sequence coverages for each identified protein
#' @import dplyr
#' @importFrom magrittr %>%
#' @importFrom stringr str_count
#' @importFrom rlang .data as_name enquo
#' @export
#'
#' @examples
#' \dontrun{
#' sequence_coverage(
#' data,
#' protein_sequence = protein_sequence,
#' peptides = pep_stripped_sequence
#' )
#' }
#'
sequence_coverage <-
  function(data, protein_sequence, peptides)
  {
    result <- data %>%
      dplyr::distinct({{protein_sequence}}, {{peptides}}) %>%
      dplyr::group_by({{protein_sequence}}) %>%
      find_peptide({{protein_sequence}}, {{peptides}}) %>%
      dplyr::mutate(sequence_length = nchar({{protein_sequence}})) %>%
      dplyr::mutate(modified_sequence = replace_identified_by_x({{protein_sequence}}, .data$start, .data$end)) %>%
      dplyr::mutate(covered = stringr::str_count(.data$modified_sequence, "x")) %>%
      dplyr::mutate(coverage = .data$covered / .data$sequence_length * 100) %>%
      dplyr::select(-c(.data$sequence_length, .data$modified_sequence, .data$covered, .data$start, .data$end, .data$aa_before, .data$last_aa, .data$aa_after, {{peptides}})) %>%
      dplyr::distinct() %>%
      dplyr::ungroup()
    
    
    result <- data %>%
      dplyr::left_join(result, by = rlang::as_name(rlang::enquo(protein_sequence)))
    
    result
  }