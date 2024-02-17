#' Protein sequence coverage
#'
#' `r lifecycle::badge('deprecated')`
#' This function was deprecated due to its name changing to `calculate_sequence_coverage()`.
#'
#' @return A new column in the \code{data} data frame containing the calculated sequence coverage
#' for each identified protein
#' @keywords internal
#' @export
sequence_coverage <- function(...) {
  # This function has been renamed and is therefore deprecated.
  lifecycle::deprecate_warn("0.2.0",
    "sequence_coverage()",
    "calculate_sequence_coverage()",
    details = "This function has been renamed."
  )
  calculate_sequence_coverage(...)
}
#' Protein sequence coverage
#'
#' Calculate sequence coverage for each identified protein.
#'
#' @param data a data frame containing at least the protein sequence and the identified peptides
#' as columns.
#' @param protein_sequence a character column in the \code{data} data frame that contains protein
#' sequences. Can be obtained by using the function \code{fetch_uniprot()}
#' @param peptides a character column in the \code{data} data frame that contains the identified
#' peptides.
#'
#' @return A new column in the \code{data} data frame containing the calculated sequence coverage
#' for each identified protein
#' @import dplyr
#' @importFrom magrittr %>%
#' @importFrom stringr str_count
#' @importFrom rlang .data as_name enquo
#' @importFrom tidyr drop_na
#' @export
#'
#' @examples
#' data <- data.frame(
#'   protein_sequence = c("abcdefghijklmnop", "abcdefghijklmnop"),
#'   pep_stripped_sequence = c("abc", "jklmn")
#' )
#'
#' calculate_sequence_coverage(
#'   data,
#'   protein_sequence = protein_sequence,
#'   peptides = pep_stripped_sequence
#' )
calculate_sequence_coverage <-
  function(data, protein_sequence, peptides) {
    result <- data %>%
      dplyr::ungroup() %>%
      # drop_na prevents function from failing if a protein group contains only NA peptide sequences.
      tidyr::drop_na({{ peptides }}) %>%
      dplyr::distinct({{ protein_sequence }}, {{ peptides }}) %>%
      dplyr::group_by({{ protein_sequence }}) %>%
      find_peptide({{ protein_sequence }}, {{ peptides }}) %>%
      dplyr::mutate(sequence_length = nchar({{ protein_sequence }})) %>%
      dplyr::mutate(modified_sequence = replace_identified_by_x({{ protein_sequence }}, .data$start, .data$end)) %>%
      dplyr::mutate(covered = stringr::str_count(.data$modified_sequence, "x")) %>%
      dplyr::mutate(coverage = .data$covered / .data$sequence_length * 100) %>%
      dplyr::select(-c(
        .data$sequence_length,
        .data$modified_sequence,
        .data$covered,
        .data$start,
        .data$end,
        .data$aa_before,
        .data$last_aa,
        .data$aa_after,
        {{ peptides }}
      )) %>%
      dplyr::distinct() %>%
      dplyr::ungroup()


    result <- data %>%
      dplyr::left_join(result, by = rlang::as_name(rlang::enquo(protein_sequence)))

    result
  }
