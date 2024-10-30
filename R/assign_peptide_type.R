#' Assign peptide type
#'
#' `r lifecycle::badge('deprecated')`
#' This function was deprecated due to its name changing to `assign_peptide_type()`.
#'
#' @return A data frame that contains the input data and an additional column with the peptide
#' type information.
#' @keywords internal
#' @export
peptide_type <- function(...) {
  # This function has been renamed and is therefore deprecated.
  lifecycle::deprecate_warn("0.2.0",
    "peptide_type()",
    "assign_peptide_type()",
    details = "This function has been renamed."
  )

  assign_peptide_type(...)
}
#' Assign peptide type
#'
#' Based on preceding and C-terminal amino acid, the peptide type of a given peptide is assigned.
#' Peptides with preceeding and C-terminal lysine or arginine are considered fully-tryptic. If a
#' peptide is located at the N- or C-terminus of a protein and fulfills the criterium to be
#' fully-tryptic otherwise, it is also considered as fully-tryptic. Peptides that only fulfill the
#' criterium on one terminus are semi-tryptic peptides. Lastly, peptides that are not fulfilling
#' the criteria for both termini are non-tryptic peptides. In addition, peptides that miss the initial
#' Methionine of a protein are considered "tryptic" at that site if there is no other peptide
#' starting at position 1 for that protein.
#'
#' @param data a data frame containing at least information about the preceding and C-terminal
#' amino acids of peptides.
#' @param aa_before a character column in the \code{data} data frame that contains the preceding amino
#' acid as one letter code.
#' @param last_aa a character column in the \code{data} data frame that contains the C-terminal amino
#' acid as one letter code.
#' @param aa_after a character column in the \code{data} data frame that contains the following amino
#' acid as one letter code.
#' @param protein_id a character column in the \code{data} data frame that contains the protein
#' accession numbers.
#' @param start a numeric column in the \code{data} data frame that contains the start position of
#' each peptide within the corresponding protein. This is used to check if the protein is consistently
#' missing the initial Methionine, making peptides starting at position 2 "tryptic" on that site.
#'
#' @return A data frame that contains the input data and an additional column with the peptide
#' type information.
#' @import dplyr
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @importFrom stringr str_detect
#' @export
#'
#' @examples
#' data <- data.frame(
#'   aa_before = c("K", "M", "", "M", "S", "M", "-"),
#'   last_aa = c("R", "K", "R", "R", "Y", "K", "K"),
#'   aa_after = c("T", "R", "T", "R", "T", "R", "T"),
#'   protein_id = c("P1", "P1", "P3", "P3", "P2", "P2", "P2"),
#'   start = c(38, 2, 1, 2, 10, 2, 1)
#' )
#'
#' assign_peptide_type(data, aa_before, last_aa, aa_after, protein_id, start)
assign_peptide_type <- function(data,
                                aa_before = aa_before,
                                last_aa = last_aa,
                                aa_after = aa_after,
                                protein_id = NULL,
                                start = start) {
  # Check if there's any peptide starting at position 1 for each protein
  start_summary <- data %>%
    dplyr::group_by({{ protein_id }}) %>%
    dplyr::summarize(has_start_1 = any({{ start }} == 1), .groups = "drop")

  peptide_data <- data %>%
    dplyr::distinct({{ aa_before }}, {{ last_aa }}, {{ aa_after }}, {{ protein_id }}, {{ start }}, .keep_all = TRUE) %>%
    dplyr::left_join(start_summary, by = rlang::as_name(rlang::enquo(protein_id))) %>%
    # Determine N-terminal trypticity
    dplyr::mutate(N_term_tryp = dplyr::if_else(
      !stringr::str_detect({{ aa_before }}, "[A-Y]") | {{ aa_before }} == "K" | {{ aa_before }} == "R",
      TRUE,
      FALSE
    )) %>%
    # Determine C-terminal trypticity
    dplyr::mutate(C_term_tryp = dplyr::if_else(
      {{ last_aa }} == "K" | {{ last_aa }} == "R" | !stringr::str_detect({{ aa_after }}, "[A-Y]"),
      TRUE,
      FALSE
    )) %>%
    # Reassign peptides to be N_term_tryp if the protein starts on position 2
    dplyr::mutate(N_term_tryp = ifelse({{ start }} == 2 & !.data$has_start_1, TRUE, .data$N_term_tryp)) %>%
    # Assign peptide type based on N-term and C-term trypticity
    dplyr::mutate(pep_type = dplyr::case_when(
      .data$N_term_tryp & .data$C_term_tryp ~ "fully-tryptic",
      .data$N_term_tryp | .data$C_term_tryp ~ "semi-tryptic",
      TRUE ~ "non-tryptic"
    )) %>%
    # Drop unnecessary columns
    dplyr::select(-c("N_term_tryp", "C_term_tryp", "has_start_1"))

  # Join back to original data to return the full result
  result <- data %>%
    dplyr::left_join(
      peptide_data %>%
        dplyr::select({{ aa_before }}, {{ last_aa }}, {{ aa_after }}, {{ protein_id }}, {{ start }}, "pep_type"),
      by = c(
        rlang::as_name(rlang::enquo(aa_before)),
        rlang::as_name(rlang::enquo(last_aa)),
        rlang::as_name(rlang::enquo(aa_after)),
        rlang::as_name(rlang::enquo(protein_id)),
        rlang::as_name(rlang::enquo(start))
      )
    )

  return(result)
}
