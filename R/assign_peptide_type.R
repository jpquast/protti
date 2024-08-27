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
#' the criteria for both termini are non-tryptic peptides.
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
#'
#' @return A data frame that contains the input data and an additional column with the peptide
#' type information.
#' @import dplyr
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @export
#'
#' @examples
#' data <- data.frame(
#'   aa_before = c("K", "S", "T", "M"),
#'   last_aa = c("R", "K", "Y", "K"),
#'   aa_after = c("T", "R", "T", "R"),
#'   protein_id = c("P1", "P1", "P2", "P2")
#' )
#'
#' assign_peptide_type(data, aa_before, last_aa, aa_after)
assign_peptide_type <- function(data,
                                aa_before = aa_before,
                                last_aa = last_aa,
                                aa_after = aa_after,
                                protein_id = protein_id) {
  summary_data <- data %>%
    dplyr::group_by({{ protein_id }}) %>%
    dplyr::summarize(missing_methionine = !any({{ aa_before }} == "M"), .groups = "drop")

  peptide_data <- data %>%
    dplyr::distinct({{ aa_before }}, {{ last_aa }}, {{ aa_after }}, {{ protein_id }}, .keep_all = TRUE) %>%
    dplyr::left_join(summary_data, by = rlang::as_name(rlang::enquo(protein_id))) %>%
    dplyr::mutate(N_term_tryp = dplyr::if_else(
      {{ aa_before }} == "" |
        {{ aa_before }} == "K" |
        {{ aa_before }} == "R",
      TRUE,
      FALSE
    )) %>%
    dplyr::mutate(C_term_tryp = dplyr::if_else(
      {{ last_aa }} == "K" |
        {{ last_aa }} == "R" |
        {{ aa_after }} == "",
      TRUE,
      FALSE
    )) %>%
    dplyr::mutate(pep_type = dplyr::case_when(
      N_term_tryp + C_term_tryp == 2 ~ "fully-tryptic",
      N_term_tryp + C_term_tryp == 1 ~ "semi-tryptic",
      N_term_tryp + C_term_tryp == 0 ~ "non-tryptic"
    )) %>%
    # Check if initial methionine is missing and reassign semi-tryptic to fully-tryptic
    dplyr::mutate(pep_type = dplyr::if_else(
      pep_type == "semi-tryptic" & missing_methionine,
      "fully-tryptic",
      pep_type
    )) %>%
    dplyr::select(-N_term_tryp, -C_term_tryp, -missing_methionine)

  result <- data %>%
    dplyr::left_join(
      peptide_data %>%
        dplyr::select({{ aa_before }}, {{ last_aa }}, {{ aa_after }}, {{ protein_id }}, pep_type),
      by = c(
        rlang::as_name(rlang::enquo(aa_before)),
        rlang::as_name(rlang::enquo(last_aa)),
        rlang::as_name(rlang::enquo(aa_after)),
        rlang::as_name(rlang::enquo(protein_id))
      )
    )
  return(result)
}
