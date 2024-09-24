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
#' @param start_pos A numeric column in the \code{data} data frame that contains the start position of
#' each peptide within the corresponding protein. This is used to check if the peptide starts at position 1
#' or position 2, which affects whether the peptide can be considered fully-tryptic.
#'
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
#' aa_before = c("K", "S", "K", "S", "T", "M"),
#' last_aa = c("R", "K","R", "K", "Y", "K"),
#' aa_after = c("T", "R", "T", "R", "T", "R"),
#' protein_id = c("P1", "P1", "P3", "P3","P2", "P2"),
#' start_pos = c(2, 2, 1, 2, 1, 2),
#' )
#'
#' assign_peptide_type(data, aa_before, last_aa, aa_after, protein_id, start_pos)
assign_peptide_type <- function(data,
                                aa_before = aa_before,
                                last_aa = last_aa,
                                aa_after = aa_after,
                                protein_id = protein_id,
                                start_pos = start_pos) {
  # Check if there's any peptide starting at position 1 for each protein
  start_summary <- data %>%
    dplyr::group_by({{ protein_id }}) %>%
    dplyr::summarize(has_start_pos_1 = any({{ start_pos }} == 1), .groups = "drop")

  peptide_data <- data %>%
    dplyr::distinct({{ aa_before }}, {{ last_aa }}, {{ aa_after }}, {{ protein_id }}, {{ start_pos }}, .keep_all = TRUE) %>%
    dplyr::left_join(start_summary, by = rlang::as_name(rlang::enquo(protein_id))) %>%

    # Determine N-terminal trypticity
    dplyr::mutate(N_term_tryp = dplyr::if_else(
      {{ aa_before }} == "" | {{ aa_before }} == "K" | {{ aa_before }} == "R",
      TRUE,
      FALSE
    )) %>%

    # Determine C-terminal trypticity
    dplyr::mutate(C_term_tryp = dplyr::if_else(
      {{ last_aa }} == "K" | {{ last_aa }} == "R" | {{ aa_after }} == "",
      TRUE,
      FALSE
    )) %>%

    # Assign peptide type based on N-term and C-term trypticity
    dplyr::mutate(pep_type = dplyr::case_when(
      N_term_tryp & C_term_tryp ~ "fully-tryptic",
      N_term_tryp | C_term_tryp ~ "semi-tryptic",
      TRUE ~ "non-tryptic"
    )) %>%

    # Reassign semi-tryptic peptides at position 2 to fully-tryptic if no start_pos == 1
    dplyr::mutate(pep_type = dplyr::if_else(
      pep_type == "semi-tryptic" & {{ start_pos }} == 2 & !.data$has_start_pos_1 & C_term_tryp,
      "fully-tryptic",
      .data$pep_type
    )) %>%

    # Drop unnecessary columns
    dplyr::select(-N_term_tryp, -C_term_tryp, -has_start_pos_1)

  # Join back to original data to return the full result
  result <- data %>%
    dplyr::left_join(
      peptide_data %>%
        dplyr::select({{ aa_before }}, {{ last_aa }}, {{ aa_after }}, {{ protein_id }}, {{ start_pos }}, pep_type),
      by = c(
        rlang::as_name(rlang::enquo(aa_before)),
        rlang::as_name(rlang::enquo(last_aa)),
        rlang::as_name(rlang::enquo(aa_after)),
        rlang::as_name(rlang::enquo(protein_id)),
        rlang::as_name(rlang::enquo(start_pos))
      )
    )

  return(result)
}
