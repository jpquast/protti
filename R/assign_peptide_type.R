#' Assign peptide type
#'
#' `r lifecycle::badge('deprecated')`
#' This function was deprecated due to its name changing to `assign_peptide_type()`.
#'
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
#'   aa_before = c("K", "S", "T"),
#'   last_aa = c("R", "K", "Y"),
#'   aa_after = c("T", "R", "T")
#' )
#'
#' assign_peptide_type(data, aa_before, last_aa, aa_after)
assign_peptide_type <- function(data,
                                aa_before = aa_before,
                                last_aa = last_aa,
                                aa_after = aa_after) {
  data %>%
    dplyr::mutate(N_term_tryp = dplyr::if_else({{ aa_before }} == "" |
      {{ aa_before }} == "K" |
      {{ aa_before }} == "R",
    TRUE,
    FALSE
    )) %>%
    dplyr::mutate(C_term_tryp = dplyr::if_else({{ last_aa }} == "K" |
      {{ last_aa }} == "R" |
      {{ aa_after }} == "",
    TRUE,
    FALSE
    )) %>%
    dplyr::mutate(pep_type = dplyr::case_when(
      .data$N_term_tryp + .data$C_term_tryp == 2 ~ "fully-tryptic",
      .data$N_term_tryp + .data$C_term_tryp == 1 ~ "semi-tryptic",
      .data$N_term_tryp + .data$C_term_tryp == 0 ~ "non-tryptic"
    )) %>%
    dplyr::select(-.data$N_term_tryp, -.data$C_term_tryp)
}
