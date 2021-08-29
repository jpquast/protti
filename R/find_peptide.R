#' Find peptide location
#'
#' The position of the given peptide sequence is searched within the given protein sequence. In
#' addition the last amino acid of the peptide and the amino acid right before are reported.
#'
#' @param data a data frame that contains at least the protein and peptide sequence.
#' @param protein_sequence a character column in the \code{data} data frame that contains the
#' protein sequence.
#' @param peptide_sequence a character column in the \code{data} data frame that contains the
#' peptide sequence.
#'
#' @return A data frame that contains the input data and four additional columns with peptide
#' start and end position, the last amino acid and the amino acid before the peptide.
#' @import dplyr
#' @import stringr
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @export
#'
#' @examples
#' # Create example data
#' data <- data.frame(
#'   protein_sequence = c("abcdefg"),
#'   peptide_sequence = c("cde")
#' )
#'
#' # Find peptide
#' find_peptide(
#'   data = data,
#'   protein_sequence = protein_sequence,
#'   peptide_sequence = peptide_sequence
#' )
find_peptide <-
  function(data, protein_sequence, peptide_sequence) {
    data %>%
      dplyr::mutate(
        start = stringr::str_locate({{ protein_sequence }}, {{ peptide_sequence }})[, 1],
        end = stringr::str_locate({{ protein_sequence }}, {{ peptide_sequence }})[, 2]
      ) %>%
      dplyr::mutate(aa_before = stringr::str_sub({{ protein_sequence }},
        start = .data$start - 1,
        end = .data$start - 1
      )) %>%
      dplyr::mutate(last_aa = stringr::str_sub({{ protein_sequence }},
        start = .data$end,
        end = .data$end
      )) %>%
      dplyr::mutate(aa_after = stringr::str_sub({{ protein_sequence }},
        start = .data$end + 1,
        end = .data$end + 1
      ))
  }
