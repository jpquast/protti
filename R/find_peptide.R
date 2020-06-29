#' Find peptide location
#'
#' The position of the given peptide sequence is searched within the given protein sequence. In addition the last amino acid of the peptide and the amino acid right before are reported.
#'
#' @param data A data frame containing at least the protein and peptide sequence.
#' @param protein_sequence The name of the column containing the protein sequence.
#' @param peptide_sequence The name of the column containing the peptide sequence.
#'
#' @return A data frame that contains the input data and four additional columns with start and end position, the last amino acid and the amino acid before.
#' @import dplyr
#' @import stringr
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @export
#'
#' @examples
#' \dontrun{
#' find_peptide(data, protein_sequence, pep_stripped_sequence)
#' }
find_peptide <- 
  function(data, protein_sequence, peptide_sequence){
  data%>%
    mutate(start = str_locate({{protein_sequence}}, {{peptide_sequence}})[,1], end = str_locate({{protein_sequence}}, {{peptide_sequence}})[,2])%>%
    mutate(aa_before = str_sub({{protein_sequence}}, start = .data$start-1, end = .data$start-1))%>%
    mutate(last_aa = str_sub({{protein_sequence}}, start = .data$end, end = .data$end))
  }
