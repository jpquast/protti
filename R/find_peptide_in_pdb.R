#' Finds peptide positions in a PDB structure based on positional matching
#'
#' Finds peptide positions in a PDB structure. Often positions of peptides in UniProt and a PDB structure are different due to different
#' lengths of structures. This function maps a peptide based on its UniProt positions onto a PDB structure. This method is superior to
#' sequence alignment of the peptide to the PDB structure sequence, since it can also match the peptide if there are truncations or
#' mismatches. This function also provides an easy way to check if a peptide is present in a PDB structure.
#'
#' @param peptide_data a data frame containing at least the input columns to this function.
#' @param peptide_sequence a character column in the \code{peptide_data} data frame that contains the sequence or any other unique identifier
#' for the peptide that should be found.
#' @param start a numeric column in the \code{peptide_data} data frame that contains start positions of peptides.
#' @param end a numeric column in the \code{peptide_data} data frame that contains end positions of peptides.
#' @param uniprot_id a character column in the \code{peptide_data} data frame that contains UniProt identifiers that correspond to the
#' peptides.
#' @param pdb_data optional data frame containing data obtained with \code{fetch_pdb()}. If not provided, information is fetched automatically.
#' If this function should be run multiple times it is faster to fetch the information once and provide it to the function. If provided,
#' make sure that the column names are identical to the ones that would be obtained by calling \code{fetch_pdb()}.
#' @param retain_columns a vector indicating if certain columns should be retained from the input data frame. Default is not retaining
#' additional columns \code{retain_columns = NULL}. Specific columns can be retained by providing their names (not in quotations marks,
#' just like other column names, but in a vector).
#'
#' @return A data frame that contains peptide positions in the corresponding PDB structures. If a peptide is not found in any structure,
#' it contains NAs values for the corresponding positional columns.
#' @import dplyr
#' @import tidyr
#' @importFrom stringr str_sub
#' @importFrom magrittr %>%
#' @importFrom rlang as_name enquo .data
#' @export
#'
#' @examples
#' \dontrun{
#' find_peptide_in_pdb(
#'   peptide_data = peptide_data,
#'   peptide_sequence = peptide,
#'   start = start,
#'   end = end,
#'   uniprot_id = pg_protein_accessions,
#'   pdb_data = pdb_data,
#'   retain_columns = c(map_value)
#' )
#' }
find_peptide_in_pdb <- function(peptide_data, peptide_sequence, start, end, uniprot_id, pdb_data = NULL, retain_columns = NULL) {
  peptide_data_prep <- peptide_data %>%
    dplyr::distinct({{ peptide_sequence }}, {{ start }}, {{ end }}, {{ uniprot_id }})

  if (missing(pdb_data)) {
    # if pdb_data is not provided by the user, it is fetched
    unis <- unique(dplyr::pull(peptide_data_prep, {{ uniprot_id }}))

    pdb_id_mapping <- fetch_uniprot(unis, columns = c("database(PDB)")) %>%
      tidyr::drop_na() %>%
      dplyr::mutate(pdb_ids = strsplit(.data$database_pdb, split = ";")) %>%
      tidyr::unnest(.data$pdb_ids)

    pdb_data <- fetch_pdb(pdb_ids = unique(pdb_id_mapping$pdb_ids))
  }

  result <- pdb_data %>%
    dplyr::distinct(.data$reference_database_accession, .data$pdb_ids, .data$auth_chain, .data$entity_beg_seq_id, .data$ref_beg_seq_id, .data$length, .data$pdb_sequence, .data$auth_to_entity_poly_seq_mapping, .data$database_chain) %>%
    dplyr::rename(length_pdb = .data$length) %>%
    dplyr::mutate({{ uniprot_id }} := .data$reference_database_accession) %>%
    dplyr::mutate(length_pdb_sequence = nchar(.data$pdb_sequence)) %>%
    dplyr::mutate(
      ref_end_seq_id = as.numeric(.data$ref_beg_seq_id) + as.numeric(.data$length_pdb) - 1,
      entity_end_seq_id = as.numeric(.data$entity_beg_seq_id) + as.numeric(.data$length_pdb) - 1
    ) %>%
    dplyr::select(.data$pdb_ids, .data$auth_chain, {{ uniprot_id }}, .data$entity_beg_seq_id, .data$entity_end_seq_id, .data$ref_beg_seq_id, .data$ref_end_seq_id, .data$pdb_sequence, .data$length_pdb_sequence, .data$auth_to_entity_poly_seq_mapping, .data$database_chain) %>%
    dplyr::right_join(peptide_data_prep, by = c(rlang::as_name(rlang::enquo(uniprot_id)))) %>%
    dplyr::mutate(peptide_in_pdb = ({{ start }} >= .data$ref_beg_seq_id & {{ start }} <= .data$ref_end_seq_id) | ({{ end }} >= .data$ref_beg_seq_id & {{ end }} <= .data$ref_end_seq_id) | ({{ start }} < .data$ref_beg_seq_id & {{ end }} > .data$ref_end_seq_id)) %>%
    dplyr::group_by(.data$pdb_ids, {{ uniprot_id }}, .data$auth_chain) %>%
    dplyr::mutate(n_peptides = dplyr::n_distinct({{ peptide_sequence }})) %>%
    dplyr::filter(.data$peptide_in_pdb) %>%
    dplyr::mutate(n_peptides_in_pdb = dplyr::n_distinct({{ peptide_sequence }})) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      peptide_start_pdb_sequence = .data$entity_beg_seq_id - .data$ref_beg_seq_id + {{ start }},
      peptide_end_pdb_sequence = .data$entity_beg_seq_id - .data$ref_beg_seq_id + {{ end }}
    ) %>%
    dplyr::mutate(fit_type = ifelse((.data$peptide_start_pdb_sequence < .data$entity_beg_seq_id | .data$peptide_end_pdb_sequence > .data$entity_end_seq_id), "partial", "fully")) %>%
    dplyr::mutate(
      peptide_end_pdb_sequence = ifelse(.data$peptide_end_pdb_sequence > .data$length_pdb_sequence, .data$length_pdb_sequence, .data$peptide_end_pdb_sequence),
      peptide_start_pdb_sequence = ifelse(.data$peptide_start_pdb_sequence < 1, 1, .data$peptide_start_pdb_sequence),
      peptide_end_pdb_sequence = ifelse(.data$peptide_end_pdb_sequence > .data$entity_end_seq_id, .data$entity_end_seq_id, .data$peptide_end_pdb_sequence),
      peptide_start_pdb_sequence = ifelse(.data$peptide_start_pdb_sequence < .data$entity_beg_seq_id, .data$entity_beg_seq_id, .data$peptide_start_pdb_sequence)
    ) %>%
    dplyr::mutate(peptide_sequence_in_pdb = stringr::str_sub(.data$pdb_sequence, start = .data$peptide_start_pdb_sequence, end = .data$peptide_end_pdb_sequence)) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      peptide_start_pdb = suppressWarnings(as.numeric(.data$auth_to_entity_poly_seq_mapping[[.data$peptide_start_pdb_sequence]])),
      peptide_end_pdb = suppressWarnings(as.numeric(.data$auth_to_entity_poly_seq_mapping[[.data$peptide_end_pdb_sequence]]))
    ) %>%
    dplyr::ungroup() %>%
    dplyr::select({{ uniprot_id }}, .data$pdb_ids, .data$auth_chain, .data$database_chain, {{ peptide_sequence }}, .data$peptide_sequence_in_pdb, .data$fit_type, {{ start }}, {{ end }}, .data$peptide_start_pdb_sequence, .data$peptide_end_pdb_sequence, .data$peptide_start_pdb, .data$peptide_end_pdb, .data$n_peptides, .data$n_peptides_in_pdb)

  # Retain also peptides in the data frame that were not found in any pdb structure

  output <- result %>%
    dplyr::right_join(peptide_data %>% dplyr::distinct({{ peptide_sequence }}, {{ start }}, {{ end }}, {{ uniprot_id }}), by = c(
      rlang::as_name(rlang::enquo(peptide_sequence)),
      rlang::as_name(rlang::enquo(start)),
      rlang::as_name(rlang::enquo(end)),
      rlang::as_name(rlang::enquo(uniprot_id))
    ))

  if (!missing(retain_columns)) {
    output <- peptide_data %>%
      dplyr::select(!!enquo(retain_columns), colnames(output)[!colnames(output) %in% c(
        "pdb_ids",
        "auth_chain",
        "database_chain",
        "peptide_sequence_in_pdb",
        "fit_type",
        "peptide_start_pdb_sequence",
        "peptide_end_pdb_sequence",
        "peptide_start_pdb",
        "peptide_end_pdb",
        "n_peptides",
        "n_peptides_in_pdb"
      )]) %>%
      dplyr::distinct() %>%
      dplyr::right_join(output, by = colnames(output)[!colnames(output) %in% c(
        "pdb_ids",
        "auth_chain",
        "database_chain",
        "peptide_sequence_in_pdb",
        "fit_type",
        "peptide_start_pdb_sequence",
        "peptide_end_pdb_sequence",
        "peptide_start_pdb",
        "peptide_end_pdb",
        "n_peptides",
        "n_peptides_in_pdb"
      )])
  }

  output
}
