#' Finds peptide positions in a PDB structure based on positional matching
#'
#' Finds peptide positions in a PDB structure. Often positions of peptides in UniProt and a PDB
#' structure are different due to different lengths of structures. This function maps a peptide
#' based on its UniProt positions onto a PDB structure. This method is superior to sequence
#' alignment of the peptide to the PDB structure sequence, since it can also match the peptide if
#' there are truncations or mismatches. This function also provides an easy way to check if a
#' peptide is present in a PDB structure.
#'
#' @param peptide_data a data frame containing at least the input columns to this function.
#' @param peptide a character column in the \code{peptide_data} data frame that contains the
#' sequence or any other unique identifier for the peptide that should be found.
#' @param start a numeric column in the \code{peptide_data} data frame that contains start positions
#' of peptides.
#' @param end a numeric column in the \code{peptide_data} data frame that contains end positions of
#' peptides.
#' @param uniprot_id a character column in the \code{peptide_data} data frame that contains UniProt
#' identifiers that correspond to the peptides.
#' @param pdb_data optional, a data frame containing data obtained with \code{fetch_pdb()}. If not
#' provided, information is fetched automatically. If this function should be run multiple times
#' it is faster to fetch the information once and provide it to the function. If provided, make
#' sure that the column names are identical to the ones that would be obtained by calling \code{fetch_pdb()}.
#' @param retain_columns a vector indicating if certain columns should be retained from the input
#' data frame. Default is not retaining additional columns \code{retain_columns = NULL}. Specific
#' columns can be retained by providing their names (not in quotations marks, just like other
#' column names, but in a vector).
#'
#' @return A data frame that contains peptide positions in the corresponding PDB structures. If a
#' peptide is not found in any structure or no structure is associated with the protein, the data
#' frame contains NAs values for the output columns. The data frame contains the following and
#' additional columns:
#'
#' * auth_asym_id: Chain identifier provided by the author of the structure in order to
#' match the identification used in the publication that describes the structure.
#' * label_asym_id: Chain identifier following the standardised convention for mmCIF files.
#' * peptide_seq_in_pdb: The sequence of the peptide mapped to the structure. If the
#' peptide only maps partially, then only the part of the sequence that maps on the structure is
#' returned.
#' * fit_type: The fit type is either "partial" or "fully" and it indicates if the complete
#' peptide or only part of it was found in the structure.
#' * start_adjusted: The adjusted start position of the peptide if the peptide was only partially
#' covered by the structure. If the peptide is fully covered it is just the provided start position.
#' * end_adjusted: The adjusted end position of the peptide if the peptide was only partially
#' covered by the structure. If the peptide is fully covered it is just the provided end position.
#' * label_seq_id_start: Contains the first residue position of the peptide in the structure
#' following the standardised convention for mmCIF files.
#' * label_seq_id_end: Contains the last residue position of the peptide in the structure
#' following the standardised convention for mmCIF files.
#' * auth_seq_id_start: Contains the first residue position of the peptide in the structure
#' based on the alternative residue identifier provided by the author of the structure in order
#' to match the identification used in the publication that describes the structure. This does
#' not need to be numeric and is therefore of type character.
#' * auth_seq_id_end: Contains the last residue position of the peptide in the structure
#' based on the alternative residue identifier provided by the author of the structure in order
#' to match the identification used in the publication that describes the structure. This does
#' not need to be numeric and is therefore of type character.
#' * auth_seq_id: Contains all positions (separated by ";") of the peptide in the structure
#' based on the alternative residue identifier provided by the author of the structure in order
#' to match the identification used in the publication that describes the structure. This does
#' not need to be numeric and is therefore of type character.
#' * n_peptides: The number of peptides from one protein that were searched for within the
#' current structure.
#' * n_peptides_in_structure: The number of peptides from one protein that were found within
#' the current structure.
#' * percentage_covered_peptides: The percentage of all provided peptides that were at least
#' partially found in the structure.
#'
#' @import dplyr
#' @import tidyr
#' @importFrom stringr str_sub str_split
#' @importFrom magrittr %>%
#' @importFrom rlang as_name enquo .data
#' @export
#'
#' @examples
#' \donttest{
#' # Create example data
#' peptide_data <- data.frame(
#'   uniprot_id = c("P0A8T7", "P0A8T7", "P60906"),
#'   peptide_sequence = c(
#'     "SGIVSFGKETKGKRRLVITPVDGSDPYEEMIPKWRQLNV",
#'     "NVFEGERVER",
#'     "AIGEVTDVVEKE"
#'   ),
#'   start = c(1160, 1197, 55),
#'   end = c(1198, 1206, 66)
#' )
#'
#' # Find peptides in protein structure
#' peptide_in_structure <- find_peptide_in_structure(
#'   peptide_data = peptide_data,
#'   peptide = peptide_sequence,
#'   start = start,
#'   end = end,
#'   uniprot_id = uniprot_id
#' )
#'
#' head(peptide_in_structure, n = 10)
#' }
find_peptide_in_structure <- function(peptide_data,
                                      peptide,
                                      start,
                                      end,
                                      uniprot_id,
                                      pdb_data = NULL,
                                      retain_columns = NULL) {
  peptide_data_prep <- peptide_data %>%
    dplyr::distinct({{ peptide }}, {{ start }}, {{ end }}, {{ uniprot_id }})

  if (missing(pdb_data)) {
    # if pdb_data is not provided by the user, it is fetched
    unis <- unique(dplyr::pull(peptide_data_prep, {{ uniprot_id }}))

    uniprot_info <- fetch_uniprot(unis, columns = c("xref_pdb"))

    # Make sure to only execute the code below if there are
    # PDB structures to be extracted and fetched.
    if (!all(is.na(uniprot_info$xref_pdb))) {
      pdb_id_mapping <- uniprot_info %>%
        tidyr::drop_na() %>%
        dplyr::mutate(pdb_ids = strsplit(.data$xref_pdb, split = ";")) %>%
        tidyr::unnest("pdb_ids")

      pdb_data <- fetch_pdb(pdb_ids = unique(pdb_id_mapping$pdb_ids))
    }
  }

  # Rescue the case that none of the provided proteins has a structure
  # or no matching structure was provided in pdb_data. Then the output
  # data frame is returned that contains NAs in all fields usually
  # added by this function.

  if (ifelse(missing(pdb_data),
    all(is.na(uniprot_info$xref_pdb)),
    !any(unique(dplyr::pull(peptide_data_prep, {{ uniprot_id }})) %in%
      pdb_data$reference_database_accession)
  )) {
    result <- peptide_data_prep %>%
      dplyr::mutate(
        pdb_ids = NA,
        auth_asym_id = NA,
        label_asym_id = NA,
        peptide_seq_in_pdb = NA,
        start_adjusted = NA,
        end_adjusted = NA,
        fit_type = NA,
        label_seq_id_start = {{ start }},
        label_seq_id_end = {{ end }},
        auth_seq_id_start = as.character({{ start }}),
        auth_seq_id_end = as.character({{ end }}),
        n_peptides = NA,
        n_peptides_in_structure = NA
      ) %>%
      dplyr::group_by({{ peptide }}, .data$auth_seq_id_start, .data$auth_seq_id_end) %>%
      dplyr::mutate(auth_seq_id = ifelse(!is.na(.data$label_seq_id_start) &
        !is.na(.data$label_seq_id_end),
      list(as.character(seq(.data$label_seq_id_start, .data$label_seq_id_end))),
      list(NA)
      )) %>%
      dplyr::mutate(auth_seq_id = paste0(.data$auth_seq_id[[1]], collapse = ";")) %>%
      dplyr::ungroup()
  } else {
    result <- pdb_data %>%
      dplyr::distinct(
        .data$reference_database_accession,
        .data$pdb_ids,
        .data$auth_asym_id,
        .data$entity_beg_seq_id,
        .data$ref_beg_seq_id,
        .data$length,
        .data$pdb_sequence,
        .data$auth_seq_id,
        .data$label_asym_id
      ) %>%
      dplyr::rename(length_pdb = "length") %>%
      dplyr::mutate({{ uniprot_id }} := .data$reference_database_accession) %>%
      dplyr::mutate(length_pdb_sequence = nchar(.data$pdb_sequence)) %>%
      dplyr::mutate(
        ref_end_seq_id = as.numeric(.data$ref_beg_seq_id) + as.numeric(.data$length_pdb) - 1,
        entity_end_seq_id = as.numeric(.data$entity_beg_seq_id) + as.numeric(.data$length_pdb) - 1
      ) %>%
      dplyr::select(
        "pdb_ids",
        "auth_asym_id",
        {{ uniprot_id }},
        "entity_beg_seq_id",
        "entity_end_seq_id",
        "ref_beg_seq_id",
        "ref_end_seq_id",
        "pdb_sequence",
        "length_pdb_sequence",
        "auth_seq_id",
        "label_asym_id"
      ) %>%
      dplyr::right_join(peptide_data_prep, by = c(rlang::as_name(rlang::enquo(uniprot_id))), relationship = "many-to-many") %>%
      dplyr::mutate(peptide_in_pdb = ({{ start }} >= .data$ref_beg_seq_id &
        {{ start }} <= .data$ref_end_seq_id) |
        ({{ end }} >= .data$ref_beg_seq_id &
          {{ end }} <= .data$ref_end_seq_id) |
        ({{ start }} < .data$ref_beg_seq_id &
          {{ end }} > .data$ref_end_seq_id)) %>%
      dplyr::group_by(.data$pdb_ids, .data$auth_asym_id) %>%
      dplyr::mutate(n_peptides = dplyr::n_distinct({{ peptide }})) %>%
      tidyr::drop_na("pdb_ids") %>%
      dplyr::mutate(n_peptides_in_structure = sum(.data$peptide_in_pdb)) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(
        label_seq_id_start = .data$entity_beg_seq_id - .data$ref_beg_seq_id + {{ start }},
        label_seq_id_end = .data$entity_beg_seq_id - .data$ref_beg_seq_id + {{ end }}
      ) %>%
      dplyr::mutate(fit_type = ifelse((.data$label_seq_id_start < .data$entity_beg_seq_id |
        .data$label_seq_id_end > .data$entity_end_seq_id),
      "partial",
      "fully"
      )) %>%
      # find the actual start and end if the coverage is only partial
      dplyr::mutate(
        start_adjusted = ifelse(.data$label_seq_id_start < .data$entity_beg_seq_id,
          {{ start }} - (.data$label_seq_id_start - .data$entity_beg_seq_id),
          {{ start }}
        ),
        end_adjusted = ifelse((.data$label_seq_id_end > .data$entity_end_seq_id),
          {{ end }} - (.data$label_seq_id_end - .data$entity_end_seq_id),
          {{ end }}
        )
      ) %>%
      dplyr::mutate(
        label_seq_id_end = ifelse(.data$label_seq_id_end > .data$length_pdb_sequence,
          .data$length_pdb_sequence,
          .data$label_seq_id_end
        ),
        label_seq_id_start = ifelse(.data$label_seq_id_start < 1,
          1,
          .data$label_seq_id_start
        ),
        label_seq_id_end = ifelse(.data$label_seq_id_end > .data$entity_end_seq_id,
          .data$entity_end_seq_id,
          .data$label_seq_id_end
        ),
        label_seq_id_start = ifelse(.data$label_seq_id_start < .data$entity_beg_seq_id,
          .data$entity_beg_seq_id,
          .data$label_seq_id_start
        )
      ) %>%
      dplyr::mutate(peptide_seq_in_pdb = stringr::str_sub(
        .data$pdb_sequence,
        start = .data$label_seq_id_start,
        end = .data$label_seq_id_end
      )) %>%
      dplyr::mutate(
        label_seq_id_start = ifelse(.data$peptide_in_pdb, .data$label_seq_id_start, Inf),
        label_seq_id_end = ifelse(.data$peptide_in_pdb, .data$label_seq_id_end, Inf),
        fit_type = ifelse(.data$peptide_in_pdb, .data$fit_type, NA),
        peptide_seq_in_pdb = ifelse(.data$peptide_in_pdb, .data$peptide_seq_in_pdb, NA)
      ) %>%
      dplyr::mutate(auth_seq_id_vector = stringr::str_split(.data$auth_seq_id, pattern = ";")) %>%
      dplyr::rowwise() %>%
      dplyr::mutate(
        auth_seq_id_start = suppressWarnings(.data$auth_seq_id_vector[.data$label_seq_id_start]),
        auth_seq_id_end = suppressWarnings(.data$auth_seq_id_vector[.data$label_seq_id_end])
      ) %>%
      dplyr::mutate(auth_seq_id = ifelse(.data$label_seq_id_start != Inf &
        .data$label_seq_id_end != Inf,
      paste0(.data$auth_seq_id_vector[seq(.data$label_seq_id_start, .data$label_seq_id_end)], collapse = ";"),
      NA
      )) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(
        label_seq_id_start = ifelse(.data$label_seq_id_start != Inf, .data$label_seq_id_start, NA),
        label_seq_id_end = ifelse(.data$label_seq_id_end != Inf, .data$label_seq_id_end, NA),
        start_adjusted = ifelse(.data$label_seq_id_start != Inf, .data$start_adjusted, NA),
        end_adjusted = ifelse(.data$label_seq_id_end != Inf, .data$end_adjusted, NA)
      ) %>%
      # calculate percentage of detected peptides that are present in structure
      # can be over 100% if the peptide is found multiple times
      dplyr::mutate(percentage_covered_peptides = .data$n_peptides_in_structure / .data$n_peptides * 100) %>%
      dplyr::select(
        {{ uniprot_id }},
        "pdb_ids",
        "auth_asym_id",
        "label_asym_id",
        {{ peptide }},
        "peptide_seq_in_pdb",
        "fit_type",
        {{ start }},
        "start_adjusted",
        {{ end }},
        "end_adjusted",
        "label_seq_id_start",
        "label_seq_id_end",
        "auth_seq_id_start",
        "auth_seq_id_end",
        "auth_seq_id",
        "n_peptides",
        "n_peptides_in_structure",
        "percentage_covered_peptides"
      )
  }
  # Retain also peptides in the data frame that were not found in any pdb structure or of
  # which no peptide was found in the available pdb structures.

  # All proteins with structures and that map at least one of their peptides
  good_result <- result %>%
    dplyr::filter(.data$n_peptides_in_structure > 0)

  # All proteins with no structures or no mapped peptides
  missing_result <- peptide_data_prep %>%
    filter(!{{ uniprot_id }} %in% unique(dplyr::pull(good_result, {{ uniprot_id }})))

  if (nrow(missing_result) > 0) {
    # only run the code below if the missing_result data frame contains data, otherwise
    # the function would return an error.
    missing_result <- missing_result %>%
      # Make structure positions equal to start and end positions if no structure is
      # available so this can be used for predictions
      dplyr::mutate(
        label_seq_id_start = {{ start }},
        label_seq_id_end = {{ end }},
        auth_seq_id_start = as.character({{ start }}),
        auth_seq_id_end = as.character({{ end }})
      ) %>%
      dplyr::group_by({{ peptide }}, .data$auth_seq_id_start, .data$auth_seq_id_end) %>%
      dplyr::mutate(auth_seq_id = ifelse(!is.na(.data$label_seq_id_start) &
        !is.na(.data$label_seq_id_end),
      list(as.character(seq(.data$label_seq_id_start, .data$label_seq_id_end))),
      list(NA)
      )) %>%
      dplyr::mutate(auth_seq_id = paste0(.data$auth_seq_id[[1]], collapse = ";")) %>%
      dplyr::ungroup()
  }
  # Add back info that does not have any structures matching the peptides
  output <- result %>%
    dplyr::bind_rows(missing_result) %>%
    dplyr::distinct()

  if (!missing(retain_columns)) {
    output <- peptide_data %>%
      dplyr::select(!!enquo(retain_columns), colnames(output)[!colnames(output) %in% c(
        "pdb_ids",
        "auth_asym_id",
        "label_asym_id",
        "peptide_seq_in_pdb",
        "fit_type",
        "label_seq_id_start",
        "label_seq_id_end",
        "auth_seq_id_start",
        "auth_seq_id_end",
        "auth_seq_id",
        "n_peptides",
        "n_peptides_in_structure",
        "start_adjusted",
        "end_adjusted",
        "percentage_covered_peptides"
      )]) %>%
      dplyr::distinct() %>%
      dplyr::right_join(output, by = colnames(output)[!colnames(output) %in% c(
        "pdb_ids",
        "auth_asym_id",
        "label_asym_id",
        "peptide_seq_in_pdb",
        "fit_type",
        "label_seq_id_start",
        "label_seq_id_end",
        "auth_seq_id_start",
        "auth_seq_id_end",
        "auth_seq_id",
        "n_peptides",
        "n_peptides_in_structure",
        "start_adjusted",
        "end_adjusted",
        "percentage_covered_peptides"
      )])
  }

  output
}
