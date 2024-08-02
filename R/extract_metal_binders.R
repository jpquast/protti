#' Extract metal-binding protein information from UniProt
#'
#' Information of metal binding proteins is extracted from UniProt data retrieved with
#' \code{fetch_uniprot} as well as QuickGO data retrieved with \code{fetch_quickgo}.
#'
#' @param data_uniprot a data frame containing at least the `ft_binding`, `cc_cofactor`,
#' `cc_catalytic_activity` and `keyword` columns.
#' @param data_quickgo a data frame containing molecular function gene ontology information for at
#' least the proteins of interest. This data should be obtained by calling \code{fetch_quickgo()}.
#' @param data_chebi optional, a data frame that can be manually obtained with \code{fetch_chebi(stars = c(2, 3))}.
#' It should contain 2 and 3 star entries. If not provided it will be fetched within the function. If the
#' function is run many times it is recommended to provide the data frame to save time.
#' @param data_chebi_relation optional, a data frame that can be manually obtained with
#' \code{fetch_chebi(relation = TRUE)}. If not provided it will be fetched within the function.
#' If the function is run many times it is recommended to provide the data frame to save time.
#' @param data_eco optional, a data frame that contains evidence and conclusion ontology data that can be
#' obtained by calling \code{fetch_eco()}. If not provided it will be fetched within the function.
#' If the function is run many times it is recommended to provide the data frame to save time.
#' @param data_eco_relation optional, a data frame that contains relational evidence and conclusion
#' ontology data that can be obtained by calling \code{fetch_eco(return_relation = TRUE)}. If not provided it
#' will be fetched within the function. If the function is run many times it is recommended to provide
#' the data frame to save time.
#' @param show_progress a logical value that specifies if progress will be shown (default is TRUE).
#'
#' @return A data frame containing information on protein metal binding state. It contains the
#' following columns:
#'
#' * \code{accession}: UniProt protein identifier.
#' * \code{most_specific_id}: ChEBI ID that is most specific for the position after combining information from all sources.
#' Can be multiple IDs separated by "," if a position appears multiple times due to multiple fitting IDs.
#' * \code{most_specific_id_name}: The name of the ID in the \code{most_specific_id} column. This information is based on
#' ChEBI.
#' * \code{ligand_identifier}: A ligand identifier that is unique per ligand per protein. It consists of the ligand ID and
#' ligand name. The ligand ID counts the number of ligands of the same type per protein.
#' * \code{ligand_position}: The amino acid position of the residue interacting with the ligand.
#' * \code{binding_mode}: Contains information about the way the amino acid residue interacts with the ligand. If it is
#' "covalent" then the residue is not in contact with the metal directly but only the cofactor that binds the metal.
#' * \code{metal_function}: Contains information about the function of the metal. E.g. "catalytic".
#' * \code{metal_id_part}: Contains a ChEBI ID that identifiers the metal part of the ligand. This is always the metal atom.
#' * \code{metal_id_part_name}: The name of the ID in the \code{metal_id_part} column. This information is based on
#' ChEBI.
#' * \code{note}: Contains notes associated with information based on cofactors.
#' * \code{chebi_id}: Contains the original ChEBI IDs the information is based on.
#' * \code{source}: Contains the sources of the information. This can consist of "binding", "cofactor", "catalytic_activity"
#' and "go_term".
#' * \code{eco}: If there is evidence the annotation is based on it is annotated with an ECO ID, which is split by source.
#' * \code{eco_type}: The ECO identifier can fall into the "manual_assertion" group for manually curated annotations or the
#' "automatic_assertion" group for automatically generated annotations. If there is no evidence it is annotated as
#' "automatic_assertion". The information is split by source.
#' * \code{evidence_source}: The original sources (e.g. literature, PDB) of evidence annotations split by source.
#' * \code{reaction}: Contains information about the chemical reaction catalysed by the protein that involves the metal.
#' Can contain the EC ID, Rhea ID, direction specific Rhea ID, direction of the reaction and evidence for the direction.
#' * \code{go_term}: Contains gene ontology terms if there are any metal related ones associated with the annotation.
#' * \code{go_name}: Contains gene ontology names if there are any metal related ones associated with the annotation.
#' * \code{assigned_by}: Contains information about the source of the gene ontology term assignment.
#' * \code{database}: Contains information about the source of the ChEBI annotation associated with gene ontology terms.
#' * `keyword`: Contains keywords if they were annotated in UniProt.
#'
#' For each protein identifier the data frame contains information on the bound ligand as well as on its position if it is known.
#' Since information about metal ligands can come from multiple sources, additional information (e.g. evidence) is nested in the returned
#' data frame. In order to unnest the relevant information the following steps have to be taken: It is
#' possible that there are multiple IDs in the "most_specific_id" column. This means that one position cannot be uniquely
#' attributed to one specific ligand even with the same ligand_identifier. Apart from the "most_specific_id" column, in
#' which those instances are separated by ",", in other columns the relevant information is separated by "||". Then
#' information should be split based on the source (not the \code{source} column, that one can be removed from the data
#' frame). There are certain columns associated with specific sources (e.g. \code{go_term} is associated
#' with the \code{"go_term"} source). Values of columns not relevant for a certain source should be replaced with \code{NA}.
#' Since a \code{most_specific_id} can have multiple \code{chebi_id}s associated with it we need to unnest the \code{chebi_id}
#' column and associated columns in which information is separated by "|". Afterwards evidence and additional information can be
#' unnested by first splitting data for ";;" and then for ";".
#' @import dplyr
#' @import tidyr
#' @import purrr
#' @import stringr
#' @importFrom stats na.omit setNames
#' @importFrom rlang .data as_name enquo
#' @importFrom magrittr %>%
#' @export
#'
#' @examples
#' \donttest{
#' # Create example data
#'
#' uniprot_ids <- c("P00393", "P06129", "A0A0C5Q309", "A0A0C9VD04")
#'
#' ## UniProt data
#' data_uniprot <- fetch_uniprot(
#'   uniprot_ids = uniprot_ids,
#'   columns = c(
#'     "ft_binding",
#'     "cc_cofactor",
#'     "cc_catalytic_activity"
#'   )
#' )
#'
#' ## QuickGO data
#' data_quickgo <- fetch_quickgo(
#'   id_annotations = uniprot_ids,
#'   ontology_annotations = "molecular_function"
#' )
#'
#' ## ChEBI data (2 and 3 star entries)
#' data_chebi <- fetch_chebi(stars = c(2, 3))
#' data_chebi_relation <- fetch_chebi(relation = TRUE)
#'
#' ## ECO data
#' eco <- fetch_eco()
#' eco_relation <- fetch_eco(return_relation = TRUE)
#'
#' # Extract metal binding information
#' metal_info <- extract_metal_binders(
#'   data_uniprot = data_uniprot,
#'   data_quickgo = data_quickgo,
#'   data_chebi = data_chebi,
#'   data_chebi_relation = data_chebi_relation,
#'   data_eco = eco,
#'   data_eco_relation = eco_relation
#' )
#'
#' metal_info
#' }
extract_metal_binders <- function(data_uniprot,
                                  data_quickgo,
                                  data_chebi = NULL,
                                  data_chebi_relation = NULL,
                                  data_eco = NULL,
                                  data_eco_relation = NULL,
                                  show_progress = TRUE) {
  metal_list <- protti::metal_list
  metal_chebi_uniprot <- protti::metal_chebi_uniprot
  metal_go_slim_subset <- protti::metal_go_slim_subset
  # Check if required R packages are installed
  if (!requireNamespace("igraph", quietly = TRUE)) {
    message("Package \"igraph\" is needed for this function to work. Please install it.", call. = FALSE)
    return(invisible(NULL))
  }
  if (!requireNamespace("stringi", quietly = TRUE)) {
    message("Package \"stringi\" is needed for this function to work. Please install it.", call. = FALSE)
    return(invisible(NULL))
  }

  # Check if data was provided, if not the connection might have failed
  if ((is.null(data_uniprot) & !missing(data_uniprot)) |
    (is.null(data_quickgo) & !missing(data_quickgo)) |
    (is.null(data_chebi) & !missing(data_chebi)) |
    (is.null(data_chebi_relation) & !missing(data_chebi_relation)) |
    (is.null(data_eco) & !missing(data_eco)) |
    (is.null(data_eco_relation) & !missing(data_eco_relation))) {
    message(strwrap('The provided "data_uniprot", "data_quickgo", "data_chebi", "data_eco" or "data_eco_relation" is NULL.
            Please check that there was no connection problem when you retrieved this data and try again.', prefix = "\n", initial = ""))
    return(invisible(NULL))
  }

  # Check if provided data has the right format
  # data_uniprot
  if (!("ft_binding" %in% colnames(data_uniprot) &
    "cc_cofactor" %in% colnames(data_uniprot) &
    "keyword" %in% colnames(data_uniprot) &
    "cc_catalytic_activity" %in% colnames(data_uniprot))) {
    stop('Please include at least the columns "ft_binding", "cc_cofactor", "cc_catalytic_activity" and "keyword" in "data_uniprot"!')
  }
  # data_quickgo
  if (!("gene_product_db" %in% colnames(data_quickgo) &
    "gene_product_id" %in% colnames(data_quickgo) &
    "go_name" %in% colnames(data_quickgo) &
    "go_term" %in% colnames(data_quickgo) &
    "go_aspect" %in% colnames(data_quickgo) &
    "eco_id" %in% colnames(data_quickgo) &
    "go_evidence_code" %in% colnames(data_quickgo) &
    "reference" %in% colnames(data_quickgo) &
    "with_from" %in% colnames(data_quickgo))) {
    stop('Please include at least the columns "gene_product_db", "gene_product_id", "go_name", "go_term", "go_aspect", "go_evidence_code", "reference" and "with_from" in "data_quickgo"!')
  }
  # Download ChEBI database if not provided
  if (missing(data_chebi)) {
    if (show_progress == TRUE) {
      message("Downloading ChEBI information ... ", appendLF = FALSE)
      start_time <- Sys.time()
    }

    data_chebi <- fetch_chebi(stars = c(2, 3))
    # If database could not be retrieved let the function return NULL.
    if (is.null(data_chebi)) {
      return(invisible(NULL))
    }

    if (show_progress == TRUE) {
      message("DONE ", paste0("(", round(as.numeric(difftime(Sys.time(), start_time, units = "secs")), digits = 2), "s)"))
    }
  }
  # Download chebi relation dataset if not provided
  if (missing(data_chebi_relation)) {
    if (show_progress == TRUE) {
      message("Downloading ChEBI relational information ... ", appendLF = FALSE)
      start_time <- Sys.time()
    }

    data_chebi_relation <- fetch_chebi(relation = TRUE)

    # If database could not be retrieved let the function return NULL.
    if (is.null(data_chebi_relation)) {
      return(invisible(NULL))
    }

    if (show_progress == TRUE) {
      message("DONE ", paste0("(", round(as.numeric(difftime(Sys.time(), start_time, units = "secs")), digits = 2), "s)"))
    }
  }

  # Retrieve ECO annotation information if not provided
  if (missing(data_eco)) {
    if (show_progress == TRUE) {
      message("Downloading ECO annotation information from QuickGO ... ", appendLF = FALSE)
      start_time <- Sys.time()
    }

    eco_data <- fetch_eco(show_progress = FALSE)

    # If database could not be retrieved let the function return NULL.
    if (is.null(eco_data)) {
      return(invisible(NULL))
    }

    if (show_progress == TRUE) {
      message("DONE ", paste0("(", round(as.numeric(difftime(Sys.time(), start_time, units = "secs")), digits = 2), "s)"))
    }
  }

  # Retrieve ECO relation information if not provided
  if (missing(data_eco_relation)) {
    if (show_progress == TRUE) {
      message("Downloading ECO relation information from QuickGO ... ", appendLF = FALSE)
      start_time <- Sys.time()
    }

    data_eco_relation <- fetch_eco(return_relation = TRUE, show_progress = FALSE)

    if (show_progress == TRUE) {
      message("DONE ", paste0("(", round(as.numeric(difftime(Sys.time(), start_time, units = "secs")), digits = 2), "s)"))
    }
  }

  if (show_progress == TRUE) {
    message("Preparing annotation data frames ... ", appendLF = FALSE)
    start_time <- Sys.time()
  }

  # Prepare ChEBI data frame for annotation and filtering:
  # We create a data frame that contains all metal related entries from ChEBI
  # These entries are identified based on the formula that contains a metal
  # protti provides a data frame that contains all UniProt related metal ChEBI IDs
  # As not all entries contain a formula, the protti data frame (metal_chebi_uniprot)
  # contains manual annotations for these entries.
  # This function uses the below created list for annotation. The list contains "metal-type" information that was
  # directly extracted from the formula and the non-formula containing metal entries from metal_chebi_uniprot.
  # Since metal_chebi_uniprot also contains formula containing metal entries we can check if there are new ChEBI IDs
  # in the result, which would indicate that things have changed since metal_chebi_uniprot was created and that
  # there might also be novel non-formula entries.
  data_chebi_filtered <- data_chebi %>%
    dplyr::filter(.data$type_name == "STANDARD") %>%
    dplyr::mutate(chebi_id = as.character(.data$id)) %>%
    dplyr::distinct(.data$chebi_id, .data$definition, .data$name, .data$formula) %>%
    dplyr::filter(stringr::str_detect(.data$formula, pattern = paste0(paste0(metal_list$symbol, collapse = "(?![:lower:])|"), "(?![:lower:])")) |
      .data$chebi_id %in% metal_chebi_uniprot$id) %>% # This recreates a version of the data frame provided by protti that contains all metal containing entries from ChEBI
    dplyr::mutate(extract_formula = stringr::str_extract_all(.data$formula, pattern = paste0(paste0(metal_list$symbol, collapse = "(?![:lower:])|"), "(?![:lower:])"))) %>%
    tidyr::unnest("extract_formula") %>%
    dplyr::mutate(metal_atom_id = ifelse(is.na(.data$extract_formula),
      stats::setNames(metal_chebi_uniprot$metal_atom_id, as.character(metal_chebi_uniprot$id))[.data$chebi_id],
      stats::setNames(metal_list$chebi_id, metal_list$symbol)[.data$extract_formula]
    )) %>%
    dplyr::select(-"extract_formula") %>%
    dplyr::group_by(.data$chebi_id) %>%
    dplyr::mutate(metal_atom_id = paste0(.data$metal_atom_id, collapse = ",")) %>%
    dplyr::distinct()

  chebi_names <- data_chebi %>%
    dplyr::filter(.data$type_name == "STANDARD") %>%
    dplyr::mutate(id = as.character(.data$id)) %>%
    dplyr::distinct(.data$id, .data$name)

  # Create two vectors that contain all IDs related to either manual or automatic assertion
  manual_eco <- data_eco_relation %>%
    find_all_subs(
      ids = c("ECO:0000352"),
      main_id = .data$main_id,
      type = "relation",
      accepted_types = "all"
    ) %>%
    unlist()

  automatic_eco <- data_eco_relation %>%
    find_all_subs(
      ids = c("ECO:0000501"),
      main_id = .data$main_id,
      type = "relation",
      accepted_types = "all"
    ) %>%
    unlist()

  if (show_progress == TRUE) {
    message("DONE ", paste0("(", round(as.numeric(difftime(Sys.time(), start_time, units = "secs")), digits = 2), "s)"))
  }

  # Extract ft_binding information from UniProt
  if (show_progress == TRUE) {
    message("Extract ft_binding information from UniProt ... ", appendLF = FALSE)
    start_time <- Sys.time()
  }

  b_uniprot <- data_uniprot %>%
    dplyr::distinct(.data$accession, .data$ft_binding) %>%
    tidyr::drop_na("ft_binding") %>%
    # Extract each position
    dplyr::mutate(ft_binding = stringr::str_extract_all(
      .data$ft_binding,
      pattern = "BINDING.+?(?=BINDING)|BINDING.+$"
    )) %>%
    tidyr::unnest("ft_binding") %>%
    dplyr::mutate(chebi_id = stringr::str_extract(.data$ft_binding, pattern = '(?<=/ligand_id=\\"ChEBI:CHEBI:)[^\\";]+(?=[\\";])')) %>%
    # Filter with the previously generated data_chebi_filtered to only keep metal ChEBI IDs
    dplyr::filter(.data$chebi_id %in% data_chebi_filtered$chebi_id) %>%
    dplyr::mutate(ligand_name = stringr::str_extract(.data$ft_binding, pattern = '(?<=/ligand=\\")[^\\";]+(?=[\\";])')) %>%
    dplyr::mutate(ligand_identifier = stringr::str_extract(
      .data$ft_binding,
      pattern = '(?<=/ligand_label=\\")[^\\";]+(?=[\\";])'
    )) %>%
    dplyr::mutate(ligand_position = stringr::str_extract(.data$ft_binding, pattern = "(?<=BINDING )[^;]+(?=;)")) %>%
    dplyr::mutate(
      isoform = stringr::str_extract(.data$ligand_position, pattern = "[^:]+(?=:)"),
      accession = ifelse(!is.na(.data$isoform),
        .data$isoform,
        .data$accession
      ),
      ligand_position = stringr::str_remove(.data$ligand_position, pattern = "[^:]+:")
    ) %>%
    dplyr::mutate(metal_id_part = stringr::str_extract(.data$ft_binding, pattern = '(?<=/ligand_part_id=\\"ChEBI:CHEBI:)[^\\";]+(?=[\\";])')) %>%
    dplyr::mutate(binding_mode = stringr::str_extract(.data$ft_binding, pattern = '(?<=/note=\\")[^\\";]+(?=[\\";])')) %>%
    dplyr::mutate(evidence = stringr::str_extract(.data$ft_binding, pattern = '(?<=/evidence=\\")[^\\";]+(?=[\\";])')) %>%
    dplyr::mutate(metal_function = stringr::str_extract(.data$ft_binding, pattern = '(?<=/ligand_note=\\")[^\\";]+(?=[\\";])')) %>%
    # Remove rows with the same position, metal etc., NA ligand_identifier and NA evidence
    dplyr::group_by(.data$accession, .data$chebi_id, .data$ligand_position) %>%
    dplyr::filter(!(is.na(.data$ligand_identifier) &
      dplyr::n() > 1 &
      is.na(.data$evidence) &
      !all(is.na(.data$evidence)))) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(evidence_split = stringr::str_split(.data$evidence, pattern = ", ")) %>%
    tidyr::unnest("evidence_split") %>%
    tidyr::separate(.data$evidence_split, into = c("eco", "evidence_source"), sep = "\\|", fill = "right") %>%
    dplyr::distinct() %>%
    dplyr::mutate(eco_type = dplyr::case_when(
      .data$eco %in% manual_eco ~ "manual_assertion",
      .data$eco %in% automatic_eco ~ "automatic_assertion"
    )) %>%
    dplyr::mutate(eco_type = ifelse(is.na(.data$eco_type), "automatic_assertion", .data$eco_type)) %>%
    # If identifier is missing use "1"
    dplyr::mutate(ligand_identifier = ifelse(is.na(.data$ligand_identifier),
      paste0("1(", .data$ligand_name, ")"),
      paste0(.data$ligand_identifier, "(", .data$ligand_name, ")")
    )) %>%
    dplyr::select(-c("ft_binding", "evidence", "isoform", "ligand_name")) %>%
    dplyr::distinct() %>%
    # Extract metal positions
    dplyr::mutate(ligand_position = stringr::str_split(.data$ligand_position, pattern = "\\.\\.")) %>%
    dplyr::mutate(ligand_position = purrr::map(
      .x = .data$ligand_position,
      .f = ~ if (length(.x) > 1) {
        `:`(as.numeric(.x[1]), as.numeric(.x[2]))
      } else {
        as.numeric(.x)
      }
    )) %>%
    tidyr::unnest("ligand_position") %>%
    # Combine the binding_mode column to prevent duplicates
    # The reason for duplicates is a wrong annotation in UniProt (P00081).
    # There are also issues with additional IDs such as E3PRJ4, Q9DHD6, P00081.
    dplyr::group_by(.data$accession, .data$ligand_identifier, .data$ligand_position, .data$eco, .data$evidence_source, .data$chebi_id) %>%
    dplyr::mutate(
      binding_mode = paste0(stats::na.omit(unique(.data$binding_mode)), collapse = ","),
      metal_id_part = paste0(stats::na.omit(unique(.data$metal_id_part)), collapse = ","),
      metal_function = paste0(stats::na.omit(unique(.data$metal_function)), collapse = ",")
    ) %>%
    dplyr::ungroup() %>%
    dplyr::distinct() %>%
    dplyr::mutate(
      binding_mode = ifelse(.data$binding_mode == "", NA, .data$binding_mode),
      metal_id_part = ifelse(.data$metal_id_part == "", NA, .data$metal_id_part),
      metal_function = ifelse(.data$metal_function == "", NA, .data$metal_function)
    ) %>%
    # Combine the evidence_source column
    dplyr::group_by(.data$accession, .data$ligand_identifier, .data$ligand_position, .data$eco, .data$chebi_id) %>%
    dplyr::mutate(evidence_source = paste0(.data$evidence_source, collapse = ";")) %>%
    dplyr::ungroup() %>%
    dplyr::distinct() %>%
    # Combine the eco, eco_type and evidence_source column for each position
    dplyr::group_by(.data$accession, .data$ligand_identifier, .data$ligand_position, .data$chebi_id) %>%
    dplyr::mutate(
      evidence_source = paste0(.data$evidence_source, collapse = ";;"),
      eco = paste0(.data$eco, collapse = ";;"),
      eco_type = paste0(.data$eco_type, collapse = ";;")
    ) %>%
    dplyr::mutate(
      binding_mode = paste0(stats::na.omit(unique(.data$binding_mode)), collapse = ","),
      metal_id_part = paste0(stats::na.omit(unique(.data$metal_id_part)), collapse = ","),
      metal_function = paste0(stats::na.omit(unique(.data$metal_function)), collapse = ",")
    ) %>%
    dplyr::ungroup() %>%
    dplyr::distinct() %>%
    dplyr::mutate(
      binding_mode = ifelse(.data$binding_mode == "", NA, .data$binding_mode),
      metal_id_part = ifelse(.data$metal_id_part == "", NA, .data$metal_id_part),
      metal_function = ifelse(.data$metal_function == "", NA, .data$metal_function)
    ) %>%
    # Update metal_id_part position using the data provided by protti
    dplyr::left_join(dplyr::distinct(data_chebi_filtered, .data$chebi_id, .data$metal_atom_id), by = "chebi_id") %>%
    dplyr::group_by(.data$accession, .data$chebi_id, .data$ligand_identifier) %>%
    dplyr::mutate(
      metal_atom_id = ifelse(all(is.na(.data$metal_id_part)),
        .data$metal_atom_id,
        .data$metal_id_part
      ),
      # I don't know why this doesn't work in one go...
      metal_id_part = ifelse(is.na(.data$metal_atom_id),
        .data$metal_id_part,
        .data$metal_atom_id
      )
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(metal_id_part = ifelse(.data$binding_mode == "covalent" & !is.na(.data$binding_mode),
      NA,
      .data$metal_id_part
    )) %>%
    dplyr::select(-"metal_atom_id") %>%
    # Combine positions
    dplyr::group_by(.data$accession, .data$chebi_id, .data$ligand_identifier) %>%
    dplyr::mutate(
      ligand_position = paste0(.data$ligand_position, collapse = ";;;"),
      eco = paste0(.data$eco, collapse = ";;;"),
      eco_type = paste0(.data$eco_type, collapse = ";;;"),
      evidence_source = paste0(.data$evidence_source, collapse = ";;;"),
      binding_mode = paste0(.data$binding_mode, collapse = ";;;"),
      metal_function = paste0(.data$metal_function, collapse = ";;;"),
      metal_id_part = paste0(.data$metal_id_part, collapse = ";;;")
    ) %>%
    dplyr::ungroup() %>%
    dplyr::distinct() %>%
    dplyr::group_by(.data$accession, .data$chebi_id) %>%
    dplyr::mutate(
      ligand_identifier = paste0(.data$ligand_identifier, collapse = ";;;;"),
      ligand_position = paste0(.data$ligand_position, collapse = ";;;;"),
      eco = paste0(.data$eco, collapse = ";;;;"),
      eco_type = paste0(.data$eco_type, collapse = ";;;;"),
      evidence_source = paste0(.data$evidence_source, collapse = ";;;;"),
      binding_mode = paste0(.data$binding_mode, collapse = ";;;;"),
      metal_function = paste0(.data$metal_function, collapse = ";;;;"),
      metal_id_part = paste0(.data$metal_id_part, collapse = ";;;;")
    ) %>%
    dplyr::ungroup() %>%
    dplyr::distinct() %>%
    dplyr::rename(
      metal_id_part_binding = "metal_id_part",
      eco_binding = "eco",
      eco_type_binding = "eco_type",
      evidence_source_binding = "evidence_source"
    ) %>%
    dplyr::mutate(source = "binding") %>%
    # make sure that accession column is of type "chr" even if data frame is empty
    dplyr::mutate(accession = as.character(.data$accession))

  if (show_progress == TRUE) {
    message("DONE ", paste0("(", round(as.numeric(difftime(Sys.time(), start_time, units = "secs")), digits = 2), "s)"))
  }

  # Extract cc_cofactor information from UniProt
  if (show_progress == TRUE) {
    message("Extract cc_cofactor information from UniProt ... ", appendLF = FALSE)
    start_time <- Sys.time()
  }

  cofactor_uniprot <- data_uniprot %>%
    dplyr::distinct(.data$accession, .data$cc_cofactor) %>%
    tidyr::drop_na("cc_cofactor") %>%
    dplyr::mutate(cofactor_split = stringr::str_extract_all(
      .data$cc_cofactor,
      pattern = "(?<=COFACTOR:).+?(?=COFACTOR|$)"
    )) %>%
    tidyr::unnest("cofactor_split") %>%
    # Extract notes
    dplyr::mutate(note = stringr::str_extract(
      .data$cofactor_split,
      pattern = "(?<=Note\\=).+?(?=;)"
    )) %>%
    tidyr::unnest("note") %>%
    # Split names
    dplyr::mutate(name_split = stringr::str_extract_all(
      .data$cofactor_split,
      pattern = "(?<=Name\\=).+?(?=Name|Note|COFACTOR|$)"
    )) %>%
    tidyr::unnest("name_split") %>%
    # Extract ChEBI IDs from cc_cofactor
    dplyr::mutate(chebi_id = stringr::str_extract(
      .data$name_split,
      pattern = "(?<=CHEBI:)\\d+"
    )) %>%
    # Filter with data_chebi_filtered to only keep metal ChEBI IDs
    dplyr::filter(.data$chebi_id %in% data_chebi_filtered$chebi_id) %>%
    # Extract evidence from cc_cofactor
    dplyr::mutate(evidence = stringr::str_extract(.data$name_split, pattern = "(?<=Evidence=).+(?=;|$)")) %>%
    dplyr::mutate(evidence = stringr::str_remove_all(.data$evidence, pattern = ";|\\{|\\}")) %>%
    dplyr::mutate(evidence_split = stringr::str_split(.data$evidence, pattern = ", ")) %>%
    tidyr::unnest("evidence_split") %>%
    tidyr::separate(.data$evidence_split, into = c("eco", "evidence_source"), sep = "\\|", fill = "right") %>%
    dplyr::mutate(eco = stringr::str_trim(.data$eco)) %>%
    dplyr::distinct() %>%
    dplyr::mutate(eco_type = dplyr::case_when(
      .data$eco %in% manual_eco ~ "manual_assertion",
      .data$eco %in% automatic_eco ~ "automatic_assertion"
    )) %>%
    dplyr::mutate(eco_type = ifelse(is.na(.data$eco_type), "automatic_assertion", .data$eco_type)) %>%
    # dplyr::mutate(note_evidence = str_extract(.data$note,
    #                                           pattern = "(?<=\\{).+(?=\\})")) %>%
    dplyr::select(-c("cc_cofactor", "cofactor_split", "name_split", "evidence")) %>%
    # Add metal_id_part position using the data provided by protti
    dplyr::left_join(dplyr::distinct(data_chebi_filtered, .data$chebi_id, .data$metal_atom_id) %>%
      dplyr::rename(metal_id_part = "metal_atom_id"), by = "chebi_id") %>%
    # Now combine data to have one row per accession and chebi_id
    # First concatenate different notes for the same accession and chebi_id
    dplyr::group_by(.data$accession, .data$chebi_id) %>%
    dplyr::mutate(
      note = paste0(unique(.data$note), collapse = ","),
      note = stringr::str_replace_all(.data$note, pattern = "\\|", replacement = "/")
    ) %>%
    dplyr::ungroup() %>%
    dplyr::distinct() %>%
    # Combine evidence_sources
    dplyr::group_by(.data$accession, .data$eco, .data$chebi_id) %>%
    dplyr::mutate(evidence_source = paste0(.data$evidence_source, collapse = ";")) %>%
    dplyr::ungroup() %>%
    dplyr::distinct() %>%
    # Combine eco and eco_type
    dplyr::group_by(.data$accession, .data$chebi_id) %>%
    dplyr::mutate(
      evidence_source = paste0(.data$evidence_source, collapse = ";;"),
      eco = paste0(.data$eco, collapse = ";;"),
      eco_type = paste0(.data$eco_type, collapse = ";;")
    ) %>%
    dplyr::ungroup() %>%
    dplyr::distinct() %>%
    dplyr::mutate(source = "cofactor") %>%
    # make sure that accession column is of type "chr" even if data frame is empty
    dplyr::mutate(accession = as.character(.data$accession))

  if (show_progress == TRUE) {
    message("DONE ", paste0("(", round(as.numeric(difftime(Sys.time(), start_time, units = "secs")), digits = 2), "s)"))
  }

  # Extract cc_catalytic_activity information from UniProt
  if (show_progress == TRUE) {
    message("Extract cc_catalytic_activity information from UniProt ... ", appendLF = FALSE)
    start_time <- Sys.time()
  }

  catalytic_activity_uniprot <- data_uniprot %>%
    dplyr::distinct(.data$accession, .data$cc_catalytic_activity) %>%
    tidyr::drop_na("cc_catalytic_activity") %>%
    dplyr::mutate(catalytic_activity_split = stringr::str_extract_all(
      .data$cc_catalytic_activity,
      pattern = "(?<=CATALYTIC ACTIVITY:).+?(?=CATALYTIC ACTIVITY|$)"
    )) %>%
    tidyr::unnest("catalytic_activity_split") %>%
    dplyr::mutate(catalytic_activity_split = stringr::str_remove(
      stringr::str_trim(.data$catalytic_activity_split),
      pattern = "Reaction="
    )) %>%
    dplyr::mutate(chebi_id = stringr::str_extract_all(
      .data$catalytic_activity_split,
      pattern = "(?<=CHEBI:)\\d+"
    )) %>%
    tidyr::unnest("chebi_id") %>%
    # Filter with data_chebi_filtered to only keep metal ChEBI IDs
    dplyr::filter(.data$chebi_id %in% data_chebi_filtered$chebi_id) %>%
    dplyr::mutate(evidence = stringr::str_extract(.data$catalytic_activity_split, pattern = "(?<=Evidence=)[^;]+(?=;)")) %>%
    dplyr::mutate(evidence = stringr::str_remove_all(.data$evidence, pattern = ";|\\{|\\}")) %>%
    dplyr::mutate(evidence_split = stringr::str_split(.data$evidence, pattern = ", ")) %>%
    tidyr::unnest("evidence_split") %>%
    tidyr::separate(.data$evidence_split, into = c("eco", "evidence_source"), sep = "\\|", fill = "right") %>%
    dplyr::mutate(eco = stringr::str_trim(.data$eco)) %>%
    dplyr::distinct() %>%
    dplyr::mutate(eco_type = dplyr::case_when(
      .data$eco %in% manual_eco ~ "manual_assertion",
      .data$eco %in% automatic_eco ~ "automatic_assertion"
    )) %>%
    dplyr::mutate(eco_type = ifelse(is.na(.data$eco_type), "automatic_assertion", .data$eco_type)) %>%
    dplyr::mutate(reaction = stringr::str_extract(.data$catalytic_activity_split, pattern = "(?<=PhysiologicalDirection=).+(?=;$)")) %>%
    dplyr::mutate(reaction = stringr::str_split(.data$reaction, pattern = "PhysiologicalDirection=")) %>%
    tidyr::unnest("reaction") %>%
    dplyr::mutate(reaction = stringr::str_replace(.data$reaction, pattern = "Xref=Rhea:", replacement = "Direction")) %>%
    dplyr::mutate(rhea = stringr::str_extract_all(.data$catalytic_activity_split, pattern = "(?<=RHEA:)\\d+")) %>%
    dplyr::mutate(rhea = paste0(" RHEA:", purrr::map2_chr(
      .x = .data$rhea,
      .y = .data$reaction,
      .f = ~ {
        if (is.na(.y)) {
          paste0(.x, collapse = ",")
        } else {
          paste0(.x[!stringr::str_detect(.y, pattern = .x)], collapse = ",")
        }
      }
    ))) %>%
    dplyr::mutate(reaction = ifelse(!is.na(.data$reaction),
      paste0("Direction:", .data$reaction),
      NA
    )) %>%
    dplyr::mutate(ec = paste0("EC:", stringr::str_extract(.data$catalytic_activity_split, pattern = "(?<=EC=)[^;]+(?=;)"), ", ", .data$rhea)) %>%
    tidyr::unite(col = "reaction", c("ec", "reaction"), na.rm = TRUE, sep = ", ") %>%
    dplyr::mutate(
      reaction = stringr::str_replace_all(.data$reaction, pattern = "=", replacement = ":"),
      reaction = stringr::str_replace_all(.data$reaction, pattern = ";", replacement = ","),
      reaction = stringr::str_replace_all(.data$reaction, pattern = "\\|", replacement = "/")
    ) %>%
    dplyr::select(-c("cc_catalytic_activity", "catalytic_activity_split", "evidence", "rhea")) %>%
    # Add metal_id_part position using the data provided by protti
    dplyr::left_join(dplyr::distinct(data_chebi_filtered, .data$chebi_id, .data$metal_atom_id) %>%
      dplyr::rename(metal_id_part = "metal_atom_id"), by = "chebi_id") %>%
    # Now combine data to have one row per accession and chebi_id
    # First concatenate different evidences for the same accession, chebi_id, eco and reaction
    dplyr::group_by(.data$accession, .data$chebi_id, .data$eco, .data$reaction) %>%
    dplyr::mutate(evidence_source = paste0(.data$evidence_source, collapse = ",")) %>%
    dplyr::ungroup() %>%
    dplyr::distinct() %>%
    # Combine evidence_source and reaction
    dplyr::group_by(.data$accession, .data$eco, .data$chebi_id) %>%
    dplyr::mutate(
      evidence_source = paste0(.data$evidence_source, collapse = ";"),
      reaction = paste0(.data$reaction, collapse = ";")
    ) %>%
    dplyr::ungroup() %>%
    dplyr::distinct() %>%
    # Combine eco and eco_type
    dplyr::group_by(.data$accession, .data$chebi_id) %>%
    dplyr::mutate(
      evidence_source = paste0(.data$evidence_source, collapse = ";;"),
      eco = paste0(.data$eco, collapse = ";;"),
      eco_type = paste0(.data$eco_type, collapse = ";;"),
      reaction = paste0(.data$reaction, collapse = ";;")
    ) %>%
    dplyr::ungroup() %>%
    dplyr::distinct() %>%
    dplyr::mutate(source = "catalytic_activity") %>%
    # make sure that accession column is of type "chr" even if data frame is empty
    dplyr::mutate(accession = as.character(.data$accession))

  if (show_progress == TRUE) {
    message("DONE ", paste0("(", round(as.numeric(difftime(Sys.time(), start_time, units = "secs")), digits = 2), "s)"))
  }

  # Extract keyword information from UniProt
  if (show_progress == TRUE) {
    message("Extract keyword information from UniProt ... ", appendLF = FALSE)
    start_time <- Sys.time()
  }

  # Define additional metal names that cannot be found just based on the metal name from the element list
  additional_metal_names <- c("Metal-binding", "Metalloprotease", "2Fe-2S", "3Fe-4S", "4Fe-4S", "4Fe-4S", "Heme")

  unlisted_metal_list <- metal_list %>%
    dplyr::mutate(chebi_ion_id = stringr::str_split(.data$chebi_ion_id, pattern = ";")) %>%
    tidyr::unnest("chebi_ion_id")

  keyword_uniprot <- data_uniprot %>%
    dplyr::distinct(.data$accession, .data$keyword) %>%
    tidyr::drop_na("keyword") %>%
    dplyr::mutate(keyword = stringr::str_split(
      .data$keyword,
      pattern = ";"
    )) %>%
    tidyr::unnest("keyword") %>%
    dplyr::filter(.data$keyword %in% c(metal_list$name, additional_metal_names)) %>%
    # annotate metal_id_part and chebi_id
    dplyr::mutate(chebi_id = stats::setNames(metal_list$chebi_ion_id, metal_list$name)[.data$keyword]) %>%
    dplyr::mutate(chebi_id = dplyr::case_when(.data$keyword == "Metal-binding" ~ "25213",
                                              .data$keyword == "Metalloprotease" ~ "60240",
                                              .data$keyword == "2Fe-2S" ~ "190135",
                                              .data$keyword == "3Fe-4S" ~ "47402",
                                              .data$keyword == "4Fe-4S" ~ "49883",
                                              .data$keyword == "Heme" ~ "30413",
                                              TRUE ~ .data$chebi_id)) %>%
    dplyr::mutate(chebi_id = str_split(.data$chebi_id, pattern = ";")) %>%
    tidyr::unnest("chebi_id") %>%
    dplyr::mutate(chebi_id = ifelse(is.na(.data$chebi_id),
                                    stats::setNames(metal_list$chebi_id, metal_list$name)[.data$keyword],
                                    .data$chebi_id)) %>%
    dplyr::mutate(metal_id_part = stats::setNames(metal_chebi_uniprot$metal_atom_id, as.character(metal_chebi_uniprot$id))[.data$chebi_id]) %>%
    dplyr::mutate(metal_id_part = ifelse(is.na(.data$metal_id_part),
                                         stats::setNames(unlisted_metal_list$chebi_id, as.character(unlisted_metal_list$chebi_ion_id))[.data$chebi_id],
                                         .data$metal_id_part)) %>%
    dplyr::mutate(source = "Keyword")

  if (show_progress == TRUE) {
    message("DONE ", paste0("(", round(as.numeric(difftime(Sys.time(), start_time, units = "secs")), digits = 2), "s)"))
  }

  # Check if there are any metal containing entries that are not yet part of the ChEBI dataset provided by protti
  # This could be an indirect indication that some of the manually added ChEBI entries (without formula but metal related)
  # are also missing.
  if (any(!(unique(c(b_uniprot$chebi_id, cofactor_uniprot$chebi_id, catalytic_activity_uniprot$chebi_id)) %in% as.character(metal_chebi_uniprot$id)))) {
    missing_chebi_ids <- b_uniprot %>%
      dplyr::distinct(.data$accession, .data$chebi_id) %>%
      dplyr::bind_rows(dplyr::distinct(cofactor_uniprot, .data$accession, .data$chebi_id)) %>%
      dplyr::bind_rows(dplyr::distinct(catalytic_activity_uniprot, .data$accession, .data$chebi_id)) %>%
      dplyr::distinct() %>%
      dplyr::filter(!(.data$chebi_id %in% as.character(metal_chebi_uniprot$id)))

    warning(
      strwrap("The following ChEBI IDs have been found in the ft_binding, cc_cofactor or cc_catalytic_activity
                      column and have not yet been manually annotated in the reference data frame provided by protti.
                      This could be an indicator that there are additional ChEBI IDs missing that do not contain a
                      formula but are that are metal related IDs.
                      Please contact the package maintainer to let potentially missing ChEBI IDs be added.",
        prefix = "\n", initial = ""
      ),
      "\n",
      paste0(utils::capture.output(missing_chebi_ids), collapse = "\n")
    )
  }

  # Extract information from QuickGO
  if (show_progress == TRUE) {
    message("Extract molecular_function information from QuickGO ... ", appendLF = FALSE)
    start_time <- Sys.time()
  }

  mf_quickgo <- data_quickgo %>%
    dplyr::filter(.data$gene_product_db == "UniProtKB" &
      .data$go_aspect == "molecular_function") %>%
    # The data that will be used for filtering is a dataset provided with protti
    # It contains all MF GO IDs that are associated with metals
    dplyr::filter(.data$go_term %in% metal_go_slim_subset$slims_from_id) %>%
    # Filter GO data to only contain protein IDs also in the UniProt input
    dplyr::filter(.data$gene_product_id %in% data_uniprot$accession) %>%
    dplyr::select(
      "gene_product_id",
      "go_term",
      "go_name",
      "eco_id",
      "reference",
      "with_from",
      "assigned_by"
    ) %>%
    dplyr::distinct() %>%
    dplyr::rename(
      eco = "eco_id",
      accession = "gene_product_id"
    ) %>%
    # join ChEBI annotations to data
    dplyr::left_join(dplyr::distinct(metal_go_slim_subset, .data$slims_from_id, .data$chebi_id, .data$database, .data$metal_atom_id),
      by = c("go_term" = "slims_from_id"),
      relationship = "many-to-many"
    ) %>%
    dplyr::rename(metal_id_part = "metal_atom_id") %>%
    dplyr::mutate(eco_type = dplyr::case_when(
      .data$eco %in% manual_eco ~ "manual_assertion",
      .data$eco %in% automatic_eco ~ "automatic_assertion"
    )) %>%
    dplyr::mutate(eco_type = ifelse(is.na(.data$eco_type), "automatic_assertion", .data$eco_type)) %>%
    tidyr::unite("reference", "with_from", col = "evidence_source", na.rm = TRUE) %>%
    dplyr::mutate(evidence_source = stringr::str_replace_all(.data$evidence_source, pattern = "\\|", replacement = ",")) %>%
    # Now combine data to have one row per accession and chebi_id
    # First combine evidence_source etc.
    dplyr::group_by(.data$accession, .data$eco, .data$chebi_id) %>%
    dplyr::mutate(
      evidence_source = paste0(.data$evidence_source, collapse = ";"),
      go_term = paste0(.data$go_term, collapse = ";"),
      go_name = paste0(.data$go_name, collapse = ";"),
      assigned_by = paste0(.data$assigned_by, collapse = ";"),
      database = paste0(.data$database, collapse = ";")
    ) %>%
    dplyr::ungroup() %>%
    dplyr::distinct() %>%
    # Combine eco, eco_type etc.
    dplyr::group_by(.data$accession, .data$chebi_id) %>%
    dplyr::mutate(
      evidence_source = paste0(.data$evidence_source, collapse = ";;"),
      eco = paste0(.data$eco, collapse = ";;"),
      eco_type = paste0(.data$eco_type, collapse = ";;"),
      go_term = paste0(.data$go_term, collapse = ";;"),
      go_name = paste0(.data$go_name, collapse = ";;"),
      assigned_by = paste0(.data$assigned_by, collapse = ";;"),
      database = paste0(.data$database, collapse = ";;")
    ) %>%
    dplyr::ungroup() %>%
    dplyr::distinct() %>%
    dplyr::mutate(source = "go_term") %>%
    # make sure that accession column is of type "chr" even if data frame is empty
    dplyr::mutate(accession = as.character(.data$accession))

  if (show_progress == TRUE) {
    message("DONE ", paste0("(", round(as.numeric(difftime(Sys.time(), start_time, units = "secs")), digits = 2), "s)"))
  }

  # Find ChEBI sub IDs for IDs from cofactor, catalytic activity and GO.
  if (show_progress == TRUE) {
    message("Find ChEBI sub IDs ... ", appendLF = FALSE)
    start_time <- Sys.time()
  }

  # Metal sub IDs
  chebi_ids <- unique(c(b_uniprot$chebi_id, cofactor_uniprot$chebi_id, catalytic_activity_uniprot$chebi_id, mf_quickgo$chebi_id, keyword_uniprot$chebi_id))
  chebi_sub_id_mapping <- tibble::tibble(chebi_id = chebi_ids) %>%
    dplyr::mutate(chebi_sub_id = purrr::map_chr(find_all_subs(data_chebi_relation, .data$chebi_id, accepted_types = c("is_a", "is_conjugate_acid_of", "is_conjugate_base_of")),
      .f = ~ paste0(.x, collapse = ",")
    ))

  if (show_progress == TRUE) {
    message("DONE ", paste0("(", round(as.numeric(difftime(Sys.time(), start_time, units = "secs")), digits = 2), "s)"))
  }

  # Combine data
  if (show_progress == TRUE) {
    message("Combine data ... ", appendLF = FALSE)
    start_time <- Sys.time()
  }

  combined <- dplyr::bind_rows(b_uniprot, cofactor_uniprot, catalytic_activity_uniprot, mf_quickgo, keyword_uniprot) %>%
    dplyr::left_join(chebi_sub_id_mapping, by = "chebi_id") %>%
    dplyr::mutate(chebi_sub_id = ifelse(.data$chebi_sub_id == "", .data$chebi_id, .data$chebi_sub_id)) %>%
    dplyr::mutate(most_specific_id = stringr::str_split(.data$chebi_sub_id, pattern = ",")) %>%
    dplyr::group_by(.data$accession) %>%
    dplyr::mutate(most_specific_id = map(
      .x = .data$most_specific_id,
      .f = ~ .x[.x %in% .data$chebi_id]
    )) %>%
    dplyr::mutate(appears = unlist(map(
      .x = .data$most_specific_id,
      .f = ~ length(.x) > 1
    ))) %>%
    dplyr::ungroup() %>%
    tidyr::unnest("most_specific_id") %>%
    dplyr::arrange(.data$source) %>%
    dplyr::group_by(.data$accession, .data$most_specific_id) %>%
    dplyr::mutate(
      chebi_id_binding = ifelse(.data$source == "binding",
        .data$chebi_id,
        NA
      ),
      chebi_id = ifelse(.data$source != "binding",
        .data$chebi_id,
        NA
      )
    ) %>%
    dplyr::group_by(.data$accession, .data$most_specific_id, .data$source) %>%
    dplyr::mutate(
      chebi_id = ifelse(.data$source != "binding", paste0(unique(.data$source), "(", paste0(.data$chebi_id, collapse = "|"), ")"), NA),
      eco = ifelse(.data$source != "binding", paste0(unique(.data$source), "(", paste0(.data$eco, collapse = "|"), ")"), NA),
      eco_type = ifelse(.data$source != "binding", paste0(unique(.data$source), "(", paste0(.data$eco_type, collapse = "|"), ")"), NA),
      evidence_source = ifelse(.data$source != "binding", paste0(unique(.data$source), "(", paste0(.data$evidence_source, collapse = "|"), ")"), NA)
    ) %>%
    dplyr::group_by(.data$accession, .data$most_specific_id) %>%
    # use only metal_id_part of binding source otherwise just leave as is
    dplyr::mutate(
      metal_id_part_binding = paste0(stats::na.omit(unique(.data$metal_id_part_binding)), collapse = "|"),
      metal_id_part = paste0(stats::na.omit(unique(.data$metal_id_part)), collapse = ","),
      chebi_id_binding = paste0(stats::na.omit(unique(.data$chebi_id_binding)), collapse = "|"),
      eco_binding = paste0(stats::na.omit(unique(.data$eco_binding)), collapse = "|"),
      eco_type_binding = paste0(stats::na.omit(unique(.data$eco_type_binding)), collapse = "|"),
      evidence_source_binding = paste0(stats::na.omit(unique(.data$evidence_source_binding)), collapse = "|")
    ) %>%
    dplyr::mutate(
      chebi_id = paste0(stats::na.omit(unique(.data$chebi_id)), collapse = ","),
      eco = paste0(stats::na.omit(unique(.data$eco)), collapse = ","),
      eco_type = paste0(stats::na.omit(unique(.data$eco_type)), collapse = ","),
      evidence_source = paste0(stats::na.omit(unique(.data$evidence_source)), collapse = ","),
      source = paste0(unique(.data$source), collapse = " & "),
      ligand_identifier = paste0(stats::na.omit(.data$ligand_identifier), collapse = "|"),
      ligand_position = paste0(stats::na.omit(.data$ligand_position), collapse = "|"),
      binding_mode = paste0(stats::na.omit(.data$binding_mode), collapse = "|"),
      metal_function = paste0(stats::na.omit(.data$metal_function), collapse = "|"),
      note = paste0(stats::na.omit(.data$note), collapse = "|"),
      reaction = paste0(stats::na.omit(.data$reaction), collapse = "|"),
      go_term = paste0(stats::na.omit(.data$go_term), collapse = "|"),
      go_name = paste0(stats::na.omit(.data$go_name), collapse = "|"),
      assigned_by = paste0(stats::na.omit(.data$assigned_by), collapse = "|"),
      database = paste0(stats::na.omit(.data$database), collapse = "|"),
      keyword = paste0(stats::na.omit(.data$keyword), collapse = "|")
    ) %>%
    dplyr::ungroup() %>%
    dplyr::filter(!.data$appears) %>%
    dplyr::select(-c("appears", "chebi_sub_id")) %>%
    dplyr::distinct() %>%
    # unpack positions and corresponding info
    dplyr::mutate(
      ligand_identifier = stringr::str_split(.data$ligand_identifier, pattern = "\\|"),
      ligand_position = stringr::str_split(.data$ligand_position, pattern = "\\|"),
      metal_id_part_binding = stringr::str_split(.data$metal_id_part_binding, pattern = "\\|"),
      binding_mode = stringr::str_split(.data$binding_mode, pattern = "\\|"),
      metal_function = stringr::str_split(.data$metal_function, pattern = "\\|"),
      eco_binding = stringr::str_split(.data$eco_binding, pattern = "\\|"),
      eco_type_binding = stringr::str_split(.data$eco_type_binding, pattern = "\\|"),
      evidence_source_binding = stringr::str_split(.data$evidence_source_binding, pattern = "\\|"),
      chebi_id_binding = stringr::str_split(.data$chebi_id_binding, pattern = "\\|")
    ) %>%
    tidyr::unnest(c(
      "ligand_identifier",
      "ligand_position",
      "metal_id_part_binding",
      "binding_mode",
      "metal_function",
      "eco_binding",
      "eco_type_binding",
      "evidence_source_binding",
      "chebi_id_binding"
    )) %>%
    dplyr::mutate(
      ligand_identifier = stringr::str_split(.data$ligand_identifier, pattern = ";;;;"),
      ligand_position = stringr::str_split(.data$ligand_position, pattern = ";;;;"),
      metal_id_part_binding = stringr::str_split(.data$metal_id_part_binding, pattern = ";;;;"),
      binding_mode = stringr::str_split(.data$binding_mode, pattern = ";;;;"),
      metal_function = stringr::str_split(.data$metal_function, pattern = ";;;;"),
      eco_binding = stringr::str_split(.data$eco_binding, pattern = ";;;;"),
      eco_type_binding = stringr::str_split(.data$eco_type_binding, pattern = ";;;;"),
      evidence_source_binding = stringr::str_split(.data$evidence_source_binding, pattern = ";;;;"),
      chebi_id_binding = stringr::str_split(.data$chebi_id_binding, pattern = ";;;;")
    ) %>%
    tidyr::unnest(c(
      "ligand_identifier",
      "ligand_position",
      "metal_id_part_binding",
      "binding_mode",
      "metal_function",
      "eco_binding",
      "eco_type_binding",
      "evidence_source_binding",
      "chebi_id_binding"
    )) %>%
    dplyr::mutate(
      ligand_identifier = stringr::str_split(.data$ligand_identifier, pattern = ";;;"),
      ligand_position = stringr::str_split(.data$ligand_position, pattern = ";;;"),
      metal_id_part_binding = stringr::str_split(.data$metal_id_part_binding, pattern = ";;;"),
      binding_mode = stringr::str_split(.data$binding_mode, pattern = ";;;"),
      metal_function = stringr::str_split(.data$metal_function, pattern = ";;;"),
      eco_binding = stringr::str_split(.data$eco_binding, pattern = ";;;"),
      eco_type_binding = stringr::str_split(.data$eco_type_binding, pattern = ";;;"),
      evidence_source_binding = stringr::str_split(.data$evidence_source_binding, pattern = ";;;"),
      chebi_id_binding = stringr::str_split(.data$chebi_id_binding, pattern = ";;;")
    ) %>%
    tidyr::unnest(c(
      "ligand_identifier",
      "ligand_position",
      "metal_id_part_binding",
      "binding_mode",
      "metal_function",
      "eco_binding",
      "eco_type_binding",
      "evidence_source_binding",
      "chebi_id_binding"
    )) %>%
    dplyr::mutate(ligand_position = as.numeric(.data$ligand_position)) %>%
    dplyr::arrange(.data$accession, .data$ligand_position) %>%
    dplyr::left_join(chebi_names,
      by = c("most_specific_id" = "id")
    ) %>%
    dplyr::rename(most_specific_id_name = "name") %>%
    dplyr::mutate(
      note = ifelse(.data$note == "NA" | .data$note == "", NA, .data$note),
      reaction = ifelse(.data$reaction == "NA" | .data$reaction == "", NA, .data$reaction),
      go_term = ifelse(.data$go_term == "NA" | .data$go_term == "", NA, .data$go_term),
      go_name = ifelse(.data$go_name == "NA" | .data$go_name == "", NA, .data$go_name),
      assigned_by = ifelse(.data$assigned_by == "NA" | .data$assigned_by == "", NA, .data$assigned_by),
      database = ifelse(.data$database == "NA" | .data$database == "", NA, .data$database),
      keyword = ifelse(.data$keyword == "NA" | .data$keyword == "", NA, .data$keyword)
    ) %>%
    dplyr::mutate(binding_temp_1 = stringr::str_split(.data$metal_id_part, pattern = ",")) %>%
    tidyr::unnest("binding_temp_1") %>%
    dplyr::left_join(chebi_names,
      by = c("binding_temp_1" = "id")
    ) %>%
    dplyr::rename(metal_id_part_name = "name") %>%
    dplyr::select(-c("binding_temp_1")) %>%
    dplyr::group_by(.data$accession, .data$chebi_id, .data$ligand_identifier, .data$ligand_position) %>%
    dplyr::mutate(metal_id_part_name = paste0(.data$metal_id_part_name, collapse = ",")) %>%
    dplyr::ungroup() %>%
    dplyr::distinct() %>%
    dplyr::mutate(binding_temp_2 = stringr::str_split(.data$metal_id_part_binding, pattern = ",")) %>%
    tidyr::unnest("binding_temp_2") %>%
    dplyr::left_join(chebi_names,
      by = c("binding_temp_2" = "id")
    ) %>%
    dplyr::rename(metal_id_part_binding_name = "name") %>%
    dplyr::select(-c("binding_temp_2")) %>%
    dplyr::group_by(.data$accession, .data$chebi_id, .data$ligand_identifier, .data$ligand_position) %>%
    dplyr::mutate(metal_id_part_binding_name = paste0(.data$metal_id_part_binding_name, collapse = ",")) %>%
    dplyr::distinct() %>%
    dplyr::group_by(.data$accession, .data$ligand_identifier, .data$ligand_position) %>%
    dplyr::mutate(
      most_specific_id = paste0(.data$most_specific_id, collapse = ","),
      most_specific_id_name = paste0(.data$most_specific_id_name, collapse = "||"),
      source = paste0(.data$source, collapse = "||"),
      chebi_id = paste0(.data$chebi_id, collapse = "||"),
      note = paste0(.data$note, collapse = "||"),
      eco = paste0(.data$eco, collapse = "||"),
      eco_type = paste0(.data$eco_type, collapse = "||"),
      evidence_source = paste0(.data$evidence_source, collapse = "||"),
      reaction = paste0(.data$reaction, collapse = "||"),
      go_term = paste0(.data$go_term, collapse = "||"),
      go_name = paste0(.data$go_name, collapse = "||"),
      assigned_by = paste0(.data$assigned_by, collapse = "||"),
      database = paste0(.data$database, collapse = "||"),
      keyword = paste0(.data$keyword, collapse = "||"),
      metal_id_part = paste0(unique(.data$metal_id_part), collapse = ","),
      metal_id_part_name = paste0(unique(.data$metal_id_part_name), collapse = ",")
    ) %>%
    dplyr::ungroup() %>%
    dplyr::distinct() %>%
    dplyr::mutate(
      note = ifelse(.data$note == "NA" | .data$note == "", NA, .data$note),
      reaction = ifelse(.data$reaction == "NA" | .data$reaction == "", NA, .data$reaction),
      go_term = ifelse(.data$go_term == "NA" | .data$go_term == "", NA, .data$go_term),
      go_name = ifelse(.data$go_name == "NA" | .data$go_name == "", NA, .data$go_name),
      assigned_by = ifelse(.data$assigned_by == "NA" | .data$assigned_by == "", NA, .data$assigned_by),
      database = ifelse(.data$database == "NA" | .data$database == "", NA, .data$database),
      keyword = ifelse(.data$keyword == "NA" | .data$keyword == "", NA, .data$keyword),
    ) %>%
    dplyr::mutate(
      metal_id_part_binding = ifelse(.data$metal_id_part_binding == "NA" | .data$metal_id_part_binding == "", NA, .data$metal_id_part_binding),
      metal_id_part_binding_name = ifelse(.data$metal_id_part_binding_name == "NA" | .data$metal_id_part_binding_name == "", NA, .data$metal_id_part_binding_name),
      metal_id_part = ifelse(.data$metal_id_part == "NA" | .data$metal_id_part == "", NA, .data$metal_id_part),
      metal_id_part_name = ifelse(.data$metal_id_part_name == "NA" | .data$metal_id_part_name == "", NA, .data$metal_id_part_name),
      binding_mode = ifelse(.data$binding_mode == "NA" | .data$binding_mode == "", NA, .data$binding_mode),
      metal_function = ifelse(.data$metal_function == "NA" | .data$metal_function == "", NA, .data$metal_function),
      ligand_identifier = ifelse(.data$ligand_identifier == "NA" | .data$ligand_identifier == "", NA, .data$ligand_identifier)
    ) %>%
    dplyr::mutate(
      chebi_id = ifelse(.data$chebi_id_binding == "",
        .data$chebi_id,
        paste0("binding(", .data$chebi_id_binding, "),", .data$chebi_id)
      ),
      eco = ifelse(.data$eco_binding == "",
        .data$eco,
        paste0("binding(", .data$eco_binding, "),", .data$eco)
      ),
      eco_type = ifelse(.data$eco_binding == "",
        .data$eco_type,
        paste0("binding(", .data$eco_type_binding, "),", .data$eco_type)
      ),
      evidence_source = ifelse(.data$eco_binding == "",
        .data$evidence_source,
        paste0("binding(", .data$evidence_source_binding, "),", .data$evidence_source)
      )
    ) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      metal_id_part = ifelse(is.na(.data$metal_id_part_binding) & !is.na(.data$ligand_position),
        NA,
        paste0(stats::na.omit(unique(c(unlist(stringr::str_split(.data$metal_id_part_binding, pattern = ",")), unlist(stringr::str_split(.data$metal_id_part, pattern = ","))))), collapse = ",")
      ),
      metal_id_part_name = ifelse(is.na(.data$metal_id_part_binding_name) & !is.na(.data$ligand_position),
        NA,
        paste0(stats::na.omit(unique(c(unlist(stringr::str_split(.data$metal_id_part_binding_name, pattern = ",")), unlist(stringr::str_split(.data$metal_id_part_name, pattern = ","))))), collapse = ",")
      )
    ) %>%
    dplyr::ungroup() %>%
    dplyr::select(-c("eco_binding", "eco_type_binding", "evidence_source_binding", "chebi_id_binding", "metal_id_part_binding", "metal_id_part_binding_name")) %>%
    dplyr::select(
      "accession",
      "most_specific_id",
      "most_specific_id_name",
      "ligand_identifier",
      "ligand_position",
      "binding_mode",
      "metal_function",
      "metal_id_part",
      "metal_id_part_name",
      "note",
      "chebi_id",
      "source",
      "eco",
      "eco_type",
      "evidence_source",
      "reaction",
      "go_term",
      "go_name",
      "assigned_by",
      "database",
      "keyword"
    )

  if (show_progress == TRUE) {
    message("DONE ", paste0("(", round(as.numeric(difftime(Sys.time(), start_time, units = "secs")), digits = 2), "s)"))
  }

  return(combined)
}
