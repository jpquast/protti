#' Creates a contact map of all atoms from a structure file
#'
#' Creates a contact map of a subset or of all atom or residue distances in a structure or AlphaFold prediction file. Contact maps are a useful tool for the identification of protein
#' regions that are in close proximity in the folded protein. Additionally, regions that are interacting closely with a small molecule or metal ion
#' can be easily identified without the need to open the structure in programs such as PyMOL or ChimeraX.
#'
#' @param data a data frame containing at least a column with PDB ID information of which the name can be provided to the \code{id} argument.
#' If only this column is provided, all atom or residue distances are calculated. Additionally, a chain column can be present in the data frame of which
#' the name can be provided to the \code{chain} argument. If chains are provided, only distances of this chain relative to the rest of the structure
#' are calculated. Multiple chains can be provided in multiple rows. If chains are provided for one structure but not for another, the rows should contain
#' NAs. Furthermore, specific residue positions can be provided in start and end columns if the selection should be further reduced. It is not
#' recommended to create full contact maps for more than a few structures due to time and memory reasons. If contact maps are created only for
#' small regions it is possible to create multiple maps at once.
#' @param id a column in the \code{data} data frame that contains PDB or UniProt IDs for structures or AlphaFold predictions of which contact maps
#' should be created. If a structure not downloaded directly from PDB is provided (i.e. a locally stored structure file) to the \code{structure_file}
#' argument, this column should contain "my_structure" as content.
#' @param chain an optional column in the \code{data} data frame that contains chain identifiers for the structure file. Identifiers defined by the structure
#' author should be used. Distances will be only calculated between the provided chains and the rest of the structure.
#' @param start_in_pdb an optional column in the \code{data} data frame that contains start positions of regions which for distances should be calculated.
#' This needs to be always provided in combination with a corresponding end position in \code{end_in_pdb} and chain in \code{chain}. The position should
#' match the positioning defined by the structure author. For PDB structures this information can be obtained from the
#' \code{find_peptide_in_structure} function. The corresponding column in the output is called \code{auth_seq_id_start}. If an AlphaFold prediction is
#' provided, UniProt positions should be used.
#' @param end_in_pdb an optional column in the \code{data} data frame that contains end positions of regions which for distances should be calculated.
#' This needs to be always provided in combination with a corresponding start position in \code{start_in_pdb} and chain in \code{chain}. The position should
#' match the positioning defined by the structure author. For PDB structures this information can be obtained from the \code{find_peptide_in_structure} function. The
#' corresponding column in the output is called \code{auth_seq_id_end}. If an AlphaFold prediction is provided, UniProt positions should be used.
#' @param distance_cutoff a numeric value specifying the distance cutoff in Angstrom. All values for pairwise comparisons are calculated but only values
#' smaller than this cutoff will be returned in the output. If a cutoff of e.g. 5 is selected then only residues with a distance of 5 Angstrom and less are returned.
#' Using a small value can reduce the size of the contact map drastically and is therefore recommended. The default value is 10.
#' @param pdb_model_number_selection a numeric vector specifying which models from the structure files should be considered for contact maps. E.g. NMR
#' models often have many models in one file. The default for this argument is c(0, 1). This means the first model of each structure file is selected for contact map
#' calculations. For AlphaFold predictions the model number is 0 (only .pdb files), therefore this case is also included here.
#' @param return_min_residue_distance a logical, specifying if the contact map should be returned for all atom distances or the minimum residue distances.
#' Minimum residue distances are smaller in size. If atom distances are not strictly needed it is recommended to set this argument to TRUE. The default is TRUE.
#' @param show_progress a logical, if \code{show_progress = TRUE}, a progress bar will be shown (default is TRUE).
#' @param export a logical indicating if contact maps should be exported as ".csv". The name of the file will be the structure ID. Default is \code{export = FALSE}.
#' @param export_location optional, a character argument specifying the path to the location in which the contact map should be saved if \code{export = TRUE}. If left empty, they will be saved in the current
#' working directory. The location should be provided in the following format "folderA/folderB".
#' @param structure_file optional, a character argument specifying the path to the location and name of a structure file in ".cif" or ".pdb" format
#' for which a contact map should be created. All other arguments can be provided as usual with the exception of the \code{id} column in the
#' \code{data} data frame, which should not contain a PDB or UniProt ID but a character vector containing only "my_structure".
#'
#' @return A list of contact maps for each PDB or UniProt ID provided in the input is returned. If the \code{export} argument is TRUE, each contact map will
#' be saved as a ".csv" file in the current working directory or the location provided to the \code{export_location} argument.
#' @import dplyr
#' @import tidyr
#' @import progress
#' @importFrom purrr map2 map map_dfr
#' @importFrom readr read_tsv write_csv
#' @importFrom stringr str_replace_all str_sub str_detect str_extract str_replace
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @export
#'
#' @examples
#' \dontrun{
#' create_structure_contact_map(
#'   data = data,
#'   id = pdb_id,
#'   chain = auth_asym_id,
#'   start_in_pdb = auth_seq_id_start,
#'   end_in_pdb = auth_seq_id_end
#' )
#' }
create_structure_contact_map <- function(data, id, chain = NULL, start_in_pdb = NULL, end_in_pdb = NULL, distance_cutoff = 10, pdb_model_number_selection = c(0, 1), return_min_residue_distance = TRUE, show_progress = TRUE, export = FALSE, export_location = NULL, structure_file = NULL) {
  ids <- unique(dplyr::pull(data, {{ id }}))

  # assign all protein IDs to be chain A since AlphaFold only contains one chain which is always A.
  if (!missing(chain)) {
    data <- data %>%
      dplyr::mutate({{ chain }} := ifelse(stringr::str_detect({{ id }}, pattern = "^([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})$"), "A", {{ chain }}))
  }

  if (ifelse(!missing(start_in_pdb) & (missing(chain) || all(is.na(dplyr::pull(data, {{ chain }})))),
    any(!is.na(dplyr::pull(data, {{ start_in_pdb }}))),
    FALSE
  ) |
    ifelse(!missing(start_in_pdb) & !missing(chain),
      any(!is.na(dplyr::pull(data, {{ start_in_pdb }})[is.na(dplyr::pull(data, {{ chain }}))])),
      FALSE
    )) {
    stop("The data contains start and end positions whithout specified chain IDs. Please always provide a chain ID for your start and end positions.")
  }

  # create individual data retain patterns depending on which information was provided to the function.
  if (missing(chain) || all(is.na(dplyr::pull(data, {{ chain }})))) {
    data_retain_pattern <- data %>%
      dplyr::ungroup() %>%
      dplyr::mutate(retain_pattern = {{ id }}) %>%
      dplyr::pull(.data$retain_pattern) %>%
      unique()
  } else {
    if (missing(start_in_pdb) || all(is.na(dplyr::pull(data, {{ chain }})))) {
      data_retain_pattern <- data %>%
        dplyr::ungroup() %>%
        dplyr::mutate(retain_pattern = stringr::str_replace_all(paste({{ id }}, {{ chain }}, sep = "_"), pattern = "_NA", replacement = "")) %>%
        dplyr::pull(.data$retain_pattern) %>%
        unique()
    } else {
      data_retain_pattern <- data %>%
        dplyr::ungroup() %>%
        dplyr::mutate(
          start = {{ start_in_pdb }},
          end = {{ end_in_pdb }}
        ) %>%
        # do this so start and end position can be the same column.
        dplyr::distinct({{ id }}, {{ chain }}, .data$start, .data$end) %>%
        group_by({{ id }}, {{ chain }}, .data$start, .data$end) %>%
        dplyr::mutate(residue = ifelse(!is.na(.data$start), list(seq(.data$start, .data$end)), list(NA))) %>%
        dplyr::ungroup() %>%
        tidyr::unnest(.data$residue) %>%
        dplyr::mutate(retain_pattern = stringr::str_replace_all(paste({{ id }}, {{ chain }}, .data$residue, sep = "_"), pattern = "_NA", replacement = "")) %>%
        dplyr::pull(.data$retain_pattern) %>%
        unique()
    }
  }

  if (!missing(structure_file)) {
    file_format <- stringr::str_sub(structure_file, start = -4, end = -1)

    if (!file_format %in% c(".cif", ".pdb")) {
      stop('Please either provide a ".cif" or ".pdb" structure file.')
    }

    file_name <- stringr::str_extract(structure_file, pattern = "[^/]+[:punct:]\\w+$")

    structure_file <- readr::read_tsv(structure_file, quote = "", col_names = FALSE, show_col_types = FALSE, progress = FALSE)

    # load .cif file if provided
    if (file_format == ".cif") {
      structure <- structure_file %>%
        dplyr::filter(stringr::str_detect(.data$X1, pattern = "^ATOM|^HETATM")) %>%
        dplyr::mutate(X2 = stringr::str_replace_all(.data$X1, pattern = "\\s+", replacement = " ")) %>%
        tidyr::separate(.data$X2,
          sep = " ",
          into = c("x1", "label_id", "type_symbol", "label_atom_id", "x2", "label_comp_id", "label_asym_id", "entity_id", "label_seq_id", "x3", "x", "y", "z", "site_occupancy", "b_iso_or_equivalent", "formal_charge", "auth_seq_id", "auth_comp_id", "auth_asym_id", "x4", "pdb_model_number"),
          extra = "drop"
        ) %>%
        dplyr::select(-c(.data$X1, .data$x1, .data$x2, .data$x3, .data$x4)) %>%
        dplyr::group_by(.data$label_asym_id, .data$label_atom_id, .data$label_comp_id) %>%
        dplyr::mutate(label_seq_id = ifelse(.data$label_seq_id == ".", 1:n(), as.numeric(.data$label_seq_id))) %>%
        dplyr::ungroup() %>%
        dplyr::mutate(
          label_id = as.numeric(.data$label_id),
          label_seq_id = as.numeric(.data$label_seq_id),
          x = as.numeric(.data$x),
          y = as.numeric(.data$y),
          z = as.numeric(.data$z),
          b_iso_or_equivalent = as.numeric(.data$b_iso_or_equivalent),
          auth_seq_id = as.numeric(.data$auth_seq_id),
          pdb_model_number = as.numeric(.data$pdb_model_number),
          id = "my_structure"
        ) %>%
        dplyr::filter(.data$pdb_model_number %in% pdb_model_number_selection) %>%
        dplyr::select(
          .data$label_id,
          .data$x,
          .data$y,
          .data$z,
          .data$label_comp_id,
          .data$label_seq_id,
          .data$label_asym_id,
          .data$auth_comp_id,
          .data$auth_seq_id,
          .data$auth_asym_id,
          .data$id
        ) %>%
        dplyr::mutate(retain_pattern = stringr::str_replace_all(paste(.data$id, .data$auth_asym_id, .data$auth_seq_id, sep = "_"), pattern = "_NA", replacement = "")) %>%
        dplyr::mutate(should_be_retained = stringr::str_detect(.data$retain_pattern, pattern = paste(paste0(data_retain_pattern, "(?=$|_)"), collapse = "|")))
    }
    # load .pdb file if provided
    if (file_format == ".pdb") {
      structure <- structure_file %>%
        dplyr::filter(stringr::str_detect(.data$X1, pattern = "^ATOM|^HETATM|^MODEL")) %>%
        dplyr::mutate(pdb_model_number = as.numeric(ifelse(stringr::str_detect(.data$X1, pattern = "^MODEL"), stringr::str_extract(.data$X1, pattern = "\\d+"), NA))) %>%
        tidyr::fill(.data$pdb_model_number, .direction = "down") %>%
        dplyr::mutate(pdb_model_number = ifelse(is.na(.data$pdb_model_number), 0, .data$pdb_model_number)) %>%
        dplyr::filter(!stringr::str_detect(.data$X1, pattern = "^MODEL")) %>%
        dplyr::mutate(
          label_id = as.numeric(stringr::str_replace_all(stringr::str_sub(.data$X1, start = 7, end = 11), pattern = "\\s+", replacement = "")),
          auth_comp_id = stringr::str_replace_all(stringr::str_sub(.data$X1, start = 18, end = 20), pattern = "\\s+", replacement = ""),
          auth_asym_id = stringr::str_replace_all(stringr::str_sub(.data$X1, start = 22, end = 22), pattern = "\\s+", replacement = ""),
          auth_seq_id = as.numeric(stringr::str_replace_all(suppressWarnings(as.numeric(stringr::str_sub(.data$X1, start = 23, end = 26))), pattern = "\\s+", replacement = "")),
          x = as.numeric(stringr::str_replace_all(stringr::str_sub(.data$X1, start = 31, end = 38), pattern = "\\s+", replacement = "")),
          y = as.numeric(stringr::str_replace_all(stringr::str_sub(.data$X1, start = 39, end = 46), pattern = "\\s+", replacement = "")),
          z = as.numeric(stringr::str_replace_all(stringr::str_sub(.data$X1, start = 47, end = 54), pattern = "\\s+", replacement = "")),
          id = "my_structure"
        ) %>%
        dplyr::filter(.data$pdb_model_number %in% pdb_model_number_selection) %>%
        dplyr::select(-c(.data$X1, .data$pdb_model_number)) %>%
        dplyr::mutate(retain_pattern = stringr::str_replace_all(paste(.data$id, .data$auth_asym_id, .data$auth_seq_id, sep = "_"), pattern = "_NA", replacement = "")) %>%
        dplyr::mutate(should_be_retained = stringr::str_detect(.data$retain_pattern, pattern = paste(paste0(data_retain_pattern, "(?=$|_)"), collapse = "|")))
    }

    structures <- list(structure)
    names(structures) <- stringr::str_replace(file_name, pattern = file_format, replacement = "")
  }

  if (missing(structure_file)) {
    if (!requireNamespace("httr", quietly = TRUE)) {
      stop("Package \"httr\" is needed for this function to work. Please install it.", call. = FALSE)
    }

    if (!curl::has_internet()) {
      message("No internet connection.")
      return(invisible(NULL))
    }

    # make ID vectors

    pdb_ids <- ids[nchar(ids) == 4]
    uniprot_ids <- ids[nchar(ids) != 4]

    # placeholders
    pdb_structures <- NULL
    alphafold_structures <- NULL

    if (length(pdb_ids) != 0) {
      if (show_progress == TRUE) {
        pb <- progress::progress_bar$new(total = length(pdb_ids), format = "Preparing structures [:bar] :current/:total (:percent) :eta")
      }

      pdb_structures <- fetch_pdb_structure(pdb_ids = pdb_ids, return_data_frame = FALSE, show_progress = show_progress) %>%
        purrr::map(.f = ~ {
          if (show_progress == TRUE) {
            pb$tick()
          }
          .x %>%
            dplyr::filter(.data$pdb_model_number %in% pdb_model_number_selection) %>%
            dplyr::select(
              .data$label_id,
              .data$x,
              .data$y,
              .data$z,
              .data$label_comp_id,
              .data$label_seq_id,
              .data$label_asym_id,
              .data$auth_comp_id,
              .data$auth_seq_id,
              .data$auth_asym_id,
              .data$pdb_id
            ) %>%
            dplyr::mutate(retain_pattern = stringr::str_replace_all(paste(.data$pdb_id, .data$auth_asym_id, .data$auth_seq_id, sep = "_"), pattern = "_NA", replacement = "")) %>%
            dplyr::mutate(should_be_retained = stringr::str_detect(.data$retain_pattern, pattern = paste(paste0(data_retain_pattern, "(?=$|_)"), collapse = "|"))) %>%
            dplyr::rename(id = .data$pdb_id)
        })
    }
    if (length(uniprot_ids) != 0) {
      if (show_progress == TRUE) {
        pb <- progress::progress_bar$new(total = length(uniprot_ids), format = "Preparing pedictions [:bar] :current/:total (:percent) :eta")
      }

      alphafold_structures <- fetch_alphafold_prediction(uniprot_ids = uniprot_ids, return_data_frame = FALSE, show_progress = show_progress) %>%
        purrr::map(.f = ~ {
          if (show_progress == TRUE) {
            pb$tick()
          }
          .x %>%
            dplyr::select(
              .data$label_id,
              .data$x,
              .data$y,
              .data$z,
              .data$label_comp_id,
              .data$label_seq_id,
              .data$label_asym_id,
              .data$auth_comp_id,
              .data$auth_seq_id,
              .data$auth_asym_id,
              .data$uniprot_id,
              .data$prediction_score,
              .data$score_quality
            ) %>%
            dplyr::mutate(retain_pattern = stringr::str_replace_all(paste(.data$uniprot_id, .data$auth_asym_id, .data$auth_seq_id, sep = "_"), pattern = "_NA", replacement = "")) %>%
            dplyr::mutate(should_be_retained = stringr::str_detect(.data$retain_pattern, pattern = paste(paste0(data_retain_pattern, "(?=$|_)"), collapse = "|"))) %>%
            dplyr::rename(id = .data$uniprot_id)
        })
    }
    structures <- c(pdb_structures, alphafold_structures)
  }

  # this data frame contains the subsetted structures. Specifically, the subsetted atom numbers.
  if (show_progress == TRUE) {
    pb <- progress::progress_bar$new(total = length(structures), format = "Subsetting structures [:bar] :current/:total (:percent) :eta")
  }
  subset_structures <- structures %>%
    purrr::map(.f = ~ {
      if (show_progress == TRUE) {
        pb$tick()
      }
      .x %>%
        dplyr::filter(should_be_retained) %>%
        dplyr::distinct(.data$label_id)
    })

  # Segments are made for each structure to prevent too long data frames when all combinations are created.
  segments <- subset_structures %>%
    purrr::map(.f = ~ {
      split(dplyr::pull(.x, .data$label_id), ceiling(seq_along(nrow(.x)) / 1000))
    })

  if (show_progress == TRUE) {
    pb <- progress::progress_bar$new(total = length(structures), format = "Calculating atom distances :current/:total (:percent)")
  }

  result_distances <- purrr::map2(
    .x = segments,
    .y = structures,
    .f = ~ {
      if (show_progress == TRUE) {
        pb$tick()
      }

      current_structure <- .y %>%
        dplyr::select(-c(.data$x, .data$y, .data$z, .data$should_be_retained, .data$retain_pattern))

      current_structure_minimum <- .y %>%
        dplyr::distinct(.data$label_id, .data$x, .data$y, .data$z)

      current_protein <- .y %>%
        dplyr::distinct(.data$id) %>%
        dplyr::pull(.data$id)

      if (show_progress == TRUE) {
        pb <- progress::progress_bar$new(total = length(.x), format = paste("Calculating distances for", current_protein, "[:bar] :current/:total (:percent)"))
      }

      purrr::map_dfr(
        .x = .x,
        .f = ~ {
          if (show_progress == TRUE) {
            pb$tick()
          }

          tidyr::crossing(var1 = .x, var2 = current_structure_minimum$label_id) %>%
            dplyr::left_join(current_structure_minimum, by = c("var1" = "label_id")) %>%
            dplyr::left_join(current_structure_minimum, by = c("var2" = "label_id")) %>%
            dplyr::mutate(distance = sqrt((.data$x.x - .data$x.y)^2 + (.data$y.x - .data$y.y)^2 + (.data$z.x - .data$z.y)^2)) %>%
            dplyr::select(.data$var1, .data$var2, .data$distance) %>%
            dplyr::filter(.data$distance <= distance_cutoff) %>%
            dplyr::left_join(current_structure %>% dplyr::select(-.data$id), by = c("var1" = "label_id")) %>%
            dplyr::left_join(current_structure, by = c("var2" = "label_id"), suffix = c("_var1", "_var2"))
        }
      )
    }
  )

  if (show_progress == TRUE) {
    pb <- progress::progress_bar$new(total = length(result_distances), format = "Calculating minimal residue distances [:bar] :current/:total (:percent) :eta")
  }
  # calculate residue distances only after the table is complete, otherwise wrong distances might be calculated.
  result <- result_distances %>%
    purrr::map(.f = ~ {
      if (show_progress == TRUE) {
        pb$tick()
      }
      residue_distance <- .x %>%
        dplyr::group_by(.data$auth_seq_id_var1, .data$auth_asym_id_var1, .data$auth_seq_id_var2, .data$auth_asym_id_var2) %>%
        dplyr::mutate(min_distance_residue = min(.data$distance)) %>%
        dplyr::ungroup() %>%
        dplyr::rename(
          label_id_var1 = .data$var1,
          label_id_var2 = .data$var2
        )

      if (return_min_residue_distance == TRUE) {
        residue_distance <- residue_distance %>%
          dplyr::select(-c(.data$label_id_var1, .data$label_id_var2, .data$distance)) %>%
          dplyr::distinct()
      }

      residue_distance
    })

  if (export == FALSE) {
    return(result)
  } else {
    # make sure export location is correct if or if not provided.
    if (missing(export_location)) {
      export_location <- ""
    } else {
      export_location <- paste0(export_location, "/")
    }
    if (show_progress == TRUE) {
      pb <- progress::progress_bar$new(total = length(result), format = "Exporting contact maps [:bar] :current/:total (:percent) :eta")
    }
    purrr::map2(
      .x = result,
      .y = names(result),
      .f = ~ {
        if (show_progress == TRUE) {
          pb$tick()
        }
        readr::write_csv(x = .x, file = paste0(export_location, .y, "_contact_map.csv"), progress = FALSE)
      }
    )
    return(invisible(NULL))
  }
}
