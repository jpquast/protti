#' Creates a contact map of all atoms from a structure file
#'
#' Creates a contact map of a subset or of all atom or residue distances in a structure file. Contact maps are a useful tool for the identification of protein
#' regions that are in close proximity in the folded protein. Additionally, regions that are interacting closely with a small molecule or metal ion
#' can be easily identified without the need to open the structure in programs such as PyMOL or ChimeraX.
#'
#' @param data a data frame containing at least a column with PDB ID information of which the name can be provided to the \code{pdb_id} argument.
#' If only this column is provided, all atom or residue distances are calculated. Additionally, a chain column can be present in the data frame of which
#' the name can be provided to the \code{chain} argument. If chains are provided, only distances of this chain relative to the rest of the structure
#' are calculated. Multiple chains can be provided in multiple rows. If chains are provided for one structure but not for another, the rows should contain
#' NAs. Furthermore, specific residue positions can be provided in start and end columns if the selection should be further reduced. It is not
#' recommended to create full contact maps for more than a few structures due to time and memory reasons. If contact maps are created only for
#' small regions it is possible to create multiple maps at once.
#' @param pdb_id a column in the \code{data} data frame that contains PDB IDs for structures of which contact maps should be created. If an own structure
#' is provided to the \code{structure_file} argument, this column should contain "my_structure" as content.
#' @param chain an optional column in the \code{data} data frame that contains chain identifiers for the structure file. The author defined identifiers
#' should be used. Distances will be only calculated between the provided chains and the rest of the structure.
#' @param start_in_pdb an optional column in the \code{data} data frame that contains start positions of regions which for distances should be calculated.
#' This needs to be always provided in combination with a corresponding end position in \code{end_in_pdb} and chain in \code{chain}. The position should
#' match the positioning based on author defined positioning. This information can be obtained from the \code{find_peptide_in_pdb} function. The
#' corresponding column in the output is called \code{peptide_start_pdb}.
#' @param end_in_pdb an optional column in the \code{data} data frame that contains end positions of regions which for distances should be calculated.
#' This needs to be always provided in combination with a corresponding start position in \code{start_in_pdb} and chain in \code{chain}. The position should
#' match the positioning based on author defined positioning. This information can be obtained from the \code{find_peptide_in_pdb} function. The
#' corresponding column in the output is called \code{peptide_start_pdb}.
#' @param distance_cutoff a numeric value specifying the distance cutoff in Angstrom. All values greater than the cutoff will be discarded. This can reduce
#' the size of the contact map drastically and is therefore recommended. The default value is 10.
#' @param pdb_model_number_selection a numeric vector specifying which models from the structure files should be considered for contact maps. E.g. NMR
#' models often have many models in one file. The default for this argument is 1. This means only the first model of each structure file is selected for contact map
#' calculations.
#' @param return_min_residue_distance a logical, specifying if the contact map should be returned for all atom distances or the minimum residue distances.
#' Minimum residue distances are smaller in size. If atom distances are not strictly needed it is recommended to set this argument to TRUE. The default is TRUE.
#' @param show_progress a logical, if \code{show_progress = TRUE}, a progress bar will be shown (default is TRUE).
#' @param export a logical indicating if contact maps should be exported as ".csv". The name of the file will be the structure ID. Default is \code{export = FALSE}.
#' @param export_location optional, a character argument specifying the path to the location in which the contact map should be saved if \code{export = TRUE}. If left empty, they will be saved in the current
#' working directory. The location should be provided in the following format "folderA/folderB".
#' @param structure_file optional, a character argument specifying the path to the location and name of a structure file in ".cif" or ".pdb" format
#' for which a contact map should be created. All other arguments can be provided as usual with the exception of the \code{pdb_id} column in the
#' \code{data} data frame, which should not contain a PDB ID but a character vector containing only "my_structure".
#'
#' @return A list of contact maps for each structure ID provided in the input is returned. If the \code{export} argument is TRUE, each contact map will
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
#'   pdb_id = pdb_id,
#'   chain = auth_chain,
#'   start_in_pdb = peptide_start_pdb,
#'   end_in_pdb = peptide_end_pdb
#' )
#' }
create_structure_contact_map <- function(data, pdb_id, chain = NULL, start_in_pdb = NULL, end_in_pdb = NULL, distance_cutoff = 10, pdb_model_number_selection = 1, return_min_residue_distance = TRUE, show_progress = TRUE, export = FALSE, export_location = NULL, structure_file = NULL) {
  pdb_ids <- unique(dplyr::pull(data, {{ pdb_id }}))

  # create individual data retain patterns depending on which information was provided to the function.
  if (missing(chain) || all(is.na(dplyr::pull(data, {{ chain }})))) {
    data_retain_pattern <- data %>%
      dplyr::ungroup() %>%
      dplyr::mutate(retain_pattern = {{ pdb_id }}) %>%
      dplyr::pull(.data$retain_pattern) %>%
      unique()
  } else {
    if (missing(start_in_pdb) || all(is.na(dplyr::pull(data, {{ chain }})))) {
      data_retain_pattern <- data %>%
        dplyr::ungroup() %>%
        dplyr::mutate(retain_pattern = stringr::str_replace_all(paste({{ pdb_id }}, {{ chain }}, sep = "_"), pattern = "_NA", replacement = "")) %>%
        dplyr::pull(.data$retain_pattern) %>%
        unique()
    } else {
      data_retain_pattern <- data %>%
        dplyr::ungroup() %>%
        dplyr::mutate(
          start = {{ start_in_pdb }},
          end = {{ end_in_pdb }}
        ) %>% # do this so start and end position can be the same column.
        dplyr::distinct({{ pdb_id }}, {{ chain }}, .data$start, .data$end) %>%
        group_by({{ pdb_id }}, {{ chain }}, .data$start, .data$end) %>%
        dplyr::mutate(residue = ifelse(!is.na(.data$start), list(seq(.data$start, .data$end)), list(NA))) %>%
        dplyr::ungroup() %>%
        tidyr::unnest(.data$residue) %>%
        dplyr::mutate(retain_pattern = stringr::str_replace_all(paste({{ pdb_id }}, {{ chain }}, .data$residue, sep = "_"), pattern = "_NA", replacement = "")) %>%
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

    structure_file <- readr::read_tsv(structure_file, col_names = FALSE, show_col_types = FALSE, progress = FALSE)

    # load .cif file if provided
    if (file_format == ".cif") {
      structure <- structure_file %>%
        dplyr::filter(stringr::str_detect(.data$X1, pattern = "^ATOM|^HETATM")) %>%
        dplyr::mutate(X2 = stringr::str_replace_all(.data$X1, pattern = "\\s+", replacement = " ")) %>%
        tidyr::separate(.data$X2,
          sep = " ",
          into = c("x1", "atom_number", "atom_type_simple", "atom_type", "x2", "residue_name_cif", "database_chain", "entity_id", "residue_number_cif", "x3", "x", "y", "z", "site_occupancy", "b_iso_or_equivalent", "formal_charge", "residue_number_pdb", "residue_name_pdb", "auth_chain", "x4", "pdb_model_number")
        ) %>%
        dplyr::select(-c(.data$X1, .data$x1, .data$x2, .data$x3, .data$x4)) %>%
        dplyr::group_by(.data$database_chain, .data$atom_type, .data$residue_name_cif) %>%
        dplyr::mutate(residue_number_cif = ifelse(.data$residue_number_cif == ".", 1:n(), as.numeric(.data$residue_number_cif))) %>%
        dplyr::ungroup() %>%
        dplyr::mutate(
          atom_number = as.numeric(.data$atom_number),
          entity_id = as.numeric(.data$entity_id),
          residue_number_cif = as.numeric(.data$residue_number_cif),
          x = as.numeric(.data$x),
          y = as.numeric(.data$y),
          z = as.numeric(.data$z),
          site_occupancy = as.numeric(.data$site_occupancy),
          b_iso_or_equivalent = as.numeric(.data$b_iso_or_equivalent),
          residue_number_pdb = as.numeric(.data$residue_number_pdb),
          pdb_model_number = as.numeric(.data$pdb_model_number),
          pdb_id = "my_structure"
        ) %>%
        dplyr::filter(.data$pdb_model_number %in% pdb_model_number_selection) %>%
        dplyr::select(
          .data$atom_number,
          .data$x,
          .data$y,
          .data$z,
          .data$residue_name_cif,
          .data$residue_number_cif,
          .data$database_chain,
          .data$residue_name_pdb,
          .data$residue_number_pdb,
          .data$auth_chain,
          .data$pdb_id
        ) %>%
        dplyr::mutate(retain_pattern = stringr::str_replace_all(paste(.data$pdb_id, .data$auth_chain, .data$residue_number_pdb, sep = "_"), pattern = "_NA", replacement = "")) %>%
        dplyr::mutate(should_be_retained = stringr::str_detect(.data$retain_pattern, pattern = paste(paste0(data_retain_pattern, "(?=$|_)"), collapse = "|")))
    }
    # load .pdb file if provided
    if (file_format == ".pdb") {
      structure <- structure_file %>%
        dplyr::filter(stringr::str_detect(.data$X1, pattern = "^ATOM|^HETATM|^MODEL")) %>%
        dplyr::mutate(pdb_model_number = as.numeric(ifelse(stringr::str_detect(.data$X1, pattern = "^MODEL"), stringr::str_extract(.data$X1, pattern = "\\d+"), NA))) %>%
        tidyr::fill(.data$pdb_model_number, .direction = "down") %>%
        dplyr::mutate(pdb_model_number = ifelse(is.na(.data$pdb_model_number), 1, .data$pdb_model_number)) %>%
        dplyr::filter(!stringr::str_detect(.data$X1, pattern = "^MODEL")) %>%
        dplyr::mutate(
          atom_number = as.numeric(stringr::str_replace_all(stringr::str_sub(.data$X1, start = 7, end = 11), pattern = "\\s+", replacement = "")),
          residue_name_pdb = stringr::str_replace_all(stringr::str_sub(.data$X1, start = 18, end = 20), pattern = "\\s+", replacement = ""),
          auth_chain = stringr::str_replace_all(stringr::str_sub(.data$X1, start = 22, end = 22), pattern = "\\s+", replacement = ""),
          residue_number_pdb = as.numeric(stringr::str_replace_all(suppressWarnings(as.numeric(stringr::str_sub(.data$X1, start = 23, end = 26))), pattern = "\\s+", replacement = "")),
          x = as.numeric(stringr::str_replace_all(stringr::str_sub(.data$X1, start = 31, end = 38), pattern = "\\s+", replacement = "")),
          y = as.numeric(stringr::str_replace_all(stringr::str_sub(.data$X1, start = 39, end = 46), pattern = "\\s+", replacement = "")),
          z = as.numeric(stringr::str_replace_all(stringr::str_sub(.data$X1, start = 47, end = 54), pattern = "\\s+", replacement = "")),
          pdb_id = "my_structure"
        ) %>%
        dplyr::filter(.data$pdb_model_number %in% pdb_model_number_selection) %>%
        dplyr::select(-c(.data$X1, .data$pdb_model_number)) %>%
        dplyr::mutate(retain_pattern = stringr::str_replace_all(paste(.data$pdb_id, .data$auth_chain, .data$residue_number_pdb, sep = "_"), pattern = "_NA", replacement = "")) %>%
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

    if (show_progress == TRUE) {
      pb <- progress::progress_bar$new(total = length(pdb_ids), format = "Preparing structures [:bar] :current/:total (:percent) :eta")
    }

    structures <- fetch_pdb_structure(pdb_ids = pdb_ids, return_data_frame = FALSE, show_progress = show_progress) %>%
      purrr::map(.f = ~ {
        if (show_progress == TRUE) {
          pb$tick()
        }
        .x %>%
          dplyr::filter(.data$pdb_model_number %in% pdb_model_number_selection) %>%
          dplyr::select(
            .data$atom_number,
            .data$x,
            .data$y,
            .data$z,
            .data$residue_name_cif,
            .data$residue_number_cif,
            .data$database_chain,
            .data$residue_name_pdb,
            .data$residue_number_pdb,
            .data$auth_chain,
            .data$pdb_id
          ) %>%
          dplyr::mutate(retain_pattern = stringr::str_replace_all(paste(.data$pdb_id, .data$auth_chain, .data$residue_number_pdb, sep = "_"), pattern = "_NA", replacement = "")) %>%
          dplyr::mutate(should_be_retained = stringr::str_detect(.data$retain_pattern, pattern = paste(paste0(data_retain_pattern, "(?=$|_)"), collapse = "|")))
      })
  }
  # this data frame contains the subsetted structures. Specifically, the subsetted atom numbers.
  if (show_progress == TRUE) {
    pb <- progress::progress_bar$new(total = length(structures), format = "Subsetting Structures [:bar] :current/:total (:percent) :eta")
  }
  subset_structures <- structures %>%
    purrr::map(.f = ~ {
      if (show_progress == TRUE) {
        pb$tick()
      }
      .x %>%
        dplyr::filter(should_be_retained) %>%
        dplyr::distinct(.data$atom_number)
    })

  # Segments are made for each structure to prevent too long data frames when all combinations are created.
  segments <- subset_structures %>%
    purrr::map(.f = ~ {
      split(dplyr::pull(.x, .data$atom_number), ceiling(seq_along(nrow(.x)) / 1000))
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
        dplyr::distinct(.data$atom_number, .data$x, .data$y, .data$z)

      current_protein <- .y %>%
        dplyr::distinct(.data$pdb_id) %>%
        dplyr::pull(.data$pdb_id)

      if (show_progress == TRUE) {
        pb <- progress::progress_bar$new(total = length(.x), format = paste("Calculating distances for", current_protein, "[:bar] :current/:total (:percent)"))
      }

      purrr::map_dfr(
        .x = .x,
        .f = ~ {
          if (show_progress == TRUE) {
            pb$tick()
          }

          tidyr::crossing(var1 = .x, var2 = current_structure_minimum$atom_number) %>%
            dplyr::left_join(current_structure_minimum, by = c("var1" = "atom_number")) %>%
            dplyr::left_join(current_structure_minimum, by = c("var2" = "atom_number")) %>%
            dplyr::mutate(distance = sqrt((.data$x.x - .data$x.y)^2 + (.data$y.x - .data$y.y)^2 + (.data$z.x - .data$z.y)^2)) %>%
            dplyr::select(.data$var1, .data$var2, .data$distance) %>%
            dplyr::filter(.data$distance <= distance_cutoff) %>%
            dplyr::left_join(current_structure %>% dplyr::select(-.data$pdb_id), by = c("var1" = "atom_number")) %>%
            dplyr::left_join(current_structure, by = c("var2" = "atom_number"), suffix = c("_var1", "_var2"))
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
        dplyr::group_by(.data$residue_number_pdb_var1, .data$auth_chain_var1, .data$residue_number_pdb_var2, .data$auth_chain_var2) %>%
        dplyr::mutate(min_distance_residue = min(.data$distance)) %>%
        dplyr::ungroup() %>%
        dplyr::rename(
          atom_number_var1 = .data$var1,
          atom_number_var2 = .data$var2
        )

      if (return_min_residue_distance == TRUE) {
        residue_distance <- residue_distance %>%
          dplyr::select(-c(.data$atom_number_var1, .data$atom_number_var2, .data$distance)) %>%
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
