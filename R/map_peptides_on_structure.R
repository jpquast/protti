#' Maps peptides onto a PDB structure
#'
#' Peptides are mapped onto PDB structures based on their positions. This is accomplished by replacing the B-factor information in the structure file with values that allow
#' highlighting of peptides, protein regions and amino acids when the structure is coloured by B-factor. In addition to simply highlighting peptides, protein regions or amino acids,
#' a continuous variable such as fold changes associated with them can be mapped onto the structure as a colour gradient.
#'
#' @param peptide_data a data frame containing input columns to this function. If structure files should be fetched automatically, please provide
#' column names to the following arguments: **uniprot_id**, **pdb_id**, **chain**, **start_in_pdb**, **end_in_pdb**, **map_value**. If a structure file is
#' provided in the \code{structure_file} argument, this data frame should only contain information associated with the provided structure. In case of
#' a user provided structure, column names should be provided to the following arguments: **uniprot_id**, **chain**, **start_in_pdb**, **end_in_pdb**,
#' **map_value**.
#' @param uniprot_id a character column in the \code{peptide_data} data frame that contains UniProt identifiers for a corresponding peptide, protein region or amino acid.
#' @param pdb_id a character column in the \code{peptide_data} data frame that contains PDB identifiers for structures in which a corresponding peptide, protein region or amino acid is found. This
#' column is not required if a structure file is provided in the \code{structure_file} argument.
#' @param chain a character column in the \code{peptide_data} data frame that contains the name of the chain from the PDB structure in which the peptide, protein region or amino acid is found. **Important:**
#' please provide the author defined chain definitions for both ".cif" and ".pdb" files. When the output of the \code{find_peptide_in_pdb} function
#' is used as the input for this function, this corresponds to the \code{auth_chain} column.
#' @param start_in_pdb a numeric column in the \code{peptide_data} data frame that contains start positions for peptides, protein regions or amino acids in the corresponding PDB structure. This information
#' can be obtained from the \code{find_peptide_in_pdb} function. The corresponding column in the output is called \code{peptide_start_pdb}. If amino acid positions should be
#' used, start and end positions are the same and the same column can be provided to both \code{start_in_pdb} and \code{end_in_pdb}.
#' @param end_in_pdb a numeric column in the \code{peptide_data} data frame that contains end positions for peptides, protein regions or amino acids in the corresponding PDB structure. This information
#' can be obtained from the \code{find_peptide_in_pdb} function. The corresponding column in the output is called \code{peptide_end_pdb}. If amino acid positions should be
#' used, start and end positions are the same and the same column can be provided to both \code{start_in_pdb} and \code{end_in_pdb}.
#' @param map_value a numeric column in the \code{peptide_data} data frame that contains a value associated with each peptide, protein region or amino acid. This value will be displayed as
#' a colour gradient when mapped onto the structure. The value can for example be the fold change, p-value or score associated with each peptide, protein region or amino acid (selection). If the selections should
#' be displayed with just one colour, the value in this column should be the same for every selection. For the mapping, values are scaled between 50 and 100. Regions in the structure
#' that do not map any selection receive a value of 0. If an amino acid position is associated with multiple mapped values, e.g. from different peptides, the
#' maximum mapped value will be displayed.
#' @param file_format a character vector containing the file format of the structure that will be fetched from the database for the PDB identifiers provided
#' in the \code{pdb_id} column. This can be either ".cif" or ".pdb". The default is \code{".cif"}. We recommend using ".cif" files since every structure contains a ".cif" file
#' but not every structure contains a ".pdb" file. Fetching and mapping onto ".cif" files takes longer than for ".pdb" files. If a structure file is provided in the
#' \code{structure_file} argument, the file format is detected automatically and does not need to be provided.
#' @param export_location optional, a character argument specifying the path to the location in which the fetched and altered structure files should be saved. If left empty, they will be saved in the current
#' working directory. The location should be provided in the following format "folderA/folderB".
#' @param structure_file optional, a character argument specifying the path to the location and name of a structure file in ".cif" or ".pdb" format. If a structure is provided the \code{peptide_data}
#' data frame should only contain mapping information for this structure.
#' @param show_progress a logical, if \code{show_progress = TRUE}, a progress bar will be shown (default is TRUE).
#'
#' @return The function exports a modified ".pdb" or ".cif" structure file. B-factors have been replaced with scaled (50-100) values provided in the \code{map_value} column.
#' @import dplyr
#' @import tidyr
#' @import progress
#' @importFrom purrr map2 map keep discard
#' @importFrom readr read_tsv write_tsv
#' @importFrom stringr str_sub str_detect str_extract str_extract_all str_replace
#' @importFrom magrittr %>%
#' @importFrom rlang as_name enquo .data
#' @export
#'
#' @examples
#' \dontrun{
#' map_peptides_on_structure(
#'   peptide_data = peptide_data,
#'   uniprot_id = pg_protein_accessions,
#'   pdb_id = pdb_ids,
#'   chain = auth_chain,
#'   start_in_pdb = peptide_start_pdb,
#'   end_in_pdb = peptide_end_pdb,
#'   map_value = diff,
#'   file_format = ".cif",
#'   export_location = "~/Desktop/Test"
#' )
#' }
map_peptides_on_structure <- function(peptide_data, uniprot_id, pdb_id, chain, start_in_pdb, end_in_pdb, map_value, file_format = ".cif", export_location = NULL, structure_file = NULL, show_progress = TRUE) {
  if (!file_format %in% c(".cif", ".pdb")) {
    stop('Please provide either ".cif" or ".pdb" to the file_format argument.')
  }

  if (missing(structure_file)) {
    if (!requireNamespace("httr", quietly = TRUE)) {
      stop("Package \"httr\" is needed for this function to work. Please install it.", call. = FALSE)
    }

    if (!curl::has_internet()) {
      message("No internet connection.")
      return(invisible(NULL))
    }

    # make sure export location is correct if or if not provided.
    if (missing(export_location)) {
      export_location <- ""
    } else {
      export_location <- paste0(export_location, "/")
    }

    peptide_data_filter <- peptide_data %>%
      dplyr::ungroup() %>%
      dplyr::mutate(
        start = {{ start_in_pdb }},
        end = {{ end_in_pdb }}
      ) %>% # do this so start and end position can be the same column.
      dplyr::distinct({{ uniprot_id }}, {{ pdb_id }}, {{ chain }}, .data$start, .data$end, {{ map_value }}) %>%
      tidyr::drop_na({{ pdb_id }}) %>%
      dplyr::mutate({{ map_value }} := round(scale_protti(c({{ map_value }}), method = "01") * 50 + 50, digits = 2)) %>% # Scale values between 50 and 100
      group_by({{ pdb_id }}, {{ chain }}, .data$start, .data$end) %>%
      dplyr::mutate(residue = list(seq(.data$start, .data$end))) %>%
      dplyr::ungroup() %>%
      tidyr::unnest(.data$residue) %>%
      dplyr::group_by({{ pdb_id }}, {{ chain }}, .data$residue) %>%
      dplyr::mutate({{ map_value }} := max({{ map_value }})) %>%
      dplyr::ungroup() %>%
      dplyr::distinct()

    # Extract the length of the B-factor replacement value. It is used later to ensure that the spacing stays constant in the structure file.
    length_replacement <- peptide_data_filter %>%
      dplyr::distinct({{ map_value }}) %>%
      dplyr::mutate(map_value = as.character({{ map_value }})) %>%
      dplyr::mutate(length = nchar(.data$map_value)) %>%
      dplyr::pull(.data$length) %>%
      max()

    # Extract PDB IDs for fetching
    pdb_preparation <- peptide_data_filter %>%
      dplyr::distinct({{ pdb_id }}, {{ uniprot_id }}) %>%
      dplyr::group_by({{ pdb_id }}) %>%
      dplyr::mutate(name = paste({{ uniprot_id }}, collapse = "_")) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(name = paste0({{ pdb_id }}, "_", .data$name)) %>%
      dplyr::select(-{{ uniprot_id }}) %>%
      dplyr::distinct()

    pdb_ids <- dplyr::pull(pdb_preparation, {{ pdb_id }})
    pdb_names <- dplyr::pull(pdb_preparation, .data$name)

    if (file_format == ".cif") {
      batches <- purrr::map(
        .x = pdb_ids,
        .f = ~ paste0("https://files.rcsb.org/download/", .x, ".cif")
      )

      names(batches) <- pdb_names

      if (show_progress == TRUE) {
        pb <- progress::progress_bar$new(total = length(batches), format = "  Fetching structures [:bar] :current/:total (:percent) :eta")
      }

      query_result <- purrr::map(
        .x = batches,
        .f = ~ {
          # query information from database
          if (!is.null(batches)) {
            query <- try_query(.x, type = "text/tab-separated-values", col_names = FALSE, quote = "", show_col_types = FALSE, progress = FALSE)
          }
          if (show_progress == TRUE & "data.frame" %in% class(query)) {
            pb$tick()
          }
          # if previous batch had a connection problem change batches to NULL, which breaks the mapping.
          if (!"data.frame" %in% class(query)) {
            batches <<- NULL
          }

          # only proceed with data if it was correctly retrieved
          if ("data.frame" %in% class(query)) {
            query
          }
        }
      )

      if (show_progress == TRUE) {
        pb <- progress::progress_bar$new(total = length(query_result), format = "  Replacing B-factors [:bar] :current/:total (:percent) :eta")
      }

      purrr::map2(
        .x = query_result,
        .y = names(query_result),
        .f = ~ {
          result <- .x %>%
            dplyr::mutate(atoms = ifelse(stringr::str_detect(.data$X1, pattern = "^ATOM|^HETATM"), .data$X1, NA)) %>%
            dplyr::rowwise() %>%
            dplyr::mutate(
              b_factor = stringr::str_extract_all(.data$atoms, "[\\w[:punct:]]+[:space:]+")[[1]][15], # extract b-factor values based on positions
              chain = stringr::str_extract_all(.data$atoms, "[\\w[:punct:]]+")[[1]][19],
              residue = suppressWarnings(as.numeric(stringr::str_extract_all(.data$atoms, "[\\w[:punct:]]+")[[1]][17]))
            ) %>%
            dplyr::ungroup() %>%
            dplyr::mutate(pdb_id = stringr::str_extract(.y, pattern = "^\\w+(?=_)")) %>%
            dplyr::mutate(b_factor = stringr::str_replace(.data$b_factor, pattern = "\\.", replacement = "\\\\\\.")) %>%
            dplyr::left_join(peptide_data_filter, by = c(
              "pdb_id" = rlang::as_name(rlang::enquo(pdb_id)),
              "chain" = rlang::as_name(rlang::enquo(chain)),
              "residue"
            )) %>%
            dplyr::mutate({{ map_value }} := ifelse(!is.na(.data$b_factor) & is.na({{ map_value }}), 0, {{ map_value }})) %>%
            dplyr::mutate({{ map_value }} := format(as.character({{ map_value }}), width = length_replacement + 1)) %>%
            dplyr::mutate(atoms_mod = stringr::str_replace(.data$atoms, pattern = paste0("(?<= )", .data$b_factor), replacement = {{ map_value }})) %>%
            dplyr::mutate(X1 = ifelse(!is.na(.data$atoms_mod), .data$atoms_mod, .data$X1)) %>%
            dplyr::select(.data$X1) %>%
            readr::write_tsv(file = paste0(export_location, .y, file_format), quote = "none", escape = "none", col_names = FALSE, progress = FALSE)

          if (show_progress == TRUE) {
            pb$tick()
          }

          result
        }
      )

      return(invisible(NULL))
    }
    if (file_format == ".pdb") {
      # same as above but for .pdb files. Indexing is different.
      batches <- purrr::map(
        .x = pdb_ids,
        .f = ~ paste0("https://files.rcsb.org/download/", .x, ".pdb")
      )

      names(batches) <- pdb_names

      if (show_progress == TRUE) {
        pb <- progress::progress_bar$new(total = length(batches), format = "  Fetching structures [:bar] :current/:total (:percent) :eta")
      }

      query_result <- purrr::map(
        .x = batches,
        .f = ~ {
          # query information from database
          if (!is.null(batches)) {
            query <- suppressMessages(try_query(.x, type = "text/tab-separated-values", col_names = FALSE, quote = "", show_col_types = FALSE, progress = FALSE))
          }
          if (show_progress == TRUE & "data.frame" %in% class(query)) {
            pb$tick()
          }

          # only proceed with data if it was correctly retrieved
          if ("data.frame" %in% class(query)) {
            query
          }
        }
      )

      not_fetched <- query_result %>%
        purrr::keep(is.null) %>%
        names()

      if (length(not_fetched) != 0) {
        message('The following structures were not fetched, likely because no ".pdb" file is available. Try using the ".cif" format for these. ', paste(not_fetched, collapse = ", "))
      }

      query_result <- query_result %>%
        purrr::discard(is.null)

      if (show_progress == TRUE) {
        pb <- progress::progress_bar$new(total = length(query_result), format = "  Replacing B-factors [:bar] :current/:total (:percent) :eta")
      }

      result <- purrr::map2(
        .x = query_result,
        .y = names(query_result),
        .f = ~ {
          result <- .x %>%
            dplyr::mutate(atoms = ifelse(stringr::str_detect(.data$X1, pattern = "^ATOM|^HETATM"), .data$X1, NA)) %>%
            dplyr::mutate(
              chain = stringr::str_sub(.data$atoms, start = 22, end = 22),
              residue = suppressWarnings(as.numeric(stringr::str_sub(.data$atoms, start = 23, end = 26)))
            ) %>%
            dplyr::mutate(pdb_id = stringr::str_extract(.y, pattern = "^\\w+(?=_)")) %>%
            dplyr::left_join(peptide_data_filter, by = c(
              "pdb_id" = rlang::as_name(rlang::enquo(pdb_id)),
              "chain" = rlang::as_name(rlang::enquo(chain)),
              "residue"
            )) %>%
            dplyr::mutate({{ map_value }} := ifelse(!is.na(.data$residue) & is.na({{ map_value }}), 0, {{ map_value }})) %>%
            dplyr::mutate({{ map_value }} := format(as.character({{ map_value }}), width = 6, justify = "right")) %>%
            dplyr::mutate({{ map_value }} := ifelse(str_detect({{ map_value }}, pattern = "NA"), NA, {{ map_value }})) %>%
            dplyr::mutate(atoms_mod = `str_sub<-`(.data$atoms, 61, 66, value = {{ map_value }})) %>%
            dplyr::mutate(X1 = ifelse(!is.na(.data$atoms_mod), .data$atoms_mod, .data$X1)) %>%
            dplyr::select(.data$X1) %>%
            readr::write_tsv(file = paste0(export_location, .y, file_format), quote = "none", escape = "none", col_names = FALSE, progress = FALSE)

          if (show_progress == TRUE) {
            pb$tick()
          }

          result
        }
      )

      return(invisible(NULL))
    }
  }

  if (!missing(structure_file)) {
    file_format <- stringr::str_sub(structure_file, start = -4, end = -1)

    if (!file_format %in% c(".cif", ".pdb")) {
      stop('Please either provide a ".cif" or ".pdb" structure file.')
    }

    file_name <- stringr::str_extract(structure_file, pattern = "[^/]+[:punct:]\\w+$")

    pdb_file <- readr::read_tsv(structure_file, col_names = FALSE, show_col_types = FALSE, progress = FALSE)

    peptide_data_filter <- peptide_data %>%
      dplyr::ungroup() %>%
      dplyr::mutate(
        start = {{ start_in_pdb }},
        end = {{ end_in_pdb }}
      ) %>%
      dplyr::distinct({{ uniprot_id }}, {{ chain }}, .data$start, .data$end, {{ map_value }}) %>%
      dplyr::mutate({{ map_value }} := round(scale_protti(c({{ map_value }}), method = "01") * 50 + 50, digits = 2)) %>% # Scale values between 50 and 100
      group_by({{ chain }}, .data$start, .data$end) %>%
      dplyr::mutate(residue = list(seq(.data$start, .data$end))) %>%
      dplyr::ungroup() %>%
      tidyr::unnest(.data$residue) %>%
      dplyr::group_by({{ chain }}, .data$residue) %>%
      dplyr::mutate({{ map_value }} := max({{ map_value }})) %>%
      dplyr::ungroup() %>%
      dplyr::distinct()

    length_replacement <- peptide_data_filter %>%
      dplyr::distinct({{ map_value }}) %>%
      dplyr::mutate(map_value = as.character({{ map_value }})) %>%
      dplyr::mutate(length = nchar(.data$map_value)) %>%
      dplyr::pull(.data$length) %>%
      max()

    if (file_format == ".cif") {
      pdb_file %>%
        dplyr::mutate(atoms = ifelse(stringr::str_detect(.data$X1, pattern = "^ATOM|^HETATM"), .data$X1, NA)) %>%
        dplyr::rowwise() %>%
        dplyr::mutate(
          b_factor = stringr::str_extract_all(.data$atoms, "[\\w[:punct:]]+[:space:]+")[[1]][15], # extract b-factor values based on positions
          chain = stringr::str_extract_all(.data$atoms, "[\\w[:punct:]]+")[[1]][19],
          residue = suppressWarnings(as.numeric(stringr::str_extract_all(.data$atoms, "[\\w[:punct:]]+")[[1]][17]))
        ) %>%
        dplyr::ungroup() %>%
        dplyr::mutate(b_factor = stringr::str_replace(.data$b_factor, pattern = "\\.", replacement = "\\\\\\.")) %>%
        dplyr::left_join(peptide_data_filter, by = c(
          "chain" = rlang::as_name(rlang::enquo(chain)),
          "residue"
        )) %>%
        dplyr::mutate({{ map_value }} := ifelse(!is.na(.data$b_factor) & is.na({{ map_value }}), 0, {{ map_value }})) %>%
        dplyr::mutate({{ map_value }} := format(as.character({{ map_value }}), width = length_replacement + 1)) %>%
        dplyr::mutate(atoms_mod = stringr::str_replace(.data$atoms, pattern = paste0("(?<= )", .data$b_factor), replacement = {{ map_value }})) %>%
        dplyr::mutate(X1 = ifelse(!is.na(.data$atoms_mod), .data$atoms_mod, .data$X1)) %>%
        dplyr::select(.data$X1) %>%
        readr::write_tsv(file = paste0(export_location, "modified_", file_name), quote = "none", escape = "none", col_names = FALSE, progress = FALSE)

      return(invisible(NULL))
    }

    if (file_format == ".pdb") {
      pdb_file %>%
        dplyr::mutate(atoms = ifelse(stringr::str_detect(.data$X1, pattern = "^ATOM|^HETATM"), .data$X1, NA)) %>%
        dplyr::mutate(
          chain = stringr::str_sub(.data$atoms, start = 22, end = 22),
          residue = suppressWarnings(as.numeric(stringr::str_sub(.data$atoms, start = 23, end = 26)))
        ) %>%
        dplyr::left_join(peptide_data_filter, by = c(
          "chain" = rlang::as_name(rlang::enquo(chain)),
          "residue"
        )) %>%
        dplyr::mutate({{ map_value }} := ifelse(!is.na(.data$residue) & is.na({{ map_value }}), 0, {{ map_value }})) %>%
        dplyr::mutate({{ map_value }} := format(as.character({{ map_value }}), width = 6, justify = "right")) %>%
        dplyr::mutate({{ map_value }} := ifelse(str_detect({{ map_value }}, pattern = "NA"), NA, {{ map_value }})) %>%
        dplyr::mutate(atoms_mod = `str_sub<-`(.data$atoms, 61, 66, value = {{ map_value }})) %>%
        dplyr::mutate(X1 = ifelse(!is.na(.data$atoms_mod), .data$atoms_mod, .data$X1)) %>%
        dplyr::select(.data$X1) %>%
        readr::write_tsv(file = paste0(export_location, "modified_", file_name), quote = "none", escape = "none", col_names = FALSE, progress = FALSE)

      return(invisible(NULL))
    }
  }
}