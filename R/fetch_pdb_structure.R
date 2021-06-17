#' Fetch PDB structure atom data from RCSB
#'
#' Fetches atom data for a PDB structure from RCSB. If you want to retrieve metadata about PDB structures, use the function \code{fetch_pdb()}.
#' The information retrieved is based on the .cif file of the structure, which may vary from the .pdb file.
#'
#' @param pdb_ids a character vector of PDB identifiers.
#' @param return_data_frame logical, if true, a data frame instead of a list is returned. It is recommended to only use this if not many
#' pdb structures are retrieved. Default is FALSE.
#' @param show_progress logical, if true, a progress bar will be shown. Default is TRUE.
#'
#' @return A list that contains atom data for each PDB structures provided. If return_data_frame is TRUE, a data frame with this
#' information is returned instead.
#' @import dplyr
#' @import progress
#' @import purrr
#' @import tidyr
#' @importFrom stringr str_replace_all str_detect
#' @importFrom curl has_internet
#' @importFrom magrittr %>%
#' @export
#'
#' @examples
#' \donttest{
#' head(fetch_pdb_structure(pdb_ids = c("6HG1", "1E9I", "6D3Q", "4JHW"), return_data_frame = TRUE))
#' }
fetch_pdb_structure <- function(pdb_ids, return_data_frame = FALSE, show_progress = TRUE) {
  if (!requireNamespace("httr", quietly = TRUE)) {
    stop("Package \"httr\" is needed for this function to work. Please install it.", call. = FALSE)
  }

  if (!curl::has_internet()) {
    message("No internet connection.")
    return(invisible(NULL))
  }

  batches <- purrr::map(
    .x = pdb_ids,
    .f = ~ paste0("https://files.rcsb.org/download/", .x, ".cif")
  )

  names(batches) <- pdb_ids

  if (show_progress == TRUE) {
    pb <- progress::progress_bar$new(total = length(batches))
  }

  query_result <- purrr::map2(
    .x = batches,
    .y = names(batches),
    .f = ~ {
      # query information from database
      if (!is.null(batches)) {
        query <- try_query(.x, type = "text/tab-separated-values", col_names = FALSE, quote = "")
      }
      if (show_progress == TRUE & "tbl" %in% class(query)) {
        pb$tick()
      }
      # if previous batch had a connection problem change batches to NULL, which breaks the mapping.
      if (!"tbl" %in% class(query)) {
        batches <<- NULL
      }
      # only proceed with data if it was correctly retrieved
      if ("tbl" %in% class(query)) {
        query %>%
          dplyr::filter(stringr::str_detect(X1, pattern = "^ATOM|^HETATM")) %>%
          dplyr::mutate(X2 = stringr::str_replace_all(X1, pattern = "\\s+", replacement = " ")) %>%
          tidyr::separate(X2,
            sep = " ",
            into = c("x1", "atom_number", "atom_type_simple", "atom_type", "x2", "residue_name", "chain", "entity_id", "residue_number", "x3", "x", "y", "z", "site_occupancy", "b_iso_or_equivalent", "formal_charge", "x4", "x5", "x6", "x7", "pdb_model_number")
          ) %>%
          dplyr::select(-c(.data$X1, .data$x1, .data$x2, .data$x3, .data$x4, .data$x5, .data$x6, .data$x7)) %>%
          dplyr::group_by(.data$chain, .data$atom_type, .data$residue_name) %>%
          dplyr::mutate(residue_number = ifelse(.data$residue_number == ".", 1:n(), as.numeric(.data$residue_number))) %>%
          dplyr::ungroup() %>%
          dplyr::mutate(
            atom_number = as.numeric(.data$atom_number),
            entity_id = as.numeric(.data$entity_id),
            residue_number = as.numeric(.data$residue_number),
            x = as.numeric(.data$x),
            y = as.numeric(.data$y),
            z = as.numeric(.data$z),
            site_occupancy = as.numeric(.data$site_occupancy),
            b_iso_or_equivalent = as.numeric(.data$b_iso_or_equivalent),
            pdb_model_number = as.numeric(.data$pdb_model_number),
            pdb_id = .y
          )
      }
    }
  )

  if (return_data_frame == FALSE) {
    return(query_result)
  } else {
    query_result_df <- purrr::map_dfr(
      .x = query_result,
      .f = ~.x
    )
    return(query_result_df)
  }
}
