#' Fetch protein data from UniProt
#'
#' Fetches protein metadata from UniProt.
#'
#' @param uniprot_ids A character vector of UniProt accession numbers
#' @param columns Metadata columns that should be imported from UniProt (all possible columns can be found \href{https://www.uniprot.org/help/uniprotkb_column_names}{here}.)
#' @param batchsize Size of batch of proteins for a single query
#' @param show_progress Logical, if true, a progress bar will be shown
#'
#' @return A data frame that contains all protein metadata specified in \code{columns} for the proteins provided. If an invalid ID
#' was provided that contains a valid UniProt ID, the valid portion of the ID is fetched and the invalid input ID is saved in a
#' column called \code{input_id}.
#' @import dplyr
#' @importFrom janitor make_clean_names
#' @import progress
#' @import purrr
#' @importFrom tidyr drop_na
#' @importFrom stringr str_detect str_extract_all
#' @importFrom tibble tibble
#' @importFrom curl has_internet
#' @export
#'
#' @examples
#' \donttest{
#' fetch_uniprot(c("P36578", "O43324", "Q00796"))
#' }
fetch_uniprot <-
  function(uniprot_ids,
           columns = c(
             "protein names",
             "length",
             "sequence",
             "genes",
             "database(GeneID)",
             "database(String)",
             "go(molecular function)",
             "go(biological process)",
             "go(cellular compartment)",
             "interactor",
             "feature(ACTIVE SITE)",
             "feature(BINDING SITE)",
             "feature(METAL BINDING)",
             "chebi(Cofactor)",
             "chebi(Catalytic activity)",
             "database(PDB)"
           ),
           batchsize = 200,
           show_progress = TRUE) {
    if (!requireNamespace("httr", quietly = TRUE)) {
      stop("Package \"httr\" is needed for this function to work. Please install it.", call. = FALSE)
    }

    if (!curl::has_internet()) {
      message("No internet connection.")
      return(invisible(NULL))
    }
    . <- NULL
    columns <- c("id", columns)
    column_names <- janitor::make_clean_names(columns)
    collapsed_columns <- paste(columns, collapse = ",")
    uniprot_ids <- stats::na.omit(uniprot_ids)
    id_test <- stringr::str_detect(uniprot_ids, pattern = "^([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})$")
    uniprot_ids_filtered <- uniprot_ids[id_test]
    non_conform_ids <- uniprot_ids[!id_test]
    # if non_conform_ids contain IDs they are extracted and fetched.
    contains_valid_id <- non_conform_ids[stringr::str_detect(non_conform_ids, pattern = "[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}")]
    valid_id_annotations <- tibble::tibble(input_id = contains_valid_id) %>%
      dplyr::mutate(id = stringr::str_extract_all(.data$input_id, pattern = "[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}")) %>%
      tidyr::unnest(.data$id)
    uniprot_ids_filtered <- c(uniprot_ids_filtered, valid_id_annotations$id)
    non_identifiable_id <- non_conform_ids[!stringr::str_detect(non_conform_ids, pattern = "[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}")]

    if (length(non_identifiable_id) != 0) {
      warning("These UniProt accession numbers did not conform to uniprot standards and were skipped from importing: ", paste(non_identifiable_id, collapse = ", "))
    }
    if (length(contains_valid_id) != 0) {
      warning('The following input IDs were found to contain valid uniprot accession numbers. \nThey were fetched and the original input ID can be found in the "input_id" column: ', paste(contains_valid_id, collapse = ", "))
    }
    if (length(uniprot_ids_filtered) == 0) {
      stop("No valid UniProt accession numbers found.")
    }

    url <- "http://www.uniprot.org/uniprot/?query="
    batches <- split(uniprot_ids_filtered, ceiling(seq_along(uniprot_ids_filtered) / batchsize))
    if (show_progress == TRUE) {
      pb <- progress::progress_bar$new(total = length(batches))
    }
    # fetch all batches
    result <- purrr::map_df(batches, function(x) {
      id_query <- paste(paste0("id:", x), collapse = "+or+")
      query_url <- utils::URLencode(paste0(
        url, id_query, "&format=tab&columns=",
        collapsed_columns
      ))
      # only try to fetch more batches if previous cycle did not encounter a connection problem.
      if (!is.null(batches)) {
        query <- try_query(query_url)
      }
      if (show_progress == TRUE & "tbl" %in% class(query)) {
        pb$tick()
      }
      # if previous batch had a connection problem change batches to NULL, which breaks the mapping.
      if (!"tbl" %in% class(query)) {
        batches <<- NULL
      }
      query
    })
    if (length(result) == 0) {
      return(invisible(NULL))
    }
    colnames(result) <- column_names

    # rescue cases in which the used ID is an old ID. The information is retrieved for the new ID and then parsed to the old ID.

    new <- result %>%
      dplyr::mutate(new = ifelse(stringr::str_detect(.[[2]], pattern = "Merged"), stringr::str_extract_all(.[[2]], pattern = "(?<=into )[A-Z0-9]+", simplify = TRUE), NA)) %>%
      dplyr::distinct(id, new) %>%
      tidyr::drop_na(new)

    new_ids <- new$new

    if (length(new_ids) == 0) {
      if (length(contains_valid_id) != 0) {
        result <- result %>%
          dplyr::left_join(valid_id_annotations, by = "id") %>%
          dplyr::relocate(.data$id, .data$input_id)
      }
      return(result)
    }
    new_id_query <- paste(paste0("id:", new_ids), collapse = "+or+")
    new_query_url <- utils::URLencode(paste0(
      url, new_id_query, "&format=tab&columns=",
      collapsed_columns
    ))

    new_result <- try_query(new_query_url)
    # Ff a problem occurs at this step NULL is returned.
    if (is.null(new_result)) {
      return(new_result)
    }
    colnames(new_result) <- column_names

    new_result <- new %>%
      dplyr::left_join(new_result, by = c("new" = "id")) %>%
      dplyr::select(-new)

    result <- result %>%
      dplyr::filter(!stringr::str_detect(.[[2]], pattern = "Merged")) %>%
      dplyr::bind_rows(new_result)

    if (length(contains_valid_id) != 0) {
      result <- result %>%
        dplyr::left_join(valid_id_annotations, by = "id") %>%
        dplyr::relocate(.data$id, .data$input_id)
    }

    result
  }
