#' Fetch protein data from UniProt
#'
#' Fetches protein metadata from UniProt.
#'
#' @param uniprot_ids a character vector of UniProt accession numbers.
#' @param columns a character vector of metadata columns that should be imported from UniProt (all
#' possible columns can be found \href{https://www.uniprot.org/help/return_fields}{here}. For 
#' cross-referenced database provide the database name with the prefix "xref_", e.g. \code{"xref_pdb"})
#' @param batchsize a numeric value that specifies the number of proteins processed in a single
#' single query. Default and max value is 200. 
#' @param show_progress a logical value that determines if a progress bar will be shown. Default
#' is TRUE.
#'
#' @return A data frame that contains all protein metadata specified in \code{columns} for the
#' proteins provided. If an invalid ID was provided that contains a valid UniProt ID, the valid
#' portion of the ID is fetched and the invalid input ID is saved in a column called \code{input_id}.
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
             "protein_name",
             "length",
             "sequence",
             "gene_names",
             "xref_geneid",
             "xref_string",
             "go_f",
             "go_p",
             "go_c",
             "cc_interaction",
             "ft_act_site",
             "ft_binding",
             "cc_cofactor",
             "cc_catalytic_activity",
             "xref_pdb"
           ),
           batchsize = 200,
           show_progress = TRUE) {
    if (!curl::has_internet()) {
      message("No internet connection.")
      return(invisible(NULL))
    }
    . <- NULL
    if(show_progress){
      message("Please be aware that some column names have changed due to UniProt updating its API! This might cause Errors in your code. You can fix it be replacing the old column names with new one.")
    }
    if (batchsize > 500){
      stop("Please provide a batchsize that is smaller or equal to 500!")
    }
    columns <- c("accession", columns)
    column_names <- janitor::make_clean_names(columns)
    collapsed_columns <- paste(columns, collapse = ",")
    uniprot_ids <- stats::na.omit(uniprot_ids)
    id_test <- stringr::str_detect(uniprot_ids,
      pattern = "^([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})$"
    )
    uniprot_ids_filtered <- uniprot_ids[id_test]
    non_conform_ids <- uniprot_ids[!id_test]
    # if non_conform_ids contain IDs they are extracted and fetched.
    contains_valid_id <- non_conform_ids[stringr::str_detect(non_conform_ids,
      pattern = "[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}"
    )]
    valid_id_annotations <- tibble::tibble(input_id = contains_valid_id) %>%
      dplyr::mutate(accession = stringr::str_extract_all(.data$input_id,
        pattern = "[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}"
      )) %>%
      tidyr::unnest(.data$accession)
    uniprot_ids_filtered <- c(uniprot_ids_filtered, valid_id_annotations$accession)
    non_identifiable_id <- non_conform_ids[!stringr::str_detect(non_conform_ids,
      pattern = "[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}"
    )]

    if (length(non_identifiable_id) != 0) {
      warning(strwrap("These UniProt accession numbers did not conform
to uniprot standards and were skipped from importing: ",
        prefix = "\n", initial = ""
      ), paste(non_identifiable_id, collapse = ", "))
    }
    if (length(contains_valid_id) != 0) {
      warning(strwrap('The following input IDs were found to contain valid uniprot accession numbers.
They were fetched and the original input ID can be found in the "input_id" column: ',
        prefix = "\n", initial = ""
      ), paste(contains_valid_id, collapse = ", "))
    }
    if (length(uniprot_ids_filtered) == 0) {
      stop("No valid UniProt accession numbers found.")
    }

    url <- "http://rest.uniprot.org/uniprotkb/stream?query="
    # This function splits IDs into batches that are retrieved batch by batch
    # The new UniProt API allows this functionality on the server site by using
    # the size and cursor parameters. We don't change this since the code below also works.
    batches <- split(uniprot_ids_filtered, ceiling(seq_along(uniprot_ids_filtered) / batchsize))
    if (show_progress == TRUE) {
      pb <- progress::progress_bar$new(total = length(batches))
    }
    # fetch all batches
    result <- purrr::map(batches, function(x) {
      id_query <- paste(paste0("accession_id:", x), collapse = "+OR+")
      query_url <- utils::URLencode(paste0(
        url, id_query, "&format=tsv&fields=",
        collapsed_columns
      ))

      query <- try_query(query_url, progress = FALSE, show_col_types = FALSE)

      if (show_progress == TRUE) {
        pb$tick()
      }

      query
    })

    # catch any IDs that have not been fetched correctly
    error_list <- result %>%
      purrr::keep(.p = ~ is.character(.x))

    if (length(error_list) != 0) {
      error_table <- tibble::tibble(
        id = paste0("IDs: ", as.numeric(names(error_list)) * batchsize - batchsize + 1, " to ", as.numeric(names(error_list)) * batchsize),
        error = unlist(error_list)
      ) %>%
        dplyr::distinct()

      message("The following IDs have not been retrieved correctly.")
      message(paste0(utils::capture.output(error_table), collapse = "\n"))
    }

    # only keep data in output and transform to data.frame
    result <- result %>%
      purrr::keep(.p = ~ !is.character(.x))

    if (length(result) == 0) {
      message("No valid information was retrieved!")
      return(invisible(NULL))
    }

    result <- result %>%
      purrr::map_dfr(
        .f = ~.x
      )

    colnames(result) <- column_names

    # rescue cases in which the used ID is an old ID. The information is
    # retrieved for the new ID and then parsed to the old ID.

    new <- result %>%
      dplyr::mutate(new = ifelse(stringr::str_detect(.[[2]], pattern = "Merged"),
        stringr::str_extract_all(.[[2]], pattern = "(?<=into )[A-Z0-9]+", simplify = TRUE),
        NA
      )) %>%
      dplyr::distinct(.data$accession, .data$new) %>%
      tidyr::drop_na(.data$new)

    new_ids <- new$new

    if (length(new_ids) == 0) {
      if (length(contains_valid_id) != 0) {
        result <- result %>%
          dplyr::left_join(valid_id_annotations, by = "accession") %>%
          dplyr::relocate(.data$accession, .data$input_id)
      }
      return(result)
    }
    new_id_query <- paste(paste0("accession_id:", new_ids), collapse = "+OR+")
    new_query_url <- utils::URLencode(paste0(
      url, new_id_query, "&format=tsv&fields=",
      collapsed_columns
    ))

    new_result <- try_query(new_query_url, progress = FALSE, show_col_types = FALSE)
    # If a problem occurs at this step NULL is returned.
    if (!methods::is(new_result, "data.frame")) {
      message(new_result)
      return(invisible(NULL))
    }
    colnames(new_result) <- column_names

    new_result <- new %>%
      dplyr::left_join(new_result, by = c("new" = "accession")) %>%
      dplyr::select(-new)

    result <- result %>%
      dplyr::filter(!stringr::str_detect(.[[2]], pattern = "Merged")) %>%
      dplyr::bind_rows(new_result)

    if (length(contains_valid_id) != 0) {
      result <- result %>%
        dplyr::left_join(valid_id_annotations, by = "accession") %>%
        dplyr::relocate(.data$accession, .data$input_id)
    }

    result
  }
