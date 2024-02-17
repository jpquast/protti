#' Fetch protein disorder and mobility information from MobiDB
#'
#' Fetches information about disordered and flexible protein regions from MobiDB.
#'
#' @param uniprot_ids optional, a character vector of UniProt identifiers for which information
#' should be fetched. This argument is mutually exclusive to the \code{organism_id} argument.
#' @param organism_id optional, a character value providing the NCBI taxonomy identifier of an organism
#' (TaxId) of an organism for which all available information should be retreived. This
#' argument is mutually exclusive to the \code{uniprot_ids} argument.
#' @param show_progress a logical value; if `TRUE` a progress bar will be shown.
#' Default is `TRUE`.
#'
#' @return A data frame that contains start and end positions for disordered and flexible protein
#' regions. The \code{feature} column contains information on the source of this
#' annotation. More information on the source can be found
#' \href{https://mobidb.bio.unipd.it/about/mobidb}{here}.
#' @import progress
#' @importFrom rlang .data
#' @importFrom purrr map_dfr keep
#' @importFrom dplyr filter mutate
#' @importFrom tidyr unnest separate
#' @importFrom stringr str_split
#' @importFrom methods is
#' @importFrom magrittr %>%
#' @export
#'
#' @examples
#' \donttest{
#' fetch_mobidb(
#'   uniprot_ids = c("P0A799", "P62707")
#' )
#' }
fetch_mobidb <- function(uniprot_ids = NULL, organism_id = NULL, show_progress = TRUE) {
  if (!curl::has_internet()) {
    message("No internet connection.")
    return(invisible(NULL))
  }
  if (!missing(uniprot_ids) & !missing(organism_id)) {
    stop(strwrap("Please only provide either a list of UniProt identifiers or one organism ID!",
      prefix = "\n", initial = ""
    ))
  }

  base_url <- "https://mobidb.bio.unipd.it/api/download?"

  if (!missing(uniprot_ids)) {
    # Check uniprot ID validity
    uniprot_ids <- stats::na.omit(uniprot_ids)
    id_test <- stringr::str_detect(uniprot_ids,
      pattern = "^([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})$"
    )
    non_conform_ids <- uniprot_ids[!id_test]
    uniprot_ids <- uniprot_ids[id_test]

    if (length(non_conform_ids) != 0) {
      warning(strwrap("These UniProt accession numbers did not conform
to uniprot standards and were skipped from fetching: ",
        prefix = "\n", initial = ""
      ), paste(non_conform_ids, collapse = ", "))
    }
    if (length(uniprot_ids) == 0) {
      stop("No valid UniProt accession numbers found.")
    }

    if (length(uniprot_ids) < 800) {
      # generate url
      url <- paste0(
        base_url,
        "acc=",
        paste0(uniprot_ids, collapse = ","),
        "&format=tsv"
      )

      url_list <- list(url)
    } else {
      uniprot_id_list <- split(uniprot_ids, ceiling(seq_along(uniprot_ids) / 800))
      url_list <- purrr::map(
        .x = uniprot_id_list,
        .f = ~ {
          paste0(
            base_url,
            "acc=",
            paste0(.x, collapse = ","),
            "&format=tsv"
          )
        }
      )
    }
  }

  if (!missing(organism_id)) {
    # generate url
    url <- paste0(
      base_url,
      "ncbi_taxon_id=",
      organism_id,
      "&format=tsv"
    )

    url_list <- list(url)
  }

  if (show_progress == TRUE) {
    pb <- progress::progress_bar$new(
      total = length(url_list),
      format = "Fetching information from MobiDB [:bar] :current/:total (:percent) :eta"
    )
  }

  query_result <- purrr::map(
    .x = url_list,
    .f = ~ {
      # query information from database
      query <- try_query(
        url = .x,
        timeout = 60
      )

      if (show_progress == TRUE) {
        pb$tick()
      }

      query
    }
  )

  # catch any IDs that have not been fetched correctly
  error_list <- query_result %>%
    purrr::keep(.p = ~ is.character(.x))

  if (length(error_list) != 0) {
    error_table <- tibble::tibble(
      id = paste0("IDs: ", as.numeric(names(error_list)) * 800 - 800 + 1, " to ", as.numeric(names(error_list)) * 800),
      error = unlist(error_list)
    ) %>%
      dplyr::distinct()

    message("The following IDs have not been retrieved correctly.")
    message(paste0(utils::capture.output(error_table), collapse = "\n"))
  }

  # only keep data in output and transform to data.frame
  query_result <- query_result %>%
    purrr::keep(.p = ~ !is.character(.x)) %>%
    purrr::map_dfr(.f = ~.x)

  if (length(query_result) == 0) {
    message("No valid information was retrieved!")
    return(invisible(NULL))
  }

  result <- query_result %>%
    dplyr::mutate(start..end = stringr::str_split(.data$start..end, pattern = ",")) %>%
    tidyr::unnest(.data$start..end) %>%
    tidyr::separate(.data$start..end, into = c("start", "end"), sep = "\\.\\.") %>%
    dplyr::select(-.data$length) %>%
    dplyr::rename(accession = .data$acc)

  result
}
