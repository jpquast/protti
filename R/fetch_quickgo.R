#' Fetch information from the QuickGO API
#'
#' Fetches gene ontology (GO) annotations, terms or slims from the QuickGO EBI database.
#' Annotations can be retrieved for specific UniProt IDs or NCBI taxonomy identifiers. When
#' terms are retrieved, a complete list of all GO terms is returned. For the generation of
#' a slim dataset you can provide GO IDs that should be considered. A slim dataset is a subset
#' GO dataset that considers all child terms of the supplied IDs.
#'
#' @param type a character value that indicates if gene ontology terms, annotations or slims
#' should be retrieved. The possible values therefore include "annotations", "terms" and "slims".
#' If annotations are retrieved, the maximum number of results is 2,000,000.
#' @param id_annotations an optional character vector that specifies UniProt IDs for which GO annotations
#' should be retrieved. This argument should only be provided if annotations are retrieved.
#' @param taxon_id_annotations an optional character value that specifies the NCBI taxonomy identifier (TaxId)
#' for an organism for which GO annotations should be retrieved.
#' This argument should only be provided if annotations are retrieved.
#' @param ontology_annotations an optional character value that specifies the ontology that should be retrieved.
#' This can either have the values "all", "molecular_function", "biological_process" or
#' "cellular_component". This argument should only be provided if annotations are retrieved.
#' @param go_id_slims an optional character vector that specifies gene ontology IDs (e.g. GO:0046872) for which
#' a slim go set should be generated. This argument should only be provided if slims are retrieved.
#' @param relations_slims an optional character vector that specifies the relations of GO IDs that should be
#' considered for the generation of the slim dataset. This argument should only be provided if slims are retrieved.
#' @param timeout a numeric value specifying the time in seconds until the download times out.
#' The default is 1200 seconds.
#' @param max_tries a numeric value that specifies the number of times the function tries to download
#' the data in case an error occurs. The default is 2.
#' @param show_progress a logical value that indicates if a progress bar will be shown.
#' Default is TRUE.
#'
#' @return A data frame that contains descriptive information about gene ontology annotations, terms or slims
#' depending on what the input "type" was.
#' @import dplyr
#' @import progress
#' @import purrr
#' @importFrom methods is
#' @importFrom tidyr unnest
#' @importFrom janitor clean_names
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @importFrom curl has_internet
#' @export
#'
#' @examples
#' \donttest{
#' # Annotations
#' annotations <- fetch_quickgo(
#'   type = "annotations",
#'   id = c("P63328", "Q4FFP4"),
#'   ontology = "molecular_function"
#' )
#'
#' head(annotations)
#'
#' # Terms
#' terms <- fetch_quickgo(type = "terms")
#'
#' head(terms)
#'
#' # Slims
#' slims <- fetch_quickgo(
#'   type = "slims",
#'   go_id_slims = c("GO:0046872", "GO:0051540")
#' )
#'
#' head(slims)
#' }
fetch_quickgo <- function(type = "annotations",
                          id_annotations = NULL,
                          taxon_id_annotations = NULL,
                          ontology_annotations = "all",
                          go_id_slims = NULL,
                          relations_slims = c("is_a", "part_of", "regulates", "occurs_in"),
                          timeout = 1200,
                          max_tries = 2,
                          show_progress = TRUE) {
  type <- match.arg(type, c("annotations", "terms", "slims"))

  # Check if there is an internet connection

  if (!curl::has_internet()) {
    message("No internet connection.")
    return(invisible(NULL))
  }

  # Retrieve annotation information
  if (type == "annotations") {
    if (!missing(go_id_slims)) {
      warning(strwrap('You provided IDs to "go_id_slims" but chose to retrieve "annotations" and not "slims".
              Therefore, IDs will not be considered.', prefix = "\n", initial = ""))
    }

    base_url <- "https://www.ebi.ac.uk/QuickGO/services/annotation/downloadSearch?"

    # ontology
    ontology_annotations <- match.arg(ontology_annotations, c("all", "molecular_function", "biological_process", "cellular_component"))
    if (ontology_annotations == "all") {
      aspect <- ""
    } else {
      aspect <- paste0("aspect=", ontology_annotations, "&")
    }
    # ids
    if (missing(id_annotations)) {
      geneProductID <- ""
    } else {
      collapsed_ids <- paste0(id_annotations, collapse = ",")
      geneProductID <- paste0("geneProductId=", collapsed_ids, "&")
    }
    # taxon id
    if (missing(taxon_id_annotations)) {
      taxonId <- ""
    } else {
      taxonId <- paste0("taxonId=", taxon_id_annotations, "&")
    }
    # download limit
    downlaod_limit <- "downloadLimit=2000000&"
    # fields
    fields <- "includeFields=goName&selectedFields=geneProductId,goName,symbol,qualifier,goId,goAspect,evidenceCode,goEvidence,reference,withFrom,taxonId,assignedBy,extensions,date"

    # final url
    url <- paste0(base_url, aspect, geneProductID, taxonId, downlaod_limit, fields)

    if (show_progress == TRUE) {
      message("Retrieving GO annotations ... ", appendLF = FALSE)
      start_time <- Sys.time()
    }

    query_result <- try_query(url, timeout = timeout, max_tries = max_tries, accept = "text/tsv")

    if (show_progress == TRUE) {
      message("DONE", paste0("(", round(as.numeric(difftime(Sys.time(), start_time, units = "secs")), digits = 2), "s)"))
    }

    if (methods::is(query_result, "character")) {
      message(query_result)
      if (stringr::str_detect(query_result, pattern = "Timeout")){
        message('Consider increasing the "timeout" or "max_tries" argument. \n')
      }
      return(invisible(NULL))
    }

    if (ontology_annotations != "all") {
      query_result_clean <- query_result %>%
        janitor::clean_names() %>%
        dplyr::mutate(go_aspect = ontology_annotations)
    } else {
      query_result_clean <- query_result %>%
        janitor::clean_names()
    }


    return(query_result_clean)
  }

  # Retrieve terms
  if (type == "terms") {
    if (!missing(go_id_slims)) {
      warning(strwrap('You provided IDs to "go_id_slims" but chose to retrieve "terms" and not "slims".
              Therefore, IDs will not be considered.', prefix = "\n", initial = ""))
    }
    if (!missing(id_annotations)) {
      warning(strwrap('You provided IDs to "id_annotations" but chose to retrieve "terms" and not "annotations".
              Therefore, IDs will not be considered.', prefix = "\n", initial = ""))
    }
    if (!missing(taxon_id_annotations)) {
      warning(strwrap('You provided an ID to "taxon_id_annotations" but chose to retrieve "terms" and not "annotations".
              Therefore, the ID will not be considered.', prefix = "\n", initial = ""))
    }

    base_url <- "https://www.ebi.ac.uk/QuickGO/services/ontology/go/terms?page="

    # first we request one batch of test data to find the maximum number of pages

    page <- "1"

    url <- paste0(base_url, page)

    test_query <- try_query(url,
      type = "application/json",
      timeout = timeout,
      max_tries = max_tries
    )

    if (methods::is(test_query, "character")) {
      message(test_query, "\n")
      if (stringr::str_detect(test_query, pattern = "Timeout")){
        message('Consider increasing the "timeout" or "max_tries" argument. \n')
      }
      return(invisible(NULL))
    }

    pages <- 1:test_query[["pageInfo"]][["total"]]

    if (show_progress == TRUE) {
      pb <- progress::progress_bar$new(
        total = length(pages),
        format = "  Fetching GO terms [:bar] (:percent) :eta",
        show_after = 0
      )
      pb$tick(0)
    }

    query_result_list <- purrr::map(
      .x = pages,
      .f = function(x) {
        query_url <- paste0(
          base_url,
          x
        )
        query <- try_query(query_url,
          type = "application/json",
          timeout = timeout,
          max_tries = max_tries,
          simplifyDataFrame = TRUE
        )

        if (show_progress == TRUE) {
          pb$tick()
        }
        query
      }
    )

    # if any page was not fetched correctly return NULL and the error message

    error_vector <- query_result_list %>%
      purrr::keep(.p = ~ is.character(.x)) %>%
      unlist()

    if (!is.null(error_vector)) {
      message(paste(unique(error_vector), collapse = ", "))

      if (any(stringr::str_detect(unique(error_vector), pattern = "Timeout"))){
        message('Consider increasing the "timeout" or "max_tries" argument. \n')
      }

      return(invisible(NULL))
    } else {
      query_result <- query_result_list %>%
        purrr::map_dfr(
          .f = ~ .x[["results"]]
        )
    }

    pre_extracted_result <- query_result %>%
      janitor::clean_names() %>%
      tidyr::unnest("definition") %>%
      dplyr::rename(main_name = "name") %>%
      dplyr::rename(definition = "text") %>%
      select(-c("xrefs")) %>%
      dplyr::rename(
        main_id = "id",
        ontology = "aspect"
      )

    # synonyms (not used yet but could be added in the future)

    synonyms <- pre_extracted_result %>%
      dplyr::distinct(.data$main_id, .data$synonyms) %>%
      tidyr::unnest("synonyms") %>%
      dplyr::rename(
        synonym = "name",
        synonym_type = "type"
      )

    # children

    children <- pre_extracted_result %>%
      dplyr::distinct(.data$main_id, .data$children) %>%
      tidyr::unnest("children") %>%
      dplyr::rename(
        child_id = "id",
        children_relation = "relation"
      ) %>%
      dplyr::group_by(.data$main_id) %>%
      dplyr::mutate(
        child_id = paste0(.data$child_id, collapse = ";"),
        children_relation = paste0(.data$children_relation, collapse = ";")
      ) %>%
      dplyr::ungroup() %>%
      dplyr::distinct()

    # history (not used yet but could be added in the future)

    history <- pre_extracted_result %>%
      dplyr::distinct(.data$main_id, .data$history) %>%
      tidyr::unnest("history")

    # relations

    relations <- pre_extracted_result %>%
      dplyr::distinct(.data$main_id, .data$x_relations) %>%
      tidyr::unnest("x_relations") %>%
      dplyr::rename(
        chebi_id = "id",
        relations_term = "term",
        database = "namespace",
        relations_url = "url",
        relations_relation = "relation"
      )

    result <- pre_extracted_result %>%
      dplyr::distinct(
        .data$main_id,
        .data$is_obsolete,
        .data$main_name,
        .data$definition,
        .data$ontology,
        .data$usage
      ) %>%
      left_join(children, by = "main_id") %>%
      left_join(relations, by = "main_id")

    return(result)
  }

  # Retrieve slimed dataset
  if (type == "slims") {
    if (!missing(id_annotations)) {
      warning(strwrap('You provided IDs to "id_annotations" but chose to retrieve "slims" and not "annotations".
              Therefore, IDs will not be considered.', prefix = "\n", initial = ""))
    }
    if (!missing(taxon_id_annotations)) {
      warning(strwrap('You provided an ID to "taxon_id_annotations" but chose to retrieve "slims" and not "annotations".
              Therefore, the ID will not be considered.', prefix = "\n", initial = ""))
    }

    base_url <- "https://www.ebi.ac.uk/QuickGO/services/ontology/go/slim?"
    # GO IDs for slim set
    if (missing(go_id_slims)) {
      stop('Please provide at least one GO ID to "go_id_slims" if you chose type = "slims"!')
    }
    go_id_slims <- paste0("slimsToIds=", paste0(go_id_slims, collapse = ","), "&")

    # Relations of slim set
    relations_slims <- match.arg(relations_slims, c("is_a", "part_of", "regulates", "occurs_in"), several.ok = TRUE)

    relations_slims <- paste0("relations=", paste0(relations_slims, collapse = ","))

    # final url
    url <- paste0(base_url, go_id_slims, relations_slims)

    query_result <- try_query(url,
      type = "application/json",
      timeout = timeout,
      max_tries = max_tries,
      simplifyDataFrame = TRUE
    )

    if (methods::is(query_result, "character")) {
      message(query_result)
      if (stringr::str_detect(query_result, pattern = "Timeout")){
        message('Consider increasing the "timeout" or "max_tries" argument. \n')
      }
      return(invisible(NULL))
    } else {
      query_result <- query_result[["results"]] %>%
        janitor::clean_names() %>%
        tidyr::unnest("slims_to_ids")
    }

    return(query_result)
  }
}
