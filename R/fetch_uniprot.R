#' Fetch protein data from UniProt
#'
#' Fetches protein metadata from UniProt.
#'
#' @param uniprot_ids A character vector of UniProt accession numbers
#' @param columns Metadata columns that should be imported from UniProt (all possible columns can be found \href{https://www.uniprot.org/help/uniprotkb_column_names}{here}.)
#' @param batchsize Size of batch of proteins for a single query
#' @param show_progress Logical, if true, a progress bar will be shown
#'
#' @return A data frame that contains all protein metadata specified in \code{columns} for the proteins provided.
#' @import dplyr
#' @import janitor
#' @import progress
#' @import purrr
#' @importFrom tidyr drop_na
#' @importFrom stringr str_detect
#' @importFrom stringr str_extract_all
#' @export
#'
#' @examples
#' \dontrun{
#' fetch_uniprot(c("P36578", "O43324", "Q00796"))
#' }
fetch_uniprot <-
  function (uniprot_ids,
            columns = c(
              "protein names",
              "genes",
              "database(GeneID)",
              "database(String)",
              "go(molecular function)",
              "interactor",
              "feature(ACTIVE SITE)",
              "feature(BINDING SITE)",
              "feature(METAL BINDING)",
              "chebi(Cofactor)",
              "chebi(Catalytic activity)",
              "feature(MODIFIED RESIDUE)",
              "feature(LIPIDATION)",
              "feature(GLYCOSYLATION)",
              "database(PDB)",
              "length",
              "sequence"
            ),
            batchsize = 200,
            show_progress = TRUE)
  {
    . = NULL
    columns <- c("id", columns)
    column_names <- janitor::make_clean_names(columns)
    collapsed_columns <- paste(columns, collapse = ",")
    uniprot_ids <- stats::na.omit(uniprot_ids)
    id_test <- stringr::str_detect(uniprot_ids, pattern = "^[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}$")
    uniprot_ids_filtered <- uniprot_ids[id_test]
    if (sum(!id_test) > 0) {
      non_conform_ids <- uniprot_ids[!id_test]
      warning("These uniprot accession numbers did not conform to uniprot standards and were skipped from importing:", paste(non_conform_ids, collapse = ", "))
    }
    if (length(uniprot_ids_filtered) == 0) {
      stop("No valid UniProt accession numbers found.")
    }
    url <- "http://www.uniprot.org/uniprot/?query="
    batches <- split(uniprot_ids_filtered, ceiling(seq_along(uniprot_ids_filtered)/batchsize))
    pb <- progress::progress_bar$new(total = length(batches), show_after = 0)
    pb$tick(0)
    result <- purrr::map_df(batches, function(x) {
      id_query <- paste(paste0("id:", x), collapse = "+or+")
      query_url <- utils::URLencode(paste0(url, id_query, "&format=tab&columns=",
                                    collapsed_columns))
      query <- try_query(query_url)
      if(show_progress == TRUE) {
        pb$tick()}
      return(query)
    })
    colnames(result) <- column_names

    # rescue cases in which the used ID is an old ID. The information is retreived for the new ID and then parsed to the old ID.

    new <- result %>%
      dplyr::mutate(new = ifelse(stringr::str_detect(.[[2]], pattern = "Merged"), stringr::str_extract_all(.[[2]], pattern = "(?<=into )[A-Z0-9]+", simplify = TRUE), NA)) %>%
      dplyr::distinct(id, new) %>%
      tidyr::drop_na(new)

    new_ids <- new$new

    if(all.equal(new_ids, logical(0)) == TRUE) {return(result)}

    new_id_query <- paste(paste0("id:", new_ids), collapse = "+or+")
    new_query_url <- utils::URLencode(paste0(url, new_id_query, "&format=tab&columns=",
                                             collapsed_columns))

    new_result <- try_query(new_query_url)
    colnames(new_result) <- column_names

    new_result <- new %>%
      dplyr::left_join(new_result, by = c("new" = "id")) %>%
      dplyr::select(-new)

    result <- result %>%
      dplyr::filter(!stringr::str_detect(.[[2]], pattern = "Merged")) %>%
      dplyr::bind_rows(new_result)

    result
  }
