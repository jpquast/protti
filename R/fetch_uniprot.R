#' Fetch protein data from UniProt
#'
#' Fetches protein metadata from UniProt.
#'
#' @param uniprot_ids A character vector of UniProt accession numbers
#' @param columns Metadata columns that should be imported from UniProt (all possible columns can be found here: https://www.uniprot.org/help/uniprotkb_column_names)
#' @param batchsize Size of batch of proteins for a single query
#' @param show_progress Logical, if true, a progress bar will be shown
#'
#' @return A data frame that contains all protein metadata specified in \code{columns} for the proteins provided.
#' @export
#'
#' @examples
#' \dontrun{
#' fetch_uniprot(c("P36578" "O43324" "Q00796"))
#' }
fetch_uniprot <-
  function (uniprot_ids,
            columns = c(
              "protein names",
              "genes",
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
            batchsize = 400,
            show_progress = TRUE)
  {
    columns <- c("id", columns)
    column_names <- janitor::make_clean_names(columns)
    collapsed_columns <- paste(columns, collapse = ",")
    uniprot_ids <- na.omit(uniprot_ids)
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
    pb <- progress::progress_bar$new(total = length(batches))
    result <- purrr::map_df(batches, function(x) {
      if(show_progress == TRUE) {
        pb$tick()}
      id_query <- paste(paste0("id:", x), collapse = "+or+")
      query_url <- URLencode(paste0(url, id_query, "&format=tab&columns=",
                                    collapsed_columns))
      try_query(query_url)
    })
    colnames(result) <- column_names
    result
  }