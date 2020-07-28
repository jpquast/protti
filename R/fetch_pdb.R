#' Fetch protein structure data from RCSB
#'
#' Fetches protein structural data from the RCSB database.
#'
#' @param pdb_ids A character vector of pdb identifiers
#' @param columns A character vector of data columns that should be imported from RCSB (all possible columns can be found here: https://www.rcsb.org/pdb/results/reportField.do)
#' @param batchsize Size of batch of proteins for a single query
#' @param show_progress Logical, if true, a progress bar will be shown
#'
#' @return A data frame that contains all structural data specified in \code{columns} for the pdb structures provided.
#' @import janitor
#' @import progress
#' @import purrr
#' @import dplyr
#' @importFrom RCurl getURL
#' @importFrom data.table fread
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @export
#'
#' @examples
#' \dontrun{
#' fetch_pdb(c("3OFN", "3ZIA", "3ZRY"))
#' }
fetch_pdb <-
  function (pdb_ids,
            columns = c(
              "compound",
              "structureId",
              "entityId",
              "structureTitle",
              "chainId",
              "chainLength",
              "sequence",
              "classification",
              "db_name",
              "ligandName",
              "InChIKey",
              "uniprotAcc",
              "taxonomy",
              "taxonomyId",
              "experimentalTechnique",
              "resolution",
              "crystallizationMethod",
              "crystallizationTempK",
              "phValue",
              "method",
              "temperature",
              "ionicStrength",
              "ph"
            ), 
            batchsize = 400,
            show_progress = TRUE)
  {
    column_names <- janitor::make_clean_names(columns)
    collapsed_columns <- paste0("&customReportColumns=", paste(columns, collapse = ","))
    pdb_ids <- stats::na.omit(pdb_ids)
    id_test <- stringr::str_detect(pdb_ids, pattern = "[A-Z0-9]{4}")
    pdb_ids_filtered <- pdb_ids[id_test]
    if (sum(!id_test) > 0) {
      non_conform_ids <- pdb_ids[!id_test]
      warning("These PDB identifiers did not conform to the standards and were skipped from importing:", paste(non_conform_ids, collapse = ", "))
    }
    if (length(pdb_ids_filtered) == 0) {
      stop("No valid PDB identifier found.")
    }
    url <- "http://www.rcsb.org/pdb/rest/customReport.xml?"
    batches <- split(pdb_ids_filtered, ceiling(seq_along(pdb_ids_filtered)/batchsize))
    pb <- progress::progress_bar$new(total = length(batches))
    result <- purrr::map_df(batches, function(x) {
      if(show_progress == TRUE) {
        pb$tick()}
      id_query <- paste0("pdbids=", paste(x, collapse = ","))
      query_url <- utils::URLencode(paste0(url, id_query, collapsed_columns, "&service=wsfile&format=csv"))
      try_query(query_url)
    })
    colnames(result) <- column_names
    pdb_uniprot_mapping_url <- RCurl::getURL("ftp.ebi.ac.uk/pub/databases/msd/sifts/csv/pdb_chain_uniprot.csv")
    pdb_uniprot_mapping <- data.table::fread(pdb_uniprot_mapping_url) %>%
      janitor::clean_names() %>%
      dplyr::mutate(compound = stringr::str_to_upper(.data$pdb)) %>%
      dplyr::select(-.data$pdb)
    result <- result %>%
      dplyr::left_join(pdb_uniprot_mapping, by = c("compound", "structure_id" = "chain"))
    result
  }