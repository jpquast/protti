#' Fetch protein disorder information from MobiDB
#'
#' Fetches information about disordered protein regions from MobiDB.
#'
#' @param organism_id a character value that specifies the NCBI taxonomy identifier of an organism
#' (TaxId). Possible inputs inlude only: "9606" (Human), "559292" (Yeast), "83333" (E. coli),
#' "10090" (Mouse), "9913" (Bovine), "7227" (Fruit fly).
#' @param protein_ids a character vector of UniProt identifiers. These need to be proteins from
#' the organism provided in \code{organism_id}.
#'
#' @return A data frame that contains start and end positions for disordered regions for each
#' protein provided. The \code{feature} column contains information on the source of this
#' annotation. More information on the source can be found
#' \href{https://mobidb.bio.unipd.it/about/mobidb}{here}.
#' @importFrom rlang .data
#' @importFrom dplyr filter mutate
#' @importFrom tidyr unnest separate
#' @importFrom stringr str_split
#' @export
#'
#' @examples
#' \donttest{
#' fetch_mobidb(
#'   organism_id = "83333",
#'   protein_ids = c("P0A799", "P62707")
#' )
#' }
fetch_mobidb <- function(organism_id, protein_ids) {
  if (!requireNamespace("httr", quietly = TRUE)) {
    message("Package \"httr\" is needed for this function to work. Please install it.", call. = FALSE)
    return(invisible(NULL))
  }
  . <- NULL

  organism_id <- match.arg(organism_id, c("9606", "559292", "83333", "10090", "9913", "7227"))

  query <- paste0(
    "https://mobidb.bio.unipd.it/api/download?ncbi_taxon_id=",
    organism_id,
    "&projection=prediction-disorder-mobidb_lite,curated-disorder-merge",
    ",derived-missing_residues-th_90,derived-mobile_residues-th_90,acc,name&format=tsv"
  )

  query_result <- httr::GET(query, config = httr::config(connecttimeout = 60))

  mobidb <- suppressMessages(httr::content(query_result,
    type = "text/tab-separated-values",
    encoding = "UTF-8"
  ))

  i <- 0
  while (("ERROR: operation exceeded time limit" %in% mobidb$acc) & i < 4) {
    message("Attempt to download data timed out. Trying again")
    query_result <- httr::GET(query, config = httr::config(connecttimeout = 60))

    mobidb <- suppressMessages(httr::content(query_result,
      type = "text/tab-separated-values",
      encoding = "UTF-8"
    ))

    i <- i + 1

    Sys.sleep(3)
  }

  if ("ERROR: operation exceeded time limit" %in% mobidb$acc) {
    stop("Query operation exceeded server time limit. /n Try running it again later.")
  }

  result <- mobidb %>%
    dplyr::filter(.data$acc %in% protein_ids) %>%
    dplyr::mutate(start..end = stringr::str_split(.data$start..end, pattern = ",")) %>%
    tidyr::unnest(.data$start..end) %>%
    tidyr::separate(.data$start..end, into = c("start", "end"), sep = "\\.\\.") %>%
    dplyr::select(-.data$length)

  if (length(unique(result$acc)) != length(protein_ids)) {
    not_included <- protein_ids[!protein_ids %in% unique(result$acc)]
    warning(paste(
      length(not_included),
      "proteins have no information about disordered regions. These include:",
      paste(not_included, collapse = ", ")
    ))
  }
  result
}
