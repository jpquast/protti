#' Fetch protein disorder information from MobiDB
#'
#' Fetches information about disordered protein regions from MobiDB. 
#'
#' @param organism_id An NCBI taxonomy identifier of an organism (TaxId). Possible inputs inlude only: "9606" (Human), "559292" (Yeast), 
#' "83333" (E. coli), "10090" (Mouse), "9913" (Bovine), "7227" (Fruit fly). 
#' @param protein_ids A character vector of UniProt identifiers. These need to be proteins from the organism provided in \code{organism_id}.
#'
#' @return A data frame that contains start and end positions for disordered regions for each protein provided. The \code{feature} column
#' contains information on the source of this annotation. More information on the source can be found \href{https://mobidb.bio.unipd.it/about/mobidb}{here}.
#' @importFrom httr GET content
#' @importFrom rlang .data 
#' @importFrom dplyr filter mutate
#' @importFrom tidyr unnest separate
#' @importFrom stringr str_split
#' @export
#'
#' @examples
#' \dontrun{
#' fetch_mobidb(
#' organism_id = "9606", 
#' protein_ids = c("P36578", "O43324", "Q00796")
#' )
#' }
fetch_mobidb <- function(organism_id, protein_ids){
  . = NULL
  
  organism_id <- match.arg(organism_id, c("9606", "559292", "83333", "10090", "9913", "7227"))

  query <- paste0("https://mobidb.bio.unipd.it/api/download?ncbi_taxon_id=", organism_id, "&projection=prediction-disorder-mobidb_lite,curated-disorder-merge,derived-missing_residues-th_90,derived-mobile_residues-th_90,acc,name,length&format=tsv")
  
  mobidb <- httr::GET(query) %>% 
    httr::content() 
  
  i  <- 0
  while (("ERROR: operation exceeded time limit" %in% mobidb$acc) & i < 4)
  {
    message("Attempt to download data timed out. Trying again")
    mobidb <- httr::GET(query) %>% 
      httr::content() 
    
    i <- i + 1
    
    Sys.sleep(3)
  }

  if("ERROR: operation exceeded time limit" %in% mobidb$acc){
    stop("Query operation exceeded server time limit. /n Try running it again later.")
  }
  
  result <- mobidb %>%
    dplyr::filter(.data$acc %in% protein_ids) %>% 
    dplyr::mutate(start..end = stringr::str_split(.data$start..end, pattern = ",")) %>%
    tidyr::unnest(.data$start..end) %>% 
    tidyr::separate(.data$start..end, into = c("start", "end"), sep = "\\.\\.")
  
  if(length(unique(result$acc)) != length(protein_ids)) {
    not_included <- protein_ids[!protein_ids %in% unique(result$acc)]
    warning(paste(length(not_included), "proteins have no information about disordered regions. These include:", paste(not_included, collapse = ", ")))
  }
  result
}