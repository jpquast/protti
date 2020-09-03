#' Fetch KEGG pathway data from KEGG
#'
#' Fetches gene IDs and corresponding pathway IDs and names for the provided organism.
#'
#' @param species A character vector providing an abreviated species name. "hsa" for human, "eco" for E. coli and "sce" for S. cerevisiae. 
#' Additional possible names can be found at https://www.genome.jp/kegg-bin/show_organism?category=Eukaryotes for eukaryotes and at 
#' https://www.genome.jp/kegg-bin/show_organism?category=Prokaryotes for prokaryotes. 
#'
#' @return A data frame that contains gene IDs with corresponding pathway IDs and names for a selected organism.
#' @importFrom dplyr left_join
#' @importFrom stringr str_replace_all
#' @importFrom magrittr %>%
#' @export
#'
#' @examples
#' \dontrun{
#' fetch_kegg(species = "hsa")
#' }
fetch_kegg <- function(species){
  # download kegg_id pathway link
  url_link <- paste("http://rest.kegg.jp/link/pathway", species, sep = "/")
  result_link <- try_query(url_link, header = FALSE)
  colnames(result_link) <- c("kegg_id", "pathway_id")
  # download pathway_id names
  url_name <- paste("http://rest.kegg.jp/list/pathway", species, sep = "/")
  result_name <- try_query(url_name, header = FALSE)
  colnames(result_name) <- c("pathway_id", "pathway_name")
  # download kegg_id to uniprot_id conversion
  url_conv <- paste("http://rest.kegg.jp/conv/uniprot", species, sep = "/")
  result_conv <- try_query(url_conv, header = FALSE)
  colnames(result_conv) <- c("kegg_id", "uniprot_id")
  result_conv$uniprot_id <- stringr::str_replace_all(result_conv$uniprot, pattern = "up:", replacement = "")
  # combine datasets
  result <- result_link %>%
    dplyr::left_join(result_name, by = "pathway_id") %>%
    dplyr::left_join(result_conv, by = "kegg_id")
  result
}