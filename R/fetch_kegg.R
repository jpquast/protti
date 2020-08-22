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
  # download gene_id pathway link
  url_link <- paste("http://rest.kegg.jp/link/pathway", species, sep = "/")
  result_link <- try_query(url_link, header = FALSE)
  colnames(result_link) <- c("gene_id", "pathway_id")
  result_link$gene_id <- stringr::str_replace_all(result_link$gene_id, pattern = "hsa:", replacement = "")
  # download pathway_id names
  url_name <- paste("http://rest.kegg.jp/list/pathway", species, sep = "/")
  result_name <- try_query(url_name, header = FALSE)
  colnames(result_name) <- c("pathway_id", "pathway_name")
  # combine both datasets
  result <- result_link %>%
    dplyr::left_join(result_name, by = "pathway_id")
  result
}