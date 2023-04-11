#' Fetch KEGG pathway data from KEGG
#'
#' Fetches gene IDs and corresponding pathway IDs and names for the provided organism.
#'
#' @param species a character value providing an abreviated species name. "hsa" for human, "eco"
#' for E. coli and "sce" for S. cerevisiae. Additional possible names can be found for
#' \href{https://www.genome.jp/kegg-bin/show_organism?category=Eukaryotes}{eukaryotes} and for
#' \href{https://www.genome.jp/kegg-bin/show_organism?category=Prokaryotes}{prokaryotes}.
#'
#' @return A data frame that contains gene IDs with corresponding pathway IDs and names for a
#' selected organism.
#' @importFrom dplyr left_join
#' @importFrom stringr str_replace_all
#' @importFrom magrittr %>%
#' @importFrom curl has_internet
#' @export
#'
#' @examples
#' \donttest{
#' kegg <- fetch_kegg(species = "hsa")
#'
#' head(kegg)
#' }
fetch_kegg <- function(species) {
  if (!curl::has_internet()) {
    message("No internet connection.")
    return(invisible(NULL))
  }
  # download kegg_id pathway link
  url_link <- paste("https://rest.kegg.jp/link/pathway", species, sep = "/")
  result_link <- try_query(url_link, col_names = FALSE, progress = FALSE, show_col_types = FALSE)
  if (methods::is(result_link, "character")) {
    message(result_link)
    return(invisible(NULL))
  }
  colnames(result_link) <- c("kegg_id", "pathway_id")
  result_link$pathway_id <- stringr::str_replace_all(result_link$pathway_id,
                                                     pattern = "path:",
                                                     replacement = ""
  )
  # download pathway_id names
  url_name <- paste("https://rest.kegg.jp/list/pathway", species, sep = "/")
  result_name <- try_query(url_name, col_names = FALSE, progress = FALSE, show_col_types = FALSE)
  if (methods::is(result_name, "character")) {
    message(result_name)
    return(invisible(NULL))
  }
  colnames(result_name) <- c("pathway_id", "pathway_name")
  
  # download kegg_id to uniprot_id conversion
  url_conv <- paste("https://rest.kegg.jp/conv/uniprot", species, sep = "/")
  result_conv <- try_query(url_conv, col_names = FALSE, progress = FALSE, show_col_types = FALSE)
  if (methods::is(result_conv, "character")) {
    message(result_conv)
    return(invisible(NULL))
  }
  colnames(result_conv) <- c("kegg_id", "uniprot_id")
  result_conv$uniprot_id <- stringr::str_replace_all(result_conv$uniprot_id,
    pattern = "up:",
    replacement = ""
  )
  # combine datasets
  result <- result_link %>%
    dplyr::left_join(result_name, by = "pathway_id") %>%
    dplyr::left_join(result_conv, by = "kegg_id")
  result
}
