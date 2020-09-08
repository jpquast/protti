#' Fetch protein disorder information from MobiDB
#'
#' Fetches information about disordered protein regions from MobiDB. 
#'
#' @param organism_id An NCBI taxonomy identifier of an organism (TaxId). Possible inputs inlude only: "9606" (Human), "559292" (Yeast), 
#' "83333" (E. coli), "10090" (Mouse), "9913" (Bovine), "7227" (Fruit fly). 
#' @param protein_ids A character vector of UniProt identifiers. These need to be proteins from the organism provided in \code{organism_id}.
#'
#' @return A data frame that contains start and end positions for disordered (D), structured (S) or conflict/context-dependent (C) regions 
#' and their source for each protein provided. A full explanation of sources and regions can be found here: https://mobidb.bio.unipd.it/schema. 
#' @importFrom jsonlite stream_in
#' @importFrom rlang .data
#' @importFrom purrr map map_df map2_df
#' @importFrom dplyr select filter mutate rename
#' @importFrom tibble as_tibble
#' @importFrom tidyr unnest
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

  organism_file <- switch(organism_id, "9606" = "disorder_UP000005640.mjson.gz", 
                          "559292" = "disorder_UP000002311.mjson.gz", 
                          "83333" = "disorder_UP000000625.mjson.gz",
                          "10090" = "disorder_UP000000589.mjson.gz",
                          "9913" = "disorder_UP000009136.mjson.gz",
                          "7227" = "disorder_UP000000803.mjson.gz")
  
  mobidb <- jsonlite::stream_in(gzcon(url(paste0("https://mobidb.bio.unipd.it/mobidb3_datasets/latest/", organism_file))))
  
  result <- mobidb %>%
    dplyr::select(.data$acc, .data$mobidb_consensus) %>%
    dplyr::filter(.data$acc %in% protein_ids) %>%
    split(.$acc) %>%
    purrr::map(.f = ~ tidyr::unnest(.x, .data$mobidb_consensus)) %>%
    purrr::map(.f = ~ list(db = .x[["db"]][[1]],
                    derived = .x[["derived"]][[1]],
                    predict = .x[["predictors"]][[1]]
                    )
    ) %>%
    purrr::map(.f = ~ list(
      derived = tibble::as_tibble(
          .x[["derived"]][.x[["derived"]]$method == "full",][["regions"]][[1]]
      ),
      db = tibble::as_tibble(
          .x[["db"]][["regions"]][[1]]
      ),
      predictors = tibble::as_tibble(
        purrr::map_df(.x[["predict"]][["regions"]], .f = ~ tibble::as_tibble(.x))
      )
    )) %>%
    purrr::map(.f = ~ purrr::map2_df(.x, 
                                     .y = names(.),
                                     .f = ~ dplyr::mutate(.x, source = .y))) %>%
    purrr::map2_df(.y = names(.),
            .f = ~ dplyr::mutate(.x, protein_id = .y)) %>%
    dplyr::rename(start = .data$V1, end = .data$V2, type = .data$V3) %>%
    dplyr::mutate(start = as.numeric(.data$start),
           end = as.numeric(.data$end))
  
  if(length(unique(result$protein_id)) != length(protein_ids)) {
    not_included <- protein_ids[!protein_ids %in% unique(result$protein_id)]
    warning(paste(length(not_included), "proteins have no information about disordered regions. These include:", paste(not_included, collapse = ", ")))
  }
  result
}