#' Find all ChEBI sub IDs of an ID
#'
#' For a given ChEBI ID, find all ChEBI sub IDs (incoming IDs) and their sub IDs. The type of
#' relationship can be selected too. This is a helper function for other functions.
#'
#' @param data a data frame that contains relational information on ChEBI IDs (id), their sub
#' IDs (incoming) and their relationship (type). This data frame can be obtained by calling
#' \code{fetch_chebi(relation = TRUE)}.
#' @param ids a character vector of ChEBI IDs for which sub IDs should be searched.
#' @param types a character vector containing the types of relationships that should be considered
#' for the search. It is possible to use "all" relationships. The default type is "is_a". A list of
#' possible relationships can be found
#' \href{https://docs.google.com/document/d/1_w-DwBdCCOh1gMeeP6yqGzcnkpbHYOa3AGSODe5epcg/edit#heading=h.hnsqoqu978s5}{here}.
#' @param exclude_parent_id a logical value that specifies if the parent ID should be included in 
#' the returned list.
#'
#' @return A list of double vectors containing the provided ID and all of its sub IDs. It
#' contains one element per input ID.
#' @importFrom dplyr select filter
#' @importFrom magrittr %>%
#' @importFrom purrr map
find_all_subs <- function(data, ids, types = "is_a", exclude_parent_id = FALSE) {
  if (!requireNamespace("igraph", quietly = TRUE)) {
    message("Package \"igraph\" is needed for this function to work. Please install it.", call. = FALSE)
    return(invisible(NULL))
  }
  if (ifelse(length(types) == 1, types == "all", FALSE)) {
    data <- data %>%
      dplyr::select(-.data$type)
  } else {
    data <- data %>%
      dplyr::filter(.data$type %in% types) %>%
      dplyr::select(-.data$type)
  }
  # Generate graph
  g <- igraph::graph_from_data_frame(data, directed = TRUE)

  result <- purrr::map(ids, function(x) {
    if (!(x %in% data$id)) {
      return(NULL)
    }
    r <- igraph::subcomponent(g, match(x, igraph::V(g)$name), "out")$name
    if(exclude_parent_id){
      r <- r[r != x]
    }
    
    as.numeric(r)
  })
  result
}
