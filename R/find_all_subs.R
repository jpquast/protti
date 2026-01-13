#' Find all sub IDs of an ID in a network
#'
#' For a given ID, find all sub IDs and their sub IDs etc. The type of
#' relationship can be selected too. This is a helper function for other functions.
#'
#' @param data a data frame that contains relational information on IDs (main_id) their sub
#' IDs (sub_id) and their relationship (type). For ChEBI this data frame can be obtained by calling
#' \code{fetch_chebi(relation = TRUE)}. For ECO data it can be obtained by calling fetch_eco(relation = TRUE).
#' @param ids a character vector of IDs for which sub IDs should be searched.
#' @param main_id a character or integer column containing IDs. Default is \code{id} for ChEBI IDs.
#' @param type a character column that contains the type of interactions. Default is \code{type} for ChEBI IDs.
#' @param accepted_types a character vector containing the accepted_types of relationships that should be considered
#' for the search. It is possible to use "all" relationships. The default type is "is_a". A list of
#' possible relationships for e.g. ChEBI IDs can be found
#' \href{https://docs.google.com/document/d/1_w-DwBdCCOh1gMeeP6yqGzcnkpbHYOa3AGSODe5epcg/edit#heading=h.hnsqoqu978s5}{here}.
#' @param exclude_parent_id a logical value that specifies if the parent ID should be included in
#' the returned list.
#'
#' @return A list of character vectors containing the provided ID and all of its sub IDs. It
#' contains one element per input ID.
#' @importFrom dplyr select filter pull
#' @importFrom magrittr %>%
#' @importFrom purrr map
#' @importFrom rlang .data
find_all_subs <- function(data,
                          ids,
                          main_id = id,
                          type = type,
                          accepted_types = "is_a",
                          exclude_parent_id = FALSE) {
  if (!requireNamespace("igraph", quietly = TRUE)) {
    message("Package \"igraph\" is needed for this function to work. Please install it.", call. = FALSE)
    return(invisible(NULL))
  }
  if (ifelse(length(accepted_types) == 1, accepted_types == "all", FALSE)) {
    data <- data %>%
      dplyr::select(-{{ type }})
  } else {
    data <- data %>%
      dplyr::filter({{ type }} %in% accepted_types) %>%
      dplyr::select(-{{ type }})
  }
  # Generate graph
  g <- igraph::graph_from_data_frame(data, directed = TRUE)

  result <- purrr::map(ids, function(x) {
    if (!(x %in% dplyr::pull(data, {{ main_id }}))) {
      return(NULL)
    }
    r <- igraph::subcomponent(g, match(x, igraph::V(g)$name), "out")$name
    if (exclude_parent_id) {
      r <- r[r != x]
    }

    r
  })
  result
}

utils::globalVariables("id")