#' Find all ChEBI sub IDs of an ID
#'
#' For a given ChEBI ID, find all ChEBI sub IDs (incoming IDs) and their sub IDs. The type of relationship can be selected too. This is a helper function for other functions.
#'
#' @param data A data frame that contains information on ChEBI IDs (id), their sub IDs (incoming) and their relationship (type). This data frame can be obtained by calling \code{fetch_chebi(relation = TRUE)}.
#' @param id A character vector of ChEBI IDs for which sub IDs should be retreived.
#' @param type A character vector containing the type of relationship that should be considered for retreival. It is possible to use "all" relationships. The default type is "is_a". A list of possible relationships can be found \href{https://docs.google.com/document/d/1_w-DwBdCCOh1gMeeP6yqGzcnkpbHYOa3AGSODe5epcg/edit#heading=h.hnsqoqu978s5}{here}.
#'
#' @return A list of character vector containing the provided ID and all of its sub IDs. It contains one element per input ID.
#' @importFrom dplyr select filter
#' @importFrom magrittr %>%
#' @importFrom purrr map
find_all_subs <- function(data, id, type = "is_a") {
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("Package \"igraph\" is needed for this function to work. Please install it.", call. = FALSE)
  }
  if (type == "all") {
    data <- data %>%
      dplyr::select(-type)
  } else {
    data <- data %>%
      dplyr::filter(type == {{ type }}) %>%
      dplyr::select(-type)
  }
  result <- purrr::map(id, function(id) {
    if (!(id %in% data$id)) {
      return(NULL)
    }
    g <- igraph::graph_from_data_frame(data, directed = TRUE)
    r <- igraph::subcomponent(g, match(id, igraph::V(g)$name), "out")$name
    as.numeric(r)
  })
  result
}
