#' Find ChEBI IDs for name patterns
#'
#' A list of ChEBI IDs that contain a specific name pattern is returned.
#'
#' @param chebi_data A data frame that contains at least information on ChEBI IDs (id) and their names (name).
#' This data frame can be obtained by calling \code{fetch_chebi()}. Ideally this should be subsetted to only contain molecules of a specific type e.g. metals.
#' This can be achieved by calling \code{find_all_subs} with a general ID such as "25213" (Metal cation) and then subset the complete ChEBI database to only include the returned sub-IDs.
#' Using a subsetted database ensures better search results. This is a helper function for other functions.
#' @param pattern A character vector that contains names or name patterns of molecules. Name patterns can be for example obtained with the \code{split_metal_name} function.
#'
#' @return A list of character vectors containing ChEBI IDs that have a name matching the supplied pattern.
#' @importFrom dplyr distinct
#' @importFrom magrittr %>%
#' @importFrom purrr map
#' @importFrom stringr str_detect regex
#' @importFrom rlang .data
find_chebis <- function(chebi_data, pattern) {
  if (!requireNamespace("stringi", quietly = TRUE)) {
    stop("Package \"stringi\" is needed for this function to work. Please install it.", call. = FALSE)
  }
  data <- chebi_data %>%
    dplyr::distinct(.data$id, .data$name)

  purrr::map(pattern, function(pattern) {
    stringi::stri_remove_empty(stats::na.omit(unique(
      ifelse(
        stringr::str_detect(data$name, pattern = stringr::regex({{ pattern }}, ignore_case = TRUE)),
        data$id,
        ""
      )
    )))
  })
}
