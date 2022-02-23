#' Find ChEBI IDs for name patterns
#'
#' Search for chebi IDs that match a specific name pattern. A list of corresponding ChEBI IDs is
#' returned.
#'
#' @param chebi_data a data frame that contains at least information on ChEBI IDs (id) and their
#' names (name). This data frame can be obtained by calling \code{fetch_chebi()}. Ideally this
#' should be subsetted to only contain molecules of a specific type e.g. metals. This can be
#' achieved by calling \code{find_all_subs} with a general ID such as "25213" (Metal cation) and
#' then subset the complete ChEBI database to only include the returned sub-IDs. Using a subsetted
#' database ensures better search results. This is a helper function for other functions.
#' @param pattern a character vector that contains names or name patterns of molecules. Name
#' patterns can be for example obtained with the \code{split_metal_name} function.
#'
#' @return A list of character vectors containing ChEBI IDs that have a name matching the supplied
#' pattern. It contains one element per pattern.
#' @importFrom dplyr distinct
#' @importFrom magrittr %>%
#' @importFrom purrr map
#' @importFrom stringr str_detect regex
#' @importFrom rlang .data
find_chebis <- function(chebi_data, pattern) {
  if (!requireNamespace("stringi", quietly = TRUE)) {
    message("Package \"stringi\" is needed for this function to work. Please install it.", call. = FALSE)
    return(invisible(NULL))
  }
  data <- chebi_data %>%
    dplyr::distinct(.data$id, .data$name)

  purrr::map(pattern, function(x) {
    stringi::stri_remove_empty(stats::na.omit(unique(
      ifelse(
        stringr::str_detect(data$name,
          pattern = stringr::regex(
            x,
            ignore_case = TRUE
          )
        ),
        data$id,
        ""
      )
    )))
  })
}
