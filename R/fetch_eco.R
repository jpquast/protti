#' Fetch evidence & conclusion ontology
#'
#' Fetches all evidence & conclusion ontology (ECO) information from the QuickGO EBI database. The ECO project is
#' maintained through a public \href{https://github.com/evidenceontology/evidenceontology}{GitHub repository}.
#'
#' According to the GitHub repository ECO is defined as follows:
#'
#' "The Evidence & Conclusion Ontology (ECO) describes types of scientific evidence within the
#' biological research domain that arise from laboratory experiments, computational methods,
#' literature curation, or other means. Researchers use evidence to support conclusions
#' that arise out of scientific research. Documenting evidence during scientific research
#' is essential, because evidence gives us a sense of why we believe what we think we know.
#' Conclusions are asserted as statements about things that are believed to be true, for
#' example that a protein has a particular function (i.e. a protein functional annotation) or
#' that a disease is associated with a particular gene variant (i.e. a phenotype-gene association).
#' A systematic and structured (i.e. ontological) classification of evidence allows us to store,
#' retreive, share, and compare data associated with that evidence using computers, which are
#' essential to navigating the ever-growing (in size and complexity) corpus of scientific
#' information."
#'
#' More information can be found in their
#' \href{https://academic.oup.com/nar/article/47/D1/D1186/5165344?login=true}{publication}.
#'
#' @param return_relation a logical value that indicates if relational information should be returned instead
#' the main descriptive information. This data can be used to check the relations of ECO terms to each other.
#' Default is FALSE.
#' @param return_history a logical value that indicates if the entry history of an ECO term should be
#' returned instead the main descriptive information.
#' Default is FALSE.
#' @param show_progress a logical value that indicates if a progress bar will be shown.
#' Default is TRUE.
#'
#' @return A data frame that contains descriptive information about each ECO term in the EBI database.
#' If either \code{return_relation} or \code{return_history} is set to \code{TRUE}, the respective information is
#' returned instead of the usual output.
#' @import dplyr
#' @import progress
#' @import purrr
#' @importFrom tidyr unnest
#' @importFrom janitor clean_names
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @importFrom curl has_internet
#' @export
#'
#' @examples
#' \donttest{
#' eco <- fetch_eco()
#'
#' head(eco)
#' }
fetch_eco <- function(return_relation = FALSE,
                      return_history = FALSE,
                      show_progress = TRUE) {
  # Make sure that not both return_history and return_relation are TRUE

  if (return_history & return_relation) {
    stop("Please only set either return_history or return_relation to TRUE and not both!")
  }

  # Check if there is an internet connection

  if (!curl::has_internet()) {
    message("No internet connection.")
    return(invisible(NULL))
  }

  base_url <- "https://www.ebi.ac.uk/QuickGO/services/ontology/eco/terms?page="

  # first we request one batch of test data to find the maximum number of pages

  page <- "1"

  url <- paste0(base_url, page)

  test_query <- try_query(url,
    type = "application/json"
  )

  if (methods::is(test_query, "character")) {
    message(test_query)
    return(invisible(NULL))
  }

  pages <- 1:test_query[["pageInfo"]][["total"]]

  if (show_progress == TRUE) {
    pb <- progress::progress_bar$new(
      total = length(pages),
      format = "  Fetching ECO IDs [:bar] (:percent) :eta",
      show_after = 0
    )
    pb$tick(0)
  }

  query_result_list <- purrr::map(
    .x = pages,
    .f = function(x) {
      query_url <- paste0(
        base_url,
        x
      )
      query <- try_query(query_url,
        type = "application/json",
        simplifyDataFrame = TRUE
      )

      if (show_progress == TRUE) {
        pb$tick()
      }
      query
    }
  )

  # if any page was not fetched correctly return NULL and the error message

  error_vector <- query_result_list %>%
    purrr::keep(.p = ~ is.character(.x)) %>%
    unlist()

  if (!is.null(error_vector)) {
    message(paste(unique(error_vector), collapse = ", "))
    return(invisible(NULL))
  } else {
    query_result <- query_result_list %>%
      purrr::map_dfr(
        .f = ~ .x[["results"]]
      )
  }

  if (return_history) {
    # Extract history information

    history <- query_result %>%
      dplyr::distinct(.data$id, .data$history) %>%
      tidyr::unnest("history") %>%
      dplyr::distinct(.data$id, .data$timestamp, .data$action, .data$category, .data$text)

    return(history)
  }

  if (return_relation) {
    # Extract relation information

    relation <- query_result %>%
      dplyr::distinct(.data$id, .data$children) %>%
      dplyr::rename(main_id = "id") %>%
      tidyr::unnest("children") %>%
      dplyr::rename(child_id = "id")

    return(relation)
  }

  # Unnest the data frame bit by bit to not lose information

  query_result_unnest_1 <- query_result %>%
    dplyr::select(-c("history", "children")) %>%
    tidyr::unnest("definition") %>%
    dplyr::rename(main_name = "name") %>%
    dplyr::rename(definition = "text")

  query_result_unnest_2 <- query_result_unnest_1 %>%
    tidyr::unnest("secondaryIds") %>%
    dplyr::distinct(.data$id, .data$secondaryIds) %>%
    dplyr::right_join(dplyr::select(query_result_unnest_1, c(-"secondaryIds")), by = "id")

  query_result_unnest_3 <- query_result_unnest_2 %>%
    tidyr::unnest("xRefs") %>%
    dplyr::distinct(.data$id, .data$dbCode, .data$dbId) %>%
    dplyr::right_join(dplyr::select(query_result_unnest_2, c(-"xRefs")), by = "id")

  result <- query_result_unnest_3 %>%
    tidyr::unnest("synonyms") %>%
    dplyr::distinct(.data$id, .data$name, .data$type) %>%
    dplyr::right_join(dplyr::select(query_result_unnest_3, c(-"synonyms")), by = "id", relationship = "many-to-many") %>%
    dplyr::select(
      "id",
      "isObsolete",
      "main_name",
      "definition",
      "comment",
      "name",
      "type",
      "dbCode",
      "dbId",
      "secondaryIds"
    ) %>%
    janitor::clean_names()

  result
}
