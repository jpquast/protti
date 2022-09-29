#' Fetch ChEBI database information
#'
#' Fetches information from the ChEBI database.
#'
#' @param relation a logical value that indicates if ChEBI Ontology data will be returned instead
#' the main compound data. This data can be used to check the relations of ChEBI ID's to each other.
#' Default is FALSE.
#' @param stars a numeric vector indicating the "star" level (confidence) for which entries should
#' be retrieved (Possible levels are 1, 2 and 3). Default is \code{c(3)} retrieving only "3-star"
#' entries, which are manually annotated by the ChEBI curator team.
#'
#' @return A data frame that contains information about each molecule in the ChEBI database.
#' @import dplyr
#' @importFrom httr GET config content
#' @importFrom tidyr pivot_wider unnest
#' @importFrom janitor clean_names
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @export
#'
#' @examples
#' \donttest{
#' chebi <- fetch_chebi()
#'
#' head(chebi)
#' }
fetch_chebi <- function(relation = FALSE, stars = c(3)) {
  if (!curl::has_internet()) {
    message("No internet connection.")
    return(invisible(NULL))
  }
  # Retrieve relational information
  if (relation == TRUE) {
    chebi_relation_result <- tryCatch(httr::GET(
      "ftp://ftp.ebi.ac.uk/pub/databases/chebi/Flat_file_tab_delimited/relation.tsv",
      httr::timeout(60)
    ),
    error = function(e) conditionMessage(e),
    warning = function(w) conditionMessage(w)
    )
    # Check again if there is an internet connection
    if (!curl::has_internet()) {
      message("No internet connection.")
    }
    # check if response is error
    if (!methods::is(chebi_relation_result, "response")) {
      message(chebi_relation_result)
      return(invisible(NULL))
    }

    if (httr::http_error(chebi_relation_result)) {
      httr::message_for_status(chebi_relation_result)
      return(invisible(NULL))
    }

    chebi_relation <- suppressMessages(httr::content(chebi_relation_result,
      type = "text/tab-separated-values",
      progress = FALSE,
      show_col_types = FALSE
    ))

    chebi_relation_clean <- chebi_relation %>%
      dplyr::filter(.data$STATUS == "C") %>%
      dplyr::select(-c(.data$ID, .data$STATUS)) %>%
      dplyr::rename(incoming = .data$FINAL_ID) %>%
      dplyr::rename(ID = .data$INIT_ID) %>%
      janitor::clean_names()

    return(chebi_relation_clean)
  }

  # Download compound data
  chebi_chemical_data_result <- tryCatch(httr::GET(
    "ftp://ftp.ebi.ac.uk/pub/databases/chebi/Flat_file_tab_delimited/chemical_data.tsv",
    httr::timeout(60)
  ),
  error = function(e) conditionMessage(e),
  warning = function(w) conditionMessage(w)
  )

  # check if response is error
  if (!methods::is(chebi_chemical_data_result, "response")) {
    message(chebi_chemical_data_result)
    return(invisible(NULL))
  }

  if (httr::http_error(chebi_chemical_data_result)) {
    httr::message_for_status(chebi_chemical_data_result)
    return(invisible(NULL))
  }

  chebi_chemical_data <- suppressMessages(httr::content(chebi_chemical_data_result,
    type = "text/tab-separated-values",
    progress = FALSE,
    show_col_types = FALSE
  ))

  chebi_compounds_download <- tryCatch(readLines(con <- gzcon(url("ftp://ftp.ebi.ac.uk/pub/databases/chebi/Flat_file_tab_delimited/compounds.tsv.gz", method = "libcurl"))),
    error = function(e) conditionMessage(e),
    warning = function(w) conditionMessage(w)
  )
  close(con)
  chebi_compounds <- utils::read.delim(textConnection(chebi_compounds_download),
    quote = "", stringsAsFactors = FALSE
  )

  if (nrow(chebi_compounds) == 1) {
    message(chebi_compounds$V1)
    return(invisible(NULL))
  }
  if (nrow(chebi_compounds) == 0) {
    message("No information could be retrieved from the database, try again later!")
    return(invisible(NULL))
  }

  chebi_names_download <- tryCatch(readLines(con <- gzcon(url("ftp://ftp.ebi.ac.uk/pub/databases/chebi/Flat_file_tab_delimited/names.tsv.gz", method = "libcurl"))),
    error = function(e) conditionMessage(e),
    warning = function(w) conditionMessage(w)
  )
  close(con)
  chebi_names <- utils::read.delim(textConnection(chebi_names_download),
    quote = "", stringsAsFactors = FALSE
  )

  if (nrow(chebi_names) == 1) {
    message(chebi_names$V1)
    return(invisible(NULL))
  }
  if (nrow(chebi_names) == 0) {
    message("No information could be retrieved from the database, try again later!")
    return(invisible(NULL))
  }

  chebi_accession_result <- tryCatch(httr::GET(
    "ftp://ftp.ebi.ac.uk/pub/databases/chebi/Flat_file_tab_delimited/database_accession.tsv",
    httr::timeout(60)
  ),
  error = function(e) conditionMessage(e),
  warning = function(w) conditionMessage(w)
  )

  # check if response is error
  if (!methods::is(chebi_accession_result, "response")) {
    message(chebi_accession_result)
    return(invisible(NULL))
  }

  if (httr::http_error(chebi_accession_result)) {
    httr::message_for_status(chebi_accession_result)
    return(invisible(NULL))
  }

  chebi_accession <- suppressMessages(httr::content(chebi_accession_result,
    type = "text/tab-separated-values",
    progress = FALSE,
    show_col_types = FALSE
  ))

  # Create one file with all information after cleaning individual source files

  chebi_compounds_clean <- chebi_compounds %>%
    dplyr::filter(.data$STAR %in% stars) %>%
    dplyr::select(-c(
      .data$SOURCE,
      .data$NAME,
      .data$STATUS,
      .data$MODIFIED_ON,
      .data$CREATED_BY
    )) %>%
    dplyr::na_if("null") %>%
    dplyr::mutate(PARENT_ID = as.numeric(.data$PARENT_ID))


  chebi_accession_clean <- chebi_accession %>%
    dplyr::distinct(.data$COMPOUND_ID, .data$TYPE, .data$ACCESSION_NUMBER) %>%
    dplyr::rename(ID = .data$COMPOUND_ID) %>%
    dplyr::rename(TYPE_ACCESSION = .data$TYPE)

  chebi_chemical_data_clean <- chebi_chemical_data %>%
    dplyr::distinct(.data$COMPOUND_ID, .data$TYPE, .data$CHEMICAL_DATA) %>%
    dplyr::rename(ID = .data$COMPOUND_ID) %>%
    tidyr::pivot_wider(
      names_from = .data$TYPE,
      values_from = .data$CHEMICAL_DATA,
      values_fn = list
    ) %>%
    tidyr::unnest(cols = c(
      .data$FORMULA,
      .data$MASS,
      .data$CHARGE,
      .data$`MONOISOTOPIC MASS`
    ))

  chebi_compounds_names_clean <- chebi_compounds %>%
    dplyr::filter(.data$STAR %in% stars) %>%
    dplyr::distinct(.data$ID, .data$NAME) %>%
    dplyr::na_if("null") %>%
    dplyr::filter(!is.na(.data$NAME)) %>%
    dplyr::mutate(TYPE_NAME = "STANDARD") %>%
    dplyr::select(.data$ID, .data$TYPE_NAME, .data$NAME)

  chebi_names_clean <- chebi_names %>%
    dplyr::distinct(.data$COMPOUND_ID, .data$NAME, .data$TYPE) %>%
    dplyr::rename(ID = .data$COMPOUND_ID, TYPE_NAME = .data$TYPE) %>%
    dplyr::bind_rows(chebi_compounds_names_clean)

  chebi <- chebi_compounds_clean %>%
    dplyr::left_join(chebi_names_clean, by = "ID") %>%
    dplyr::left_join(chebi_accession_clean, by = "ID") %>%
    dplyr::left_join(chebi_chemical_data_clean, by = "ID")

  # Add info to old compound IDs

  parent_ids <- chebi %>%
    dplyr::filter(!is.na(.data$PARENT_ID)) %>%
    dplyr::pull(var = .data$PARENT_ID) %>%
    unique()

  parent_info <- chebi %>%
    dplyr::filter(.data$ID %in% parent_ids) %>%
    dplyr::select(c(
      .data$ID,
      .data$NAME,
      .data$TYPE_NAME,
      .data$DEFINITION,
      .data$TYPE_ACCESSION,
      .data$ACCESSION_NUMBER,
      .data$FORMULA,
      .data$MASS,
      .data$CHARGE,
      .data$`MONOISOTOPIC MASS`
    ))

  parent_complete <- chebi %>%
    dplyr::filter(!is.na(.data$PARENT_ID)) %>%
    dplyr::select(-c(
      .data$NAME,
      .data$TYPE_NAME,
      .data$DEFINITION,
      .data$TYPE_ACCESSION,
      .data$ACCESSION_NUMBER,
      .data$FORMULA,
      .data$MASS,
      .data$CHARGE,
      .data$`MONOISOTOPIC MASS`
    )) %>%
    dplyr::left_join(parent_info, by = c("PARENT_ID" = "ID"))

  chebi <- chebi %>%
    dplyr::filter(is.na(.data$PARENT_ID)) %>%
    dplyr::bind_rows(parent_complete) %>%
    dplyr::arrange(.data$ID) %>%
    janitor::clean_names()

  chebi
}
