#' Fetch ChEBI database information
#'
#' Fetches information from the ChEBI database.
#'
#' @param relation a logical value that indicates if ChEBI Ontology data will be returned instead
#' the main compound data. This data can be used to check the relations of ChEBI ID's to each other.
#' Default is FALSE.
#' @param stars a numeric vector indicating the "star" level (confidence) for which entries should
#' be retrieved (Possible levels are 1, 2 and 3). Default is `c(3)` retrieving only "3-star"
#' entries, which are manually annotated by the ChEBI curator team.
#' @param timeout a numeric value specifying the time in seconds until the download of an organism
#' archive times out. The default is 60 seconds.
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
fetch_chebi <- function(relation = FALSE, stars = c(3), timeout = 60) {
  if (!curl::has_internet()) {
    message("No internet connection.")
    return(invisible(NULL))
  }
  # Retrieve relational information
  if (relation == TRUE) {
    con <- NULL
    chebi_relation_download <- tryCatch(
      readLines(
        con <- gzcon(
          url(
            "https://ftp.ebi.ac.uk/pub/databases/chebi/flat_files/relation.tsv.gz",
            method = "libcurl"
          )
        )
      ),
      error = function(e) conditionMessage(e),
      warning = function(w) conditionMessage(w),
      finally = if (!is.null(con)) close(con)
    )

    # Check again if there is an internet connection
    if (!curl::has_internet()) {
      message("No internet connection.")
    }

    chebi_relation <- utils::read.delim(
      textConnection(chebi_relation_download),
      sep = "\t",
      quote = "\"",
      fill = TRUE,
      comment.char = "",
      stringsAsFactors = FALSE
    )

    if (nrow(chebi_relation) == 1) {
      message(chebi_relation$V1)
      return(invisible(NULL))
    }
    if (nrow(chebi_relation) == 0) {
      message("No information could be retrieved from the database, try again later!")
      return(invisible(NULL))
    }

    con <- NULL
    chebi_relation_type_download <- tryCatch(
      readLines(
        con <- gzcon(
          url(
            "https://ftp.ebi.ac.uk/pub/databases/chebi/flat_files/relation_type.tsv.gz",
            method = "libcurl"
          )
        )
      ),
      error = function(e) conditionMessage(e),
      warning = function(w) conditionMessage(w),
      finally = if (!is.null(con)) close(con)
    )

    chebi_relation_type <- utils::read.delim(
      textConnection(chebi_relation_type_download),
      sep = "\t",
      quote = "\"",
      fill = TRUE,
      comment.char = "",
      stringsAsFactors = FALSE
    )

    if (nrow(chebi_relation_type) == 1) {
      message(chebi_relation$V1)
      return(invisible(NULL))
    }
    if (nrow(chebi_relation_type) == 0) {
      message("No information could be retrieved from the database, try again later!")
      return(invisible(NULL))
    }

    chebi_relation_type_clean <- chebi_relation_type %>%
      dplyr::select(c("id", "code")) %>%
      dplyr::rename(
        type = "code",
        relation_type_id = "id"
      )

    chebi_relation_clean <- chebi_relation %>%
      dplyr::filter(.data$status_id == 1) %>%
      dplyr::select(-c("id", "status_id", "evidence_source_id", "evidence_accession")) %>%
      dplyr::rename(incoming = "final_id") %>%
      dplyr::rename(id = "init_id") %>%
      left_join(chebi_relation_type_clean, by = "relation_type_id") %>%
      dplyr::select(-c("relation_type_id")) %>%
      janitor::clean_names()

    return(chebi_relation_clean)
  }

  # Download compound data
  con <- NULL
  chebi_chemical_data_download <- tryCatch(
    readLines(
      con <- gzcon(
        url(
          "https://ftp.ebi.ac.uk/pub/databases/chebi/flat_files/chemical_data.tsv.gz",
          method = "libcurl"
        )
      )
    ),
    error = function(e) conditionMessage(e),
    warning = function(w) conditionMessage(w),
    finally = if (!is.null(con)) close(con)
  )

  chebi_chemical_data <- utils::read.delim(
    textConnection(chebi_chemical_data_download),
    sep = "\t",
    quote = "\"",
    fill = TRUE,
    comment.char = "",
    stringsAsFactors = FALSE
  )

  if (nrow(chebi_chemical_data) == 1) {
    message(chebi_chemical_data$V1)
    return(invisible(NULL))
  }
  if (nrow(chebi_chemical_data) == 0) {
    message("No information could be retrieved from the database, try again later!")
    return(invisible(NULL))
  }


  con <- NULL
  chebi_compounds_download <- tryCatch(readLines(con <- gzcon(url("https://ftp.ebi.ac.uk/pub/databases/chebi/flat_files/compounds.tsv.gz", method = "libcurl"))),
    error = function(e) conditionMessage(e),
    warning = function(w) conditionMessage(w),
    finally = if (!is.null(con)) close(con)
  )

  chebi_compounds <- utils::read.delim(
    textConnection(chebi_compounds_download),
    sep = "\t",
    quote = "\"",
    fill = TRUE,
    comment.char = "",
    stringsAsFactors = FALSE
  )

  if (nrow(chebi_compounds) == 1) {
    message(chebi_compounds$V1)
    return(invisible(NULL))
  }
  if (nrow(chebi_compounds) == 0) {
    message("No information could be retrieved from the database, try again later!")
    return(invisible(NULL))
  }

  con <- NULL
  chebi_names_download <- tryCatch(readLines(con <- gzcon(url("https://ftp.ebi.ac.uk/pub/databases/chebi/flat_files/names.tsv.gz", method = "libcurl"))),
    error = function(e) conditionMessage(e),
    warning = function(w) conditionMessage(w),
    finally = if (!is.null(con)) close(con)
  )

  chebi_names <- utils::read.delim(
    textConnection(chebi_names_download),
    sep = "\t",
    quote = "\"",
    fill = TRUE,
    comment.char = "",
    stringsAsFactors = FALSE
  )

  if (nrow(chebi_names) == 1) {
    message(chebi_names$V1)
    return(invisible(NULL))
  }
  if (nrow(chebi_names) == 0) {
    message("No information could be retrieved from the database, try again later!")
    return(invisible(NULL))
  }

  con <- NULL
  chebi_accession_download <- tryCatch(
    readLines(
      con <- gzcon(
        url(
          "https://ftp.ebi.ac.uk/pub/databases/chebi/flat_files/database_accession.tsv.gz",
          method = "libcurl"
        )
      )
    ),
    error = function(e) conditionMessage(e),
    warning = function(w) conditionMessage(w),
    finally = if (!is.null(con)) close(con)
  )

  chebi_accession <- utils::read.delim(
    textConnection(chebi_accession_download),
    quote = "",
    stringsAsFactors = FALSE
  )

  if (nrow(chebi_accession) == 1) {
    message(chebi_accession$V1)
    return(invisible(NULL))
  }
  if (nrow(chebi_accession) == 0) {
    message("No information could be retrieved from the database, try again later!")
    return(invisible(NULL))
  }


  # Create one file with all information after cleaning individual source files

  chebi_compounds_clean <- chebi_compounds %>%
    dplyr::filter(.data$stars %in% stars) %>%
    dplyr::select(-c(
      "source",
      "name",
      "ascii_name",
      "status_id",
      "modified_on",
      "merge_type",
      "release_date"
    )) %>%
    dplyr::mutate(parent_id = as.numeric(.data$parent_id))

  chebi_accession_clean <- chebi_accession %>%
    dplyr::distinct(.data$compound_id, .data$type, .data$accession_number) %>%
    dplyr::rename(id = "compound_id") %>%
    dplyr::rename(type_accession = "type")

  chebi_chemical_data_clean <- chebi_chemical_data %>%
    dplyr::distinct(.data$compound_id, .data$formula, .data$mass, .data$charge, .data$monoisotopic_mass) %>%
    dplyr::rename(id = "compound_id")

  chebi_compounds_names_clean <- chebi_compounds %>%
    dplyr::filter(.data$stars %in% stars) %>%
    dplyr::distinct(.data$id, .data$ascii_name) %>%
    dplyr::rename(name = "ascii_name") %>%
    dplyr::mutate(type_name = "STANDARD")

  chebi_names_clean <- chebi_names %>%
    dplyr::distinct(.data$compound_id, .data$ascii_name, .data$type) %>%
    dplyr::rename(
      id = "compound_id",
      type_name = "type",
      name = "ascii_name"
    ) %>%
    dplyr::bind_rows(chebi_compounds_names_clean) %>%
    dplyr::group_by(.data$id, .data$name) %>%
    dplyr::mutate(type_name = ifelse(dplyr::n() > 1 & any(.data$type_name == "STANDARD") & .data$type_name == "SYNONYM", "STANDARD", .data$type_name)) %>%
    dplyr::ungroup() %>%
    dplyr::distinct()

  chebi <- chebi_compounds_clean %>%
    dplyr::left_join(chebi_names_clean, by = "id") %>%
    dplyr::left_join(chebi_accession_clean, by = "id", relationship = "many-to-many") %>%
    dplyr::left_join(chebi_chemical_data_clean, by = "id", relationship = "many-to-many")

  # Add info to old compound IDs

  parent_ids <- chebi %>%
    dplyr::filter(!is.na(.data$parent_id)) %>%
    dplyr::pull(var = .data$parent_id) %>%
    unique()

  parent_info <- chebi %>%
    dplyr::filter(.data$id %in% parent_ids) %>%
    dplyr::select(c(
      "id",
      "name",
      "definition",
      "type_accession",
      "accession_number",
      "formula",
      "mass",
      "charge",
      "monoisotopic_mass"
    ))

  parent_complete <- chebi %>%
    dplyr::filter(!is.na(.data$parent_id)) %>%
    dplyr::select(-c(
      "name",
      "definition",
      "type_accession",
      "accession_number",
      "formula",
      "mass",
      "charge",
      "monoisotopic_mass"
    )) %>%
    dplyr::left_join(parent_info, by = c("parent_id" = "id"), relationship = "many-to-many")

  chebi <- chebi %>%
    dplyr::filter(is.na(.data$parent_id)) %>%
    dplyr::bind_rows(parent_complete) %>%
    dplyr::arrange(.data$id) %>%
    janitor::clean_names()

  chebi
}
