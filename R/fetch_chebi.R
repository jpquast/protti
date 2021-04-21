#' Fetch ChEBI database information
#'
#' Fetches all information from the ChEBI database.
#'
#' @param relation a logical, if TRUE, ChEBI Ontology data will be returned instead the main compound data. This data can be used to check the relations of ChEBI ID's to eachother.
#'
#' @return A data frame that contains all information about each molecule in the ChEBI database. Only "3-star" observations are included in the result. These are entries manually annotated by the ChEBI curator team.
#' @import dplyr
#' @importFrom tidyr pivot_wider unnest
#' @importFrom janitor clean_names
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @export
#'
#' @examples
#' \donttest{
#' fetch_chebi()
#' }
fetch_chebi <- function(relation = FALSE) {
  if (!requireNamespace("httr", quietly = TRUE)) {
    stop("Package \"httr\" is needed for this function to work. Please install it.", call. = FALSE)
  }
  # Retrieve relational information
  if (relation == TRUE) {
    chebi_relation_result <- httr::GET("ftp://ftp.ebi.ac.uk/pub/databases/chebi/Flat_file_tab_delimited/relation.tsv", config = httr::config(connecttimeout = 60))
    chebi_relation <- suppressMessages(httr::content(chebi_relation_result, type = "text/tab-separated-values"))
    
    chebi_relation_clean <- chebi_relation %>%
      dplyr::filter(.data$STATUS == "C") %>%
      dplyr::select(-c(.data$ID, .data$STATUS)) %>%
      dplyr::rename(incoming = .data$FINAL_ID) %>%
      dplyr::rename(ID = .data$INIT_ID) %>%
      janitor::clean_names()

    return(chebi_relation_clean)
  }

  # Download compound data
  chebi_chemical_data_result <- httr::GET("ftp://ftp.ebi.ac.uk/pub/databases/chebi/Flat_file_tab_delimited/chemical_data.tsv", config = httr::config(connecttimeout = 60)) 
  chebi_chemical_data <- suppressMessages(httr::content(chebi_chemical_data_result, type = "text/tab-separated-values"))

  chebi_compounds_download <- readLines(gzcon(url("ftp://ftp.ebi.ac.uk/pub/databases/chebi/Flat_file_tab_delimited/compounds.tsv.gz")))
  chebi_compounds <- utils::read.delim(textConnection(chebi_compounds_download), quote = "", stringsAsFactors = FALSE)

  chebi_names_download <- readLines(gzcon(url("ftp://ftp.ebi.ac.uk/pub/databases/chebi/Flat_file_tab_delimited/names.tsv.gz")))
  chebi_names <- utils::read.delim(textConnection(chebi_names_download), quote = "", stringsAsFactors = FALSE)

  chebi_accession_result <- httr::GET("ftp://ftp.ebi.ac.uk/pub/databases/chebi/Flat_file_tab_delimited/database_accession.tsv", config = httr::config(connecttimeout = 60))
  chebi_accession <- suppressMessages(httr::content(chebi_accession_result, type = "text/tab-separated-values"))

  # Create one file with all information after cleaning individual source files

  chebi_compounds_clean <- chebi_compounds %>%
    dplyr::filter(.data$STAR == 3) %>%
    dplyr::select(-c(.data$SOURCE, .data$NAME, .data$STATUS, .data$MODIFIED_ON, .data$CREATED_BY)) %>%
    dplyr::na_if("null") %>%
    dplyr::mutate(PARENT_ID = as.numeric(.data$PARENT_ID))


  chebi_accession_clean <- chebi_accession %>%
    dplyr::distinct(.data$COMPOUND_ID, .data$TYPE, .data$ACCESSION_NUMBER) %>%
    dplyr::rename(ID = .data$COMPOUND_ID) %>%
    dplyr::rename(TYPE_ACCESSION = .data$TYPE)

  chebi_chemical_data_clean <- chebi_chemical_data %>%
    dplyr::distinct(.data$COMPOUND_ID, .data$TYPE, .data$CHEMICAL_DATA) %>%
    dplyr::rename(ID = .data$COMPOUND_ID) %>%
    tidyr::pivot_wider(names_from = .data$TYPE, values_from = .data$CHEMICAL_DATA, values_fn = list) %>%
    tidyr::unnest(cols = c(.data$FORMULA, .data$MASS, .data$CHARGE, .data$`MONOISOTOPIC MASS`))

  chebi_compounds_names_clean <- chebi_compounds %>%
    dplyr::filter(.data$STAR == 3) %>%
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
    dplyr::select(c(.data$ID, .data$NAME, .data$TYPE_NAME, .data$DEFINITION, .data$TYPE_ACCESSION, .data$ACCESSION_NUMBER, .data$FORMULA, .data$MASS, .data$CHARGE, .data$`MONOISOTOPIC MASS`))

  parent_complete <- chebi %>%
    dplyr::filter(!is.na(.data$PARENT_ID)) %>%
    dplyr::select(-c(.data$NAME, .data$TYPE_NAME, .data$DEFINITION, .data$TYPE_ACCESSION, .data$ACCESSION_NUMBER, .data$FORMULA, .data$MASS, .data$CHARGE, .data$`MONOISOTOPIC MASS`)) %>%
    dplyr::left_join(parent_info, by = c("PARENT_ID" = "ID"))

  chebi <- chebi %>%
    dplyr::filter(is.na(.data$PARENT_ID)) %>%
    dplyr::bind_rows(parent_complete) %>%
    dplyr::arrange(.data$ID) %>%
    janitor::clean_names()

  chebi
}
