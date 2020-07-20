#' Fetch ChEBI database information
#'
#' Fetches all information from the ChEBI database.
#'
#' @return A data frame that contains all information about each molecule in the ChEBI database. Only "3-star" observations are included in the result. These are entries manually annotated by the ChEBI curator team.
#' @import dplyr
#' @importFrom tidyr pivot_wider
#' @importFrom tidyr unnest
#' @importFrom janitor clean_names
#' @importFrom RCurl getURL
#' @importFrom data.table fread
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @export
#'
#' @examples
#' \dontrun{
#' fetch_chebi()
#' }
fetch_chebi <- function() {
  chebi_chemical_data_download <- RCurl::getURL("ftp://ftp.ebi.ac.uk/pub/databases/chebi/Flat_file_tab_delimited/chemical_data.tsv")
  chebi_chemical_data <- data.table::fread(chebi_chemical_data_download)
  
  chebi_compounds_download <- readLines(gzcon(url("ftp://ftp.ebi.ac.uk/pub/databases/chebi/Flat_file_tab_delimited/compounds.tsv.gz")))
  chebi_compounds <- utils::read.delim(textConnection(chebi_compounds_download), quote = "")
  
  chebi_names_download <- readLines(gzcon(url("ftp://ftp.ebi.ac.uk/pub/databases/chebi/Flat_file_tab_delimited/names.tsv.gz")))
  chebi_names <- utils::read.delim(textConnection(chebi_names_download), quote = "")
  
  chebi_accession_download <- RCurl::getURL("ftp://ftp.ebi.ac.uk/pub/databases/chebi/Flat_file_tab_delimited/database_accession.tsv")
  chebi_accession <- data.table::fread(chebi_accession_download)
  
 # Create one file with all information after cleaning individual source files
  
  chebi_compounds_clean <- chebi_compounds %>%
    dplyr::filter(.data$STAR == 3) %>%
    dplyr::select(-c(.data$SOURCE, .data$NAME, .data$STATUS, .data$MODIFIED_ON, .data$CREATED_BY)) %>%
    dplyr::na_if("null") %>%
    dplyr::mutate(PARENT_ID = as.numeric(.data$PARENT_ID))
  
  
  chebi_accession_clean <- chebi_accession %>%
    dplyr::distinct(.data$COMPOUND_ID, .data$TYPE, .data$ACCESSION_NUMBER) %>%
    dplyr::rename(ID = .data$COMPOUND_ID)
    
  
  chebi_chemical_data_clean <- chebi_chemical_data %>%
    dplyr::distinct(.data$COMPOUND_ID, .data$TYPE, .data$CHEMICAL_DATA) %>%
    dplyr::rename(ID = .data$COMPOUND_ID) %>%
    tidyr::pivot_wider(names_from = .data$TYPE, values_from = .data$CHEMICAL_DATA, values_fn = list) %>%
    tidyr::unnest(cols = c(.data$FORMULA, .data$MASS, .data$CHARGE, .data$`MONOISOTOPIC MASS`))
  
  chebi_names_clean <- chebi_names %>%
    dplyr::distinct(.data$COMPOUND_ID, .data$NAME) %>%
    dplyr::rename(ID = .data$COMPOUND_ID)
  
  chebi <- chebi_compounds_clean %>%
    dplyr::left_join(chebi_names_clean, by = "ID") %>%
    dplyr::left_join(chebi_accession_clean, by = "ID") %>%
    dplyr::left_join(chebi_chemical_data_clean, by = "ID")
  
  # Add info to old compound IDs 
  
  parent_ids <- chebi %>%
    dplyr::filter(.data$PARENT_ID != "null") %>%
    dplyr::pull(var = .data$PARENT_ID) %>%
    unique()
  
  parent_info <- chebi %>%
    dplyr::filter(.data$ID %in% parent_ids) %>%
    dplyr::select(c(.data$ID, .data$NAME, .data$DEFINITION, .data$TYPE, .data$ACCESSION_NUMBER, .data$FORMULA, .data$MASS, .data$CHARGE, .data$`MONOISOTOPIC MASS`))
  
  parent_complete <- chebi %>%
    dplyr::filter(!is.na(.data$PARENT_ID)) %>%
    dplyr::select(-c(.data$NAME, .data$DEFINITION, .data$TYPE, .data$ACCESSION_NUMBER, .data$FORMULA, .data$MASS, .data$CHARGE, .data$`MONOISOTOPIC MASS`)) %>%
    dplyr::left_join(parent_info, by = c("PARENT_ID" = "ID")) 
  
  chebi <- chebi %>%
    dplyr::filter(is.na(.data$PARENT_ID)) %>%
    dplyr::bind_rows(parent_complete) %>%
    dplyr::arrange(.data$ID) %>%
    janitor::clean_names()
  
  chebi
}