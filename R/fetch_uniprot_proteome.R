#' Fetch proteome data from UniProt
#'
#' Fetches proteome data from UniProt for the provided organism ID.
#'
#' @param organism_id a numeric vector of one NCBI taxonomy identifier (TaxId) for and organism of choice.
#' @param columns a character vector of metadata columns that should be imported from UniProt (all possible columns can be found here: https://www.uniprot.org/help/uniprotkb_column_names). Note: Not more than one or two columns should be selected otherwise the function will not be able to efficiently retrieve the information. If more information is needed, \code{fetch_uniprot()} can be used with the IDs retrieved by this function.
#' @param reviewed a logical. If true, only reviewed protein entries will be retrieved.
#'
#' @return A data frame that contains all protein metadata specified in \code{columns} for the organism of choice.
#' @importFrom janitor make_clean_names
#' @export
#'
#' @examples
#' \donttest{
#' head(fetch_uniprot_proteome(9606))
#' }
fetch_uniprot_proteome <-
  function(organism_id,
           columns = c("id"),
           reviewed = TRUE) {
    if (!requireNamespace("httr", quietly = TRUE)) {
      stop("Package \"httr\" is needed for this function to work. Please install it.", call. = FALSE)
    }
    if (length(organism_id) == 0) {
      stop("No valid organism ID found.")
    }
    if (length(columns) > 3) {
      warning("We suggest to use the fetch_uniprot function to fetch more than three columns.")
    }
    url <- "http://www.uniprot.org/uniprot/?query="
    column_names <- janitor::make_clean_names(columns)
    collapsed_columns <- paste(columns, collapse = ",")
    reviewed <- paste0("reviewed:", ifelse(reviewed == TRUE, "yes", "no"))
    organism_id <- paste0("organism:", organism_id)
    query_url <-
      utils::URLencode(paste0(
        url,
        reviewed,
        "+AND+",
        organism_id,
        "&format=tab&columns=",
        collapsed_columns
      ))
    result <- try_query(query_url)
    if (is.null(result)) {
      return(result)
    }
    colnames(result) <- column_names
    result
  }
