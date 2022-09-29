#' Fetch proteome data from UniProt
#'
#' Fetches proteome data from UniProt for the provided organism ID.
#'
#' @param organism_id a numeric value that specifies the NCBI taxonomy identifier (TaxId) for an
#' organism.
#' @param columns a character vector of metadata columns that should be imported from UniProt (all
#' possible columns can be found \href{https://www.uniprot.org/help/return_fields}{here}. For
#' cross-referenced database provide the database name with the prefix "xref_", e.g. \code{"xref_pdb"}).
#' Note: Not more than one or two columns should be selected otherwise the function will not be
#' able to efficiently retrieve the information. If more information is needed, \code{fetch_uniprot()}
#' can be used with the IDs retrieved by this function.
#' @param reviewed a logical value that determines if only reviewed protein entries will be retrieved.
#'
#' @return A data frame that contains all protein metadata specified in \code{columns} for the
#' organism of choice.
#' @importFrom janitor make_clean_names
#' @export
#'
#' @examples
#' \donttest{
#' head(fetch_uniprot_proteome(9606))
#' }
fetch_uniprot_proteome <-
  function(organism_id,
           columns = c("accession"),
           reviewed = TRUE) {
    if (length(organism_id) == 0) {
      stop("No valid organism ID found.")
    }
    message("Please note that some column names have changed due to UniProt updating its API! This might cause errors in your code. You can fix it by replacing the old column names with new ones.")
    if (length(columns) > 4) {
      warning(strwrap("We suggest to use the fetch_uniprot function to fetch more than four columns.",
        prefix = "\n", initial = ""
      ))
    }
    url <- "http://rest.uniprot.org/uniprotkb/stream?query="
    column_names <- janitor::make_clean_names(columns)
    collapsed_columns <- paste(columns, collapse = ",")
    reviewed <- paste0("reviewed:", ifelse(reviewed == TRUE, "true", "false"))
    organism_id <- paste0("organism_id:", organism_id)
    query_url <-
      utils::URLencode(paste0(
        url,
        reviewed,
        "+AND+",
        organism_id,
        "&format=tsv&fields=",
        collapsed_columns
      ))
    result <- try_query(query_url, progress = FALSE, show_col_types = FALSE)
    # result can either be a data.frame or it is a character string with the error message
    if (!methods::is(result, "data.frame")) {
      return(invisible(result))
    }
    colnames(result) <- column_names
    result
  }
