#' Fetch gene ontology information from geneontology.org
#'
#' Fetches gene ontology data from geneontology.org for the provided organism ID.
#'
#' @param organism_id An NCBI taxonomy identifier of an organism (TaxId). Possible inputs inlude only: "9606" (Human), "559292" (Yeast) and
#' "83333" (E. coli).
#'
#' @return A data frame that contains gene ontology mappings to UniProt or SGD IDs. The original file is a .GAF file. A detailed description of
#' all columns can be found here: http://geneontology.org/docs/go-annotation-file-gaf-format-2.1/
#' @export
#'
#' @examples
#' \donttest{
#' fetch_go("9606")
#' }
fetch_go <- function(organism_id) {
  organism_id <- match.arg(organism_id, c("9606", "559292", "83333"))

  organism_url <- switch(organism_id, "9606" = "http://current.geneontology.org/annotations/goa_human.gaf.gz",
    "559292" = "http://current.geneontology.org/annotations/sgd.gaf.gz",
    "83333" = "http://current.geneontology.org/annotations/ecocyc.gaf.gz"
  )
  go_download <- readLines(gzcon(url(organism_url)))
  go <- utils::read.delim(textConnection(go_download), quote = "", stringsAsFactors = FALSE, comment.char = "!", header = FALSE)
  colnames(go) <- c(
    "db", "db_id", "symbol", "qualifier", "go_id", "db_reference", "evidence", "with_from", "ontology", "name", "synonyme",
    "type", "taxon", "date", "assigned_by", "annotation_extension", "gene_product_form_id"
  )
  return(go)
}
