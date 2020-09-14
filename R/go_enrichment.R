#' Perform gene ontology enrichment analysis
#'
#' Analyses enrichment of gene ontology terms associated with proteins in the fraction of significant proteins compared to all detected proteins. 
#' The enrichment analysis is based on the topGO R package.
#'
#' @param data A data frame that contains at least the input variables.
#' @param protein_id The name of the column containing the protein accession numbers.
#' @param is_significant The name of the column containing a logical indicating if the corresponding protein has a significantly changing peptide. 
#' The input data frame may contain peptide level information with significance information. The function is able to extract protein level information from this. 
#' @param ontology_type A character vector specifying the type of ontology that should be used. Possible values 
#' are molecular function (MF), biological process (BF), cellular component (CC).
#' @param organism_id An NCBI taxonomy identifier of an organism (TaxId). Possible inputs inlude only: "9606" (Human), "559292" (Yeast) and 
#' "83333" (E. coli). Is only necessary if gene_ontology data is not provided in \code{go_data}.
#' @param algorithm A character vector specifying the type of algorithm that should be used to group GO terms. Default is \code{"elim"}, which
#' simplifies GO terms. \code{"classic"} uses all available terms without simplification. Further options can be found here \url{https://bioconductor.org/packages/release/bioc/vignettes/topGO/inst/doc/topGO.pdf}. 
#' @param statistic A character vector specifying the type of test that should be performed. Default is \code{"fisher"}. Alternatively,
#' a Kolmogorov-Smirnov test can be used. Further options can be found here \url{https://bioconductor.org/packages/release/bioc/vignettes/topGO/inst/doc/topGO.pdf}.
#' @param go_data A list. If GO annotations should not be fetched automatically, a named list can be provided. Each element should be named
#' with an ID that corresponds to the IDs in the data. The content of each element should be a character vector with associated GO terms.
#' @param plot A logical indicating whether the result should be plotted or returned as a table.
#' 
#' @return A bar plot displaying negative log10 adjusted p-values for the top 10 enriched gene ontology terms. 
#' Bars are colored according to the direction of the enrichment. If \code{plot = FALSE}, a data frame is returned.
#'
#' @import dplyr
#' @import ggplot2
#' @import topGO
#' @importFrom janitor clean_names
#' @importFrom stringr str_replace
#' @importFrom tidyr drop_na
#' @importFrom rlang .data 
#' @importFrom magrittr %>%
#' @importFrom purrr map
#' @export
#'
#' @examples
#' \dontrun{
#' go_enrichment(
#' data,
#' protein_id = pg_protein_accessions,
#' is_significant = significant
#' ontology_type = "MF",
#' organism_id = "83333",
#' algorithm = "elim",
#' statistic = "fisher"
#' )
#' }
go_enrichment <- function(data, protein_id, is_significant, ontology_type, organism_id, algorithm = "elim", statistic = "fisher", go_data = NULL, plot = TRUE){
  . = NULL # to avoid note about no global variable binding. Usually this can be avoided with .data$ but not in nesting in complete function.
  
  if(length(unique(dplyr::pull(data, {{protein_id}}))) != nrow(data)){
    data <- data %>%
      dplyr::distinct({{protein_id}}, {{is_significant}}) %>%
      dplyr::group_by({{protein_id}}) %>%
      dplyr::mutate({{is_significant}} := ifelse(sum({{is_significant}}) > 0, TRUE, FALSE)) %>% # do this to remove accidental double annotations
      dplyr::distinct()
  }
  
  if(is.null(go_data)){
    go <- fetch_go(organism_id)
    # if the provided protein_ids are not from uniprot but SGD, they are converted to uniprot 
    if(!any(dplyr::pull(data, {{protein_id}}) %in% go$db_id)){
      if(organism_id != "559292") stop("The provided protein IDs are not UniProt or SGD IDs or the corresponding organism is wrong. 
                                     Please provide IDs in the correct format and check if you used the right organism ID.")
      annotation <- fetch_uniprot_proteome(559292, columns = c("id", "database(SGD)"), reviewed = TRUE) %>% 
        dplyr::mutate(database_sgd = stringr::str_replace(.data$database_sgd, pattern = ";", replacement = ""))
      
      go <- go %>% 
        dplyr::left_join(annotation, by = c("db_id" = "database_sgd")) %>% 
        dplyr::select(-.data$db_id) %>% 
        dplyr::rename(db_id = .data$id)
    }
    
    go_list <- go %>% 
      dplyr::select(.data$db_id, .data$go_id) %>% 
      dplyr::group_by(.data$db_id) %>% 
      dplyr::summarise(go_id = list(.data$go_id), .groups = "drop") %>% 
      tidyr::drop_na() %>% 
      split(.$db_id) %>% 
      purrr::map(.f = ~ unlist(.x$go_id))
    
  }
  
  if(!is.null(go_data)) go_list <- go_data
  
  protein_list <- factor(as.integer(dplyr::pull(data, {{is_significant}})))
  names(protein_list) <- dplyr::pull(data, {{protein_id}})
  
  go_data <- methods::new("topGOdata", 
                 ontology = ontology_type, 
                 allGenes = protein_list,
                 annot = annFUN.gene2GO,
                 gene2GO = go_list)
  
  result <- topGO::runTest(go_data, algorithm = algorithm, statistic = statistic)
  
  result_table <- topGO::GenTable(go_data, 
                                  pval = result, 
                                  orderBy = "pval", 
                                  topNodes = length(result@score), 
                                  numChar = 100) %>% 
    dplyr::mutate(elim = algorithm == "elim") %>% 
    dplyr::mutate(pval = as.numeric(.data$pval)) %>% 
    dplyr::mutate(adj_pval = ifelse(elim, .data$pval, stats::p.adjust(.data$pval, method = "BH"))) %>%
    dplyr::select(-elim) %>% 
    janitor::clean_names() %>% 
    dplyr::mutate(direction = ifelse(.data$expected < .data$significant, "Up", "Down"))
  
  if(plot == FALSE) return(result_table)
  
  enrichment_plot <- result_table %>%
    dplyr::mutate(neg_log_qval = -log10(.data$adj_pval)) %>%
    dplyr::slice(1:10) %>%
    ggplot2::ggplot(ggplot2::aes(stats::reorder(.data$term, .data$neg_log_qval), .data$neg_log_qval, fill = .data$direction)) +
    ggplot2::geom_col(col = "black", size = 1.5) +
    ggplot2::scale_fill_manual(values = c(Down = "#56B4E9", Up = "#E76145")) +
    ggplot2::scale_y_continuous(breaks = seq(0, 100, 1)) +
    ggplot2::coord_flip() +
    ggplot2::labs(title = "Gene ontology enrichment of significant proteins", y = "-log10 adjusted p-value") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        size = 20
      ),
      axis.text.x = ggplot2::element_text(
        size = 15
      ),
      axis.text.y = ggplot2::element_text(
        size = 15
      ),
      axis.title.x = ggplot2::element_text(
        size = 15
      ),
      axis.title.y = ggplot2::element_blank()
    )
  
  return(enrichment_plot)
}