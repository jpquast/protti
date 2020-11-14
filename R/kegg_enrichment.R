#' Perform KEGG pathway enrichment analysis
#'
#' Analyses enrichment of KEGG pathways associated with proteins in the fraction of significant proteins compared to all detected proteins. A Fisher's 
#' exact test is performed to test significance of enrichment. 
#'
#' @param data A data frame that contains at least the input variables.
#' @param protein_id The name of the column containing the protein accession numbers.
#' @param is_significant  The name of the column containing a logical indicating if the corresponding protein has a significantly changing peptide. 
#' The input data frame may contain peptide level information with significance information. The function is able to extract protein level information from this. 
#' @param pathway_id The name of the column containing KEGG pathway identifiers. These can be obtained from KEGG using \code{fetch_kegg}. 
#' @param pathway_name The name of the column containing KEGG pathway names. These can be obtained from KEGG using \code{fetch_kegg}. 
#' @param plot A logical indicating whether the result should be plotted or returned as a table.
#' 
#' @return A bar plot displaying negative log10 adjusted p-values for the top 10 enriched pathways. Bars are coloured according to the direction of the 
#' enrichment. If \code{plot = FALSE}, a data frame is returned.
#'
#' @import dplyr
#' @import ggplot2
#' @importFrom stringr str_split str_extract
#' @importFrom tidyr unnest nesting complete pivot_wider separate
#' @importFrom rlang .data := !! ensym as_name enquo
#' @importFrom tibble column_to_rownames
#' @importFrom magrittr %>%
#' @importFrom purrr map map2_df
#' @export
#'
#' @examples
#' \dontrun{
#' kegg_enrichment(
#' data,
#' protein_id = pg_protein_accessions,
#' is_significant = significant,
#' pathway_id = pathway_id,
#' pathway_name = pathway_name
#' )
#' }
kegg_enrichment <- function(data, protein_id, is_significant, pathway_id, pathway_name, plot = TRUE){
  . = NULL
  n_sig = NULL
  kegg_term = NULL # to avoid node about no global variable binding. Usually this can be avoided with .data$ but not in nesting in complete function.
  
  data <- data %>%
    tidyr::drop_na() %>%
    dplyr::mutate({{pathway_name}} := stringr::str_extract({{pathway_name}}, ".*(?=\\s\\-\\s)")) %>%
    tidyr::unite(col = "kegg_term", {{pathway_id}}, {{pathway_name}}, sep = ";") %>%
    dplyr::group_by({{protein_id}}) %>%
    tidyr::nest(kegg_term = .data$kegg_term) %>%
    dplyr::distinct({{protein_id}}, {{is_significant}}, .data$kegg_term) %>%
    dplyr::group_by({{protein_id}}) %>%
    dplyr::mutate({{is_significant}} := ifelse(sum({{is_significant}}) > 0, TRUE, FALSE)) %>% # do this to remove accidental double annotations
    dplyr::distinct()
  
  if(length(unique(pull(data, {{protein_id}}))) != nrow(data)) stop("The data frame contains more rows than unique proteins. Make sure that there
                                                                    are no double annotations.")
   cont_table <- data %>%
     dplyr::group_by({{is_significant}}) %>%
     dplyr::mutate(n_sig = dplyr::n()) %>%
     tidyr::unnest(.data$kegg_term) %>%
     dplyr::mutate(kegg_term = stringr::str_trim(.data$kegg_term)) %>%
     dplyr::group_by(.data$kegg_term, {{is_significant}}) %>%
     dplyr::mutate(n_has_pathway = dplyr::n()) %>% # count number of proteins with process for sig and non-sig proteins
     dplyr::distinct(.data$kegg_term, {{is_significant}}, .data$n_sig, .data$n_has_pathway) %>%
     dplyr::ungroup() %>%
     tidyr::complete(kegg_term, tidyr::nesting(!!rlang::ensym(is_significant), n_sig), fill = list(n_has_pathway = 0))
   
  fisher_test <- cont_table  %>%
    split(dplyr::pull(., kegg_term)) %>%
    purrr::map(.f = ~ dplyr::select(.x, -kegg_term) %>%
                 tibble::column_to_rownames(var = rlang::as_name(enquo(is_significant))) %>%
                 as.matrix() %>%
                 fisher.test()
    ) %>%
    purrr::map2_df(.y = names(.),
                   .f = ~ tibble::tibble(pval = .x$p.value,
                                         kegg_term = .y)
    )

  result_table <- cont_table %>%
    dplyr::left_join(fisher_test, by = "kegg_term") %>%
    dplyr::mutate(adj_pval = stats::p.adjust(.data$pval, method = "BH")) %>%
    dplyr::group_by(.data$kegg_term) %>%
    dplyr::mutate(n_detected_proteins = sum(.data$n_sig),
                  n_detected_proteins_in_pathway = sum(.data$n_has_pathway),
                  n_significant_proteins = ifelse({{is_significant}} == TRUE, .data$n_sig, NA),
                  n_significant_proteins_in_pathway = ifelse({{is_significant}} == TRUE, .data$n_has_pathway, NA)
    ) %>%
    tidyr::drop_na() %>%
    dplyr::select(-c({{is_significant}}, .data$n_sig, .data$n_has_pathway)) %>%
    dplyr::mutate(n_proteins_expected = round(.data$n_significant_proteins/.data$n_detected_proteins * .data$n_detected_proteins_in_pathway, digits = 2)) %>%
    dplyr::mutate(direction = ifelse(.data$n_proteins_expected < .data$n_significant_proteins_in_pathway, "Up", "Down")) %>%
    tidyr::separate(.data$kegg_term, c("pathway_id", "pathway_name"), ";") %>%
    dplyr::arrange(.data$pval)

   if(plot == FALSE) return(result_table)

  enrichment_plot <- result_table %>%
    dplyr::mutate(neg_log_qval = -log10(.data$adj_pval)) %>%
    dplyr::slice(1:10) %>%
    ggplot2::ggplot(ggplot2::aes(stats::reorder(.data$pathway_name, .data$neg_log_qval), .data$neg_log_qval, fill = .data$direction)) +
    ggplot2::geom_col(col = "black", size = 1.5) +
    ggplot2::scale_fill_manual(values = c("#56B4E9", "#E76145")) +
    ggplot2::scale_y_continuous(breaks = seq(0, 100, 2)) +
    ggplot2::coord_flip() +
    ggplot2::labs(title = "KEGG pathway enrichment of significant proteins", y = "-log10 adjusted p-value") +
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
