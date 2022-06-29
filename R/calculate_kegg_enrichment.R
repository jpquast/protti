#' Perform KEGG pathway enrichment analysis
#'
#' `r lifecycle::badge('deprecated')`
#' This function was deprecated due to its name changing to `calculate_kegg_enrichment()`.
#'
#' @return A bar plot displaying negative log10 adjusted p-values for the top 10 enriched pathways.
#' Bars are coloured according to the direction of the enrichment. If \code{plot = FALSE}, a data
#' frame is returned.
#' @keywords internal
#' @export
kegg_enrichment <- function(...) {
  # This function has been renamed and is therefore deprecated.
  lifecycle::deprecate_warn("0.2.0",
    "kegg_enrichment()",
    "calculate_kegg_enrichment()",
    details = "This function has been renamed."
  )

  calculate_kegg_enrichment(...)
}
#' Perform KEGG pathway enrichment analysis
#'
#' Analyses enrichment of KEGG pathways associated with proteins in the fraction of significant
#' proteins compared to all detected proteins. A Fisher's exact test is performed to test
#' significance of enrichment.
#'
#' @param data a data frame that contains at least the input variables.
#' @param protein_id a character column in the \code{data} data frame that contains the protein
#' accession numbers.
#' @param is_significant a logical column in the \code{data} data frame that indicates if the
#' corresponding protein has a significantly changing peptide. The input data frame may contain
#' peptide level information with significance information. The function is able to extract
#' protein level information from this.
#' @param pathway_id a character column in the \code{data} data frame that contains KEGG pathway
#' identifiers. These can be obtained from KEGG using \code{fetch_kegg}.
#' @param pathway_name a character column in the \code{data} data frame that contains KEGG pathway
#' names. These can be obtained from KEGG using \code{fetch_kegg}.
#' @param plot a logical value indicating whether the result should be plotted or returned as a
#' table.
#' @param plot_cutoff a character value indicating if the plot should contain the top 10 most
#' significant proteins (p-value or adjusted p-value), or if a significance cutoff should be used
#' to determine the number of GO terms in the plot. This information should be provided with the
#' type first followed by the threshold separated by a space. Example are
#' \code{plot_cutoff = "adj_pval top10"}, \code{plot_cutoff = "pval 0.05"} or
#' \code{plot_cutoff = "adj_pval 0.01"}. The threshold can be chosen freely.
#'
#' @return A bar plot displaying negative log10 adjusted p-values for the top 10 enriched pathways.
#' Bars are coloured according to the direction of the enrichment. If \code{plot = FALSE}, a data
#' frame is returned.
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
#' \donttest{
#' # Load libraries
#' library(dplyr)
#'
#' set.seed(123) # Makes example reproducible
#'
#' # Create example data
#' kegg_data <- fetch_kegg(species = "eco") 
#' 
#' if(!is.null(kegg_data)){ # only proceed if information was retrieved
#' data <- kegg_data %>%
#'   group_by(uniprot_id) %>%
#'   mutate(significant = rep(sample(
#'     x = c(TRUE, FALSE),
#'     size = 1,
#'     replace = TRUE,
#'     prob = c(0.2, 0.8)
#'   ),
#'   n = n()
#'   ))
#'
#' # Plot KEGG enrichment
#' calculate_kegg_enrichment(
#'   data,
#'   protein_id = uniprot_id,
#'   is_significant = significant,
#'   pathway_id = pathway_id,
#'   pathway_name = pathway_name,
#'   plot = TRUE,
#'   plot_cutoff = "pval 0.05"
#' )
#'
#' # Calculate KEGG enrichment
#' kegg <- calculate_kegg_enrichment(
#'   data,
#'   protein_id = uniprot_id,
#'   is_significant = significant,
#'   pathway_id = pathway_id,
#'   pathway_name = pathway_name,
#'   plot = FALSE
#' )
#'
#' head(kegg, n = 10)
#' }
#' }
calculate_kegg_enrichment <- function(data,
                                      protein_id,
                                      is_significant,
                                      pathway_id = pathway_id,
                                      pathway_name = pathway_name,
                                      plot = TRUE,
                                      plot_cutoff = "adj_pval top10") {
  . <- NULL
  n_sig <- NULL
  # to avoid note about no global variable binding. Usually this can be
  # avoided with .data$ but not in nesting in complete function.
  kegg_term <- NULL

  data <- data %>%
    dplyr::ungroup() %>%
    dplyr::distinct({{ protein_id }}, {{ is_significant }}, {{ pathway_id }}, {{ pathway_name }}) %>%
    tidyr::drop_na() %>%
    dplyr::mutate({{ pathway_name }} := stringr::str_extract({{ pathway_name }}, ".*(?=\\s\\-\\s)")) %>%
    tidyr::unite(col = "kegg_term", {{ pathway_id }}, {{ pathway_name }}, sep = ";") %>%
    dplyr::group_by({{ protein_id }}) %>%
    tidyr::nest(kegg_term = .data$kegg_term) %>%
    dplyr::distinct({{ protein_id }}, {{ is_significant }}, .data$kegg_term) %>%
    dplyr::group_by({{ protein_id }}) %>%
    dplyr::mutate({{ is_significant }} := ifelse(sum({{ is_significant }}, na.rm = TRUE) > 0,
      TRUE,
      FALSE
    )) %>%
    # do this to remove accidental double annotations
    dplyr::distinct()

  if (sum(dplyr::pull(data, {{ is_significant }})) == 0) {
    stop(strwrap("None of the significant proteins has any associated pathway
in the KEGG database. No pathway enrichment could be computed.",
      prefix = "\n", initial = ""
    ))
  }

  if (length(unique(pull(data, {{ protein_id }}))) != nrow(data)) {
    stop(strwrap("The data frame contains more rows than unique proteins. Make sure that
there are no double annotations.\nThere could be for example proteins annotated as significant
and not significant.", prefix = "\n", initial = ""))
  }
  cont_table <- data %>%
    dplyr::group_by({{ is_significant }}) %>%
    dplyr::mutate(n_sig = dplyr::n()) %>%
    tidyr::unnest(.data$kegg_term) %>%
    dplyr::mutate(kegg_term = stringr::str_trim(.data$kegg_term)) %>%
    dplyr::group_by(.data$kegg_term, {{ is_significant }}) %>%
    dplyr::mutate(n_has_pathway = dplyr::n()) %>%
    # count number of proteins with process for sig and non-sig proteins
    dplyr::distinct(.data$kegg_term, {{ is_significant }}, .data$n_sig, .data$n_has_pathway) %>%
    dplyr::ungroup() %>%
    tidyr::complete(kegg_term, tidyr::nesting(!!rlang::ensym(is_significant), n_sig), fill = list(n_has_pathway = 0))

  fisher_test <- cont_table %>%
    split(dplyr::pull(., kegg_term)) %>%
    purrr::map(.f = ~ dplyr::select(.x, -kegg_term) %>%
      tibble::column_to_rownames(var = rlang::as_name(enquo(is_significant))) %>%
      as.matrix() %>%
      fisher.test()) %>%
    purrr::map2_df(
      .y = names(.),
      .f = ~ tibble::tibble(
        pval = .x$p.value,
        kegg_term = .y
      )
    )

  result_table <- cont_table %>%
    dplyr::left_join(fisher_test, by = "kegg_term") %>%
    dplyr::mutate(adj_pval = stats::p.adjust(.data$pval, method = "BH")) %>%
    dplyr::group_by(.data$kegg_term) %>%
    dplyr::mutate(
      n_detected_proteins = sum(.data$n_sig),
      n_detected_proteins_in_pathway = sum(.data$n_has_pathway),
      n_significant_proteins = ifelse({{ is_significant }} == TRUE, .data$n_sig, NA),
      n_significant_proteins_in_pathway = ifelse({{ is_significant }} == TRUE, .data$n_has_pathway, NA)
    ) %>%
    tidyr::drop_na() %>%
    dplyr::select(-c({{ is_significant }}, .data$n_sig, .data$n_has_pathway)) %>%
    dplyr::mutate(n_proteins_expected = round(
      .data$n_significant_proteins / .data$n_detected_proteins * .data$n_detected_proteins_in_pathway,
      digits = 2
    )) %>%
    dplyr::mutate(direction = ifelse(.data$n_proteins_expected < .data$n_significant_proteins_in_pathway, "Up", "Down")) %>%
    tidyr::separate(.data$kegg_term, c("pathway_id", "pathway_name"), ";") %>%
    dplyr::arrange(.data$pval)

  if (plot == FALSE) {
    return(result_table)
  }

  # add cutoff for plot
  if (stringr::str_detect(plot_cutoff, pattern = "top10")) {
    split_cutoff <- stringr::str_split(plot_cutoff, pattern = " ", simplify = TRUE)
    type <- split_cutoff[1]
    plot_input <- result_table %>%
      dplyr::ungroup() %>%
      dplyr::mutate(neg_log_sig = -log10(!!rlang::ensym(type))) %>%
      dplyr::slice(1:10)
  } else {
    split_cutoff <- stringr::str_split(plot_cutoff, pattern = " ", simplify = TRUE)
    type <- split_cutoff[1]
    threshold <- as.numeric(split_cutoff[2])
    plot_input <- result_table %>%
      dplyr::ungroup() %>%
      dplyr::mutate(neg_log_sig = -log10(!!rlang::ensym(type))) %>%
      dplyr::filter(!!rlang::ensym(type) <= threshold)
  }

  enrichment_plot <- plot_input %>%
    ggplot2::ggplot(ggplot2::aes(stats::reorder(.data$pathway_name, .data$neg_log_sig), .data$neg_log_sig, fill = .data$direction)) +
    ggplot2::geom_col(col = "black", size = 1.5) +
    ggplot2::scale_fill_manual(values = c(Down = "#56B4E9", Up = "#E76145")) +
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
