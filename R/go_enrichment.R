#' Perform gene ontology enrichment analysis
#'
#' Analyses enrichment of gene ontology terms associated with proteins in the fraction of significant proteins compared to all detected proteins.
#' A two-sided Fisher's exact test is performed to test significance of enrichment or depletion. GO annotations can be provided to this
#' function either through UniProt \code{go_annotations_uniprot}, through a table obtained with \code{fetch_go} in the \code{go_data} argument
#' or GO annotations are fetched automatically by the function by providing \code{ontology_type} and \code{organism_id}.
#'
#' @param data A data frame that contains at least the input variables.
#' @param protein_id The name of the column containing the protein accession numbers.
#' @param is_significant The name of the column containing a logical indicating if the corresponding protein has a significantly changing peptide.
#' The input data frame may contain peptide level information with significance information. The function is able to extract protein level information from this.
#' @param go_annotations_uniprot (Recommended) The name of the column containing gene ontology annotations obtained from UniProt using \code{fetch_uniprot}.
#' These annotations are already separated into the desired ontology type so the argument \code{ontology_type} is not required.
#' @param ontology_type Optional, A character vector specifying the type of ontology that should be used. Possible values
#' are molecular function (MF), biological process (BP), cellular component (CC). This argument is not required if GO annotations
#' are provided from UniProt in \code{go_annotations_uniprot}. It is required if annotations are provided through \code{go_data} or
#' automatically fetched.
#' @param organism_id Optional, An NCBI taxonomy identifier of an organism (TaxId). Possible inputs include only: "9606" (Human), "559292" (Yeast) and
#' "83333" (E. coli). Is only necessary if GO data is not provided either by \code{go_annotations_uniprot} or in \code{go_data}.
#' @param go_data Optional, a data frame that can be obtained with \code{fetch_go}. If you provide data yourself make sure column
#' names for protein ID (db_id) and GO ID (go_id) are the same as for data obtained with \code{fetch_go}.
#' @param plot A logical indicating whether the result should be plotted or returned as a table.
#' @param plot_cutoff A character vector indicating if the plot should contain the top 10 most significant proteins (p-value or adjusted p-value),
#' or if a significance cutoff should be used to determine the number of GO terms in the plot. This information should be provided with the
#' type first followed by the threshold separated by a space. Example are \code{plot_cutoff = "adj_pval top10"}, \code{plot_cutoff = "pval 0.05"}
#' or \code{plot_cutoff = "adj_pval 0.01"}. The threshold can be chosen freely.
#'
#' @return A bar plot displaying negative log10 adjusted p-values for the top 10 enriched or depleted gene ontology terms. Alternatively,
#' plot cutoffs can be chosen individually with the \code{plot_cutoff} argument.
#' Bars are colored according to the direction of the enrichment. If \code{plot = FALSE}, a data frame is returned. P-values
#' are adjusted with Benjamini-Hochberg.
#'
#' @import dplyr
#' @import ggplot2
#' @importFrom stringr str_replace str_split str_detect
#' @importFrom tidyr drop_na
#' @importFrom rlang .data !! ensym
#' @importFrom magrittr %>%
#' @importFrom purrr map
#' @export
#'
#' @examples
#' \dontrun{
#' go_enrichment(
#'   data,
#'   protein_id = pg_protein_accessions,
#'   is_significant = significant,
#'   go_annotations_uniprot = go_molecular_function
#' )
#' }
go_enrichment <- function(data, protein_id, is_significant, go_annotations_uniprot = NULL, ontology_type, organism_id = NULL, go_data = NULL, plot = TRUE, plot_cutoff = "adj_pval top10") {
  . <- NULL # to avoid note about no global variable binding. Usually this can be avoided with .data$ but not in nesting in complete function.
  n_sig <- NULL

  if (length(unique(dplyr::pull(data, {{ protein_id }}))) != nrow(data)) {
    data <- data %>%
      dplyr::ungroup() %>%
      dplyr::distinct({{ protein_id }}, {{ is_significant }}, {{ go_annotations_uniprot }}) %>%
      dplyr::group_by({{ protein_id }}) %>%
      dplyr::mutate({{ is_significant }} := ifelse(sum({{ is_significant }}, na.rm = TRUE) > 0, TRUE, FALSE)) %>% # do this to remove accidental double annotations
      dplyr::distinct()
  }
  if (!missing(go_annotations_uniprot)) {
    if (!stringr::str_detect(dplyr::pull(data, {{ go_annotations_uniprot }})[dplyr::pull(data, {{ go_annotations_uniprot }}) != "" & !is.na(dplyr::pull(data, {{ go_annotations_uniprot }}))][1], pattern = "\\[GO:")) {
      stop(paste("The column", rlang::as_name(rlang::enquo(go_annotations_uniprot)), "does not contain the right GO annotation format.
                 Please use the format provided by UniProt or provide data in the go_data argument."))
    }
    input <- data %>%
      dplyr::mutate({{ go_annotations_uniprot }} := stringr::str_split({{ go_annotations_uniprot }}, ";")) %>%
      tidyr::unnest({{ go_annotations_uniprot }}) %>%
      dplyr::mutate({{ go_annotations_uniprot }} := stringr::str_trim({{ go_annotations_uniprot }})) %>%
      dplyr::mutate({{ go_annotations_uniprot }} := ifelse({{ go_annotations_uniprot }} == "", NA, {{ go_annotations_uniprot }})) %>%
      dplyr::rename(go_id = {{ go_annotations_uniprot }}) %>%
      dplyr::rename(protein_id = {{ protein_id }})
  }

  if (missing(go_data) & missing(go_annotations_uniprot)) {
    if (missing(organism_id)) stop("Please provide the organism_id argument!")
    go_data <- fetch_go(organism_id)
    # if the provided protein_ids are not from uniprot but SGD, they are converted to uniprot
    if (!any(dplyr::pull(data, {{ protein_id }}) %in% go_data$db_id)) {
      if (organism_id != "559292") stop("The provided protein IDs are not UniProt or SGD IDs or the corresponding organism is wrong. 
                                       Please provide IDs in the correct format and check if you used the right organism ID.")
      annotation <- fetch_uniprot_proteome(559292, columns = c("id", "database(SGD)"), reviewed = TRUE) %>%
        dplyr::mutate(database_sgd = stringr::str_replace(.data$database_sgd, pattern = ";", replacement = ""))

      go_data <- go_data %>%
        dplyr::left_join(annotation, by = c("db_id" = "database_sgd")) %>%
        dplyr::select(-.data$db_id) %>%
        dplyr::rename(db_id = .data$id)
    }
  }

  if (!missing(go_data) & missing(go_annotations_uniprot)) {
    if (missing(ontology_type)) stop("Please provide the ontology_type argument!")
    ontology_type <- switch(ontology_type, "MF" = "F",
      "BP" = "P",
      "CC" = "C"
    )

    go_data <- go_data %>%
      dplyr::filter(.data$ontology == ontology_type) %>% # filter for right ontology
      dplyr::distinct(.data$go_id, .data$db_id) %>%
      dplyr::right_join(data, by = c("db_id" = rlang::as_name(rlang::enquo(protein_id)))) %>%
      dplyr::rename(protein_id = .data$db_id)
  }

  if (!missing(go_annotations_uniprot)) go_data <- input

  cont_table <- go_data %>%
    tidyr::drop_na(.data$go_id, {{ is_significant }}) %>%
    dplyr::group_by({{ is_significant }}) %>%
    dplyr::mutate(n_sig = dplyr::n_distinct(.data$protein_id)) %>%
    dplyr::group_by(.data$go_id, {{ is_significant }}) %>%
    dplyr::mutate(n_has_process = dplyr::n_distinct(.data$protein_id)) %>% # count number of proteins with process for sig and non-sig proteins
    dplyr::distinct(.data$go_id, {{ is_significant }}, .data$n_sig, .data$n_has_process) %>%
    dplyr::ungroup() %>%
    tidyr::complete(.data$go_id, tidyr::nesting(!!rlang::ensym(is_significant), n_sig), fill = list(n_has_process = 0))

  fisher_test <- cont_table %>%
    split(dplyr::pull(., .data$go_id)) %>%
    purrr::map(.f = ~ dplyr::select(.x, -.data$go_id) %>%
      tibble::column_to_rownames(var = rlang::as_name(enquo(is_significant))) %>%
      as.matrix() %>%
      fisher.test()) %>%
    purrr::map2_df(
      .y = names(.),
      .f = ~ tibble::tibble(
        pval = .x$p.value,
        go_id = .y
      )
    )

  result_table <- cont_table %>%
    dplyr::left_join(fisher_test, by = "go_id") %>%
    dplyr::mutate(adj_pval = stats::p.adjust(.data$pval, method = "BH")) %>%
    dplyr::group_by(.data$go_id) %>%
    dplyr::mutate(
      n_detected_proteins = sum(.data$n_sig),
      n_detected_proteins_in_process = sum(.data$n_has_process),
      n_significant_proteins = ifelse({{ is_significant }} == TRUE, .data$n_sig, NA),
      n_significant_proteins_in_process = ifelse({{ is_significant }} == TRUE, .data$n_has_process, NA)
    ) %>%
    tidyr::drop_na() %>%
    dplyr::select(-c({{ is_significant }}, .data$n_sig, .data$n_has_process)) %>%
    dplyr::mutate(n_proteins_expected = round(.data$n_significant_proteins / .data$n_detected_proteins * .data$n_detected_proteins_in_process, digits = 2)) %>%
    dplyr::mutate(direction = ifelse(.data$n_proteins_expected < .data$n_significant_proteins_in_process, "Up", "Down")) %>%
    dplyr::arrange(.data$pval)

  if (stringr::str_detect(result_table$go_id[1], pattern = "\\[GO:")) {
    result_table <- suppressWarnings(tidyr::separate(result_table, .data$go_id, c("term", "go_id"), "(?=\\[GO:)"))
  }

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

  # plot
  enrichment_plot <-
    {
      if ("term" %in% colnames(plot_input)) {
        ggplot2::ggplot(plot_input, ggplot2::aes(stats::reorder(.data$term, .data$neg_log_sig), .data$neg_log_sig, fill = .data$direction))
      } else {
        ggplot2::ggplot(plot_input, ggplot2::aes(stats::reorder(.data$go_id, .data$neg_log_sig), .data$neg_log_sig, fill = .data$direction))
      }
    } +
    ggplot2::geom_col(col = "black", size = 1.5) +
    ggplot2::scale_fill_manual(values = c(Down = "#56B4E9", Up = "#E76145")) +
    ggplot2::scale_y_continuous(breaks = seq(0, 100, 1)) +
    ggplot2::coord_flip() +
    {
      if (type == "adj_pval") ggplot2::labs(title = "Gene ontology enrichment of significant proteins", y = "-log10 adjusted p-value")
    } +
    {
      if (type == "pval") ggplot2::labs(title = "Gene ontology enrichment of significant proteins", y = "-log10 p-value")
    } +
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
