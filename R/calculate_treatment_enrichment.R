#' Check treatment enrichment
#'
#' `r lifecycle::badge('deprecated')`
#' This function was deprecated due to its name changing to `calculate_treatment_enrichment()`.
#'
#' @keywords internal
#' @export
treatment_enrichment <- function(...) {
  # This function has been renamed and is therefore deprecated.
  lifecycle::deprecate_warn("0.2.0",
    "treatment_enrichment()",
    "calculate_treatment_enrichment()",
    details = "This function has been renamed."
  )

  calculate_treatment_enrichment(...)
}
#' Check treatment enrichment
#'
#' Check for an enrichment of proteins interacting with the treatment in significantly changing proteins as compared to all proteins.
#'
#' @param data A dataframe contains at least the input variables.
#' @param protein_id The name of the column containing the protein accession numbers.
#' @param is_significant The name of the column containing a logical indicating if the corresponding protein has a significantly changing peptide. The input data frame
#' may contain peptide level information with significance information. The function is able to extract protein level information from this.
#' @param binds_treatment The name of the column containing a logical indicating if the corresponding protein binds to the treatment. This information can be obtained
#' from different databases, e.g Uniprot.
#' @param treatment_name A character vector of the treatment name. It will be included in the plot title.
#' @param plot A logical indicating whether the result should be plotted or returned as a table.
#'
#' @return A bar plot displaying the percentage of all detect proteins and all significant proteins that bind to the treatment. A Fisher's exact test is performed to
#' calculate the significance of the enrichment in significant proteins compared to all proteins. The result is reported as a p-value. If \code{plot = FALSE} a contingency
#' table in long format is returned.
#'
#' @import dplyr
#' @import ggplot2
#' @importFrom tidyr pivot_wider pivot_longer complete
#' @importFrom rlang .data as_name enquo ensym !!
#' @importFrom tibble column_to_rownames
#' @importFrom magrittr %>%
#' @export
#'
#' @examples
#' \dontrun{
#' calculate_treatment_enrichment(
#'   data,
#'   protein_id = pg_protein_accessions,
#'   is_significant = significant,
#'   binds_treatment = binds_metals,
#'   treatment = "Metals"
#' )
#' }
calculate_treatment_enrichment <- function(data, protein_id, is_significant, binds_treatment, treatment_name, plot = TRUE) {
  data <- data %>%
    dplyr::distinct({{ protein_id }}, {{ is_significant }}, {{ binds_treatment }}) %>%
    dplyr::group_by({{ protein_id }}) %>%
    dplyr::mutate({{ is_significant }} := ifelse(sum({{ is_significant }}, na.rm = TRUE) > 0, TRUE, FALSE)) %>%
    dplyr::distinct()

  cont_table <- data %>%
    dplyr::group_by({{ binds_treatment }}, {{ is_significant }}) %>%
    dplyr::summarize(n = dplyr::n_distinct(!!rlang::ensym(protein_id)), .groups = "drop") %>%
    tidyr::complete({{ binds_treatment }}, {{ is_significant }}, fill = list(n = 0))


  fisher_test <- cont_table %>%
    tidyr::pivot_wider(names_from = {{ is_significant }}, values_from = .data$n) %>%
    tibble::column_to_rownames(var = rlang::as_name(rlang::enquo(binds_treatment))) %>%
    as.matrix() %>%
    stats::fisher.test()

  cont_table <- cont_table %>%
    dplyr::mutate(pval = fisher_test$p.value)

  if (plot == FALSE) {
    return(cont_table)
  }

  enrichment_plot <- cont_table %>%
    dplyr::mutate(total = sum(.data$n)) %>%
    dplyr::mutate(
      name = dplyr::case_when(
        {{ binds_treatment }} == FALSE &
          {{ is_significant }} == FALSE ~ "non_sig_non_interactor",
        {{ binds_treatment }} == FALSE &
          {{ is_significant }} == TRUE ~ "sig_non_interactor",
        {{ binds_treatment }} == TRUE &
          {{ is_significant }} == FALSE ~ "non_sig_interactor",
        {{ binds_treatment }} == TRUE &
          {{ is_significant }} == TRUE ~ "sig_interactor"
      )
    ) %>%
    dplyr::select(-c({{ binds_treatment }}, {{ is_significant }})) %>%
    tidyr::pivot_wider(names_from = .data$name, values_from = .data$n) %>%
    dplyr::mutate(total_interactor = .data$non_sig_interactor + .data$sig_interactor) %>%
    dplyr::mutate(total_sig = .data$sig_non_interactor + .data$sig_interactor) %>%
    dplyr::mutate(
      `All detected proteins` = .data$total_interactor / .data$total * 100,
      `Significant proteins` = .data$sig_interactor / .data$total_sig * 100
    ) %>%
    tidyr::pivot_longer(
      cols = c(.data$`All detected proteins`, .data$`Significant proteins`),
      names_to = "name",
      values_to = "value"
    ) %>%
    dplyr::mutate(count = ifelse(.data$name == "All detected proteins", .data$total_interactor, .data$sig_interactor)) %>%
    ggplot2::ggplot(ggplot2::aes(.data$name, .data$value)) +
    ggplot2::geom_col(fill = "cornflowerblue", col = "black", size = 1.2) +
    ggplot2::labs(title = paste0("Proteins interacting with ", treatment_name, " (p-value: ", ifelse(cont_table$pval < 0.01, formatC(cont_table$pval, format = "e", digits = 1), round(cont_table$pval, digits = 2)), ")"), x = "", y = paste("Interact with", treatment_name, "[%]")) +
    ggplot2::geom_text(aes(label = paste("n =", count)), position = position_stack(vjust = 0.5), size = 8) +
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
      axis.title.y = ggplot2::element_text(
        size = 15
      )
    )
  if (plot == TRUE) {
    return(enrichment_plot)
  }
}
