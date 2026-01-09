#' Check treatment enrichment
#'
#' `r lifecycle::badge('deprecated')`
#' This function was deprecated due to its name changing to `calculate_treatment_enrichment()`.
#'
#' @return A bar plot displaying the percentage of all detect proteins and all significant proteins
#' that bind to the treatment. A Fisher's exact test is performed to calculate the significance of
#' the enrichment in significant proteins compared to all proteins. The result is reported as a
#' p-value. If \code{plot = FALSE} a contingency table in long format is returned.
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
#' Check for an enrichment of proteins interacting with the treatment in significantly changing
#' proteins as compared to all proteins.
#'
#' @param data a data frame contains at least the input variables.
#' @param protein_id a character column in the \code{data} data frame that contains the protein
#' accession numbers.
#' @param is_significant a logical column in the \code{data} data frame that indicates if the
#' corresponding protein has a significantly changing peptide. The input data frame may contain
#' peptide level information with significance information. The function is able to extract protein
#' level information from this.
#' @param binds_treatment a logical column in the \code{data} data frame that indicates if the
#' corresponding protein binds to the treatment. This information can be obtained from different
#' databases, e.g. UniProt.
#' @param group optional, character column in the \code{data} data frame that contains information by
#' which the analysis should be grouped. The analysis will be performed separately for each of the
#' groups. This is most likely a column that labels separate comparisons of different conditions.
#' In protti the `assign_missingness()` function creates such a column automatically.
#' @param treatment_name a character value that indicates the treatment name. It will be included
#' in the plot title.
#' @param plot a logical value indicating whether the result should be plotted or returned as a
#' table.
#' @param fill_colours a character vector that specifies the fill colours of the plot.
#' @param fill_by_group a logical value that specifies if the bars in the plot should be filled by group
#' if the group argument is provided. Default is `FALSE`.
#' @param facet_n_col a numeric value that specifies the number of columns the facet plot should have if
#' a `group` column was provided.
#'
#' @return A bar plot displaying the percentage of all detected proteins and all significant proteins
#' that bind to the treatment. A Fisher's exact test is performed to calculate the significance of
#' the enrichment in significant proteins compared to all proteins. The result is reported as a
#' p-value. If \code{plot = FALSE} a contingency table in long format is returned.
#'
#' @import dplyr
#' @import ggplot2
#' @importFrom tidyr pivot_wider pivot_longer complete
#' @importFrom rlang .data as_name enquo ensym !!
#' @importFrom tibble column_to_rownames
#' @importFrom magrittr %>%
#' @importFrom purrr map2_dfr
#' @export
#'
#' @examples
#' # Create example data
#' data <- data.frame(
#'   protein_id = c(paste0("protein", 1:50)),
#'   significant = c(
#'     rep(TRUE, 20),
#'     rep(FALSE, 30)
#'   ),
#'   binds_treatment = c(
#'     rep(TRUE, 10),
#'     rep(FALSE, 10),
#'     rep(TRUE, 5),
#'     rep(FALSE, 25)
#'   ),
#'   group = c(
#'     rep("A", 5),
#'     rep("B", 15),
#'     rep("A", 15),
#'     rep("B", 15)
#'   )
#' )
#'
#' # Plot treatment enrichment
#' calculate_treatment_enrichment(
#'   data,
#'   protein_id = protein_id,
#'   is_significant = significant,
#'   binds_treatment = binds_treatment,
#'   treatment_name = "Rapamycin",
#'   plot = TRUE
#' )
#'
#' # Plot treatment enrichment by group
#' calculate_treatment_enrichment(
#'   data,
#'   protein_id = protein_id,
#'   group = group,
#'   is_significant = significant,
#'   binds_treatment = binds_treatment,
#'   treatment_name = "Rapamycin",
#'   plot = TRUE,
#'   fill_by_group = TRUE
#' )
#'
#' # Calculate treatment enrichment
#' enrichment <- calculate_treatment_enrichment(
#'   data,
#'   protein_id = protein_id,
#'   is_significant = significant,
#'   binds_treatment = binds_treatment,
#'   plot = FALSE
#' )
#'
#' enrichment
calculate_treatment_enrichment <- function(data,
                                           protein_id,
                                           is_significant,
                                           binds_treatment,
                                           group = NULL,
                                           treatment_name,
                                           plot = TRUE,
                                           fill_colours = protti::protti_colours,
                                           fill_by_group = FALSE,
                                           facet_n_col = 2) {
  # to avoid note about no global variable binding.
  . <- NULL

  group_missing <- missing(group)

  # group by the "group" argument if provided
  if (!group_missing) {
    data <- data %>%
      dplyr::ungroup() %>%
      dplyr::distinct({{ protein_id }}, {{ is_significant }}, {{ binds_treatment }}, {{ group }}) %>%
      dplyr::group_by({{ protein_id }}, {{ group }}) %>%
      dplyr::mutate({{ is_significant }} := ifelse(sum({{ is_significant }}, na.rm = TRUE) > 0,
        TRUE,
        FALSE
      )) %>%
      dplyr::ungroup() %>%
      dplyr::distinct()

    # Create contingency table
    cont_table <- data %>%
      dplyr::group_by({{ binds_treatment }}, {{ is_significant }}, {{ group }}) %>%
      dplyr::summarize(n = dplyr::n_distinct(!!rlang::ensym(protein_id)), .groups = "drop") %>%
      dplyr::group_by({{ group }}) %>%
      tidyr::complete({{ binds_treatment }}, {{ is_significant }}, fill = list(n = 0)) %>%
      dplyr::ungroup()

    fisher_test <- cont_table %>%
      split(dplyr::pull(., {{ group }})) %>%
      purrr::map2_dfr(
        .y = names(.),
        .f = ~ {
          ftest <- .x %>%
            dplyr::select(-{{ group }}) %>%
            tidyr::pivot_wider(names_from = {{ is_significant }}, values_from = .data$n) %>%
            tibble::column_to_rownames(var = rlang::as_name(rlang::enquo(binds_treatment))) %>%
            as.matrix() %>%
            stats::fisher.test()

          data.frame(pval = ftest$p.value) %>%
            dplyr::mutate({{ group }} := .y)
        }
      )

    cont_table <- cont_table %>%
      dplyr::left_join(fisher_test, by = rlang::as_name(rlang::enquo(group))) %>%
      dplyr::arrange({{ group }})
  } else {
    data <- data %>%
      dplyr::ungroup() %>%
      dplyr::distinct({{ protein_id }}, {{ is_significant }}, {{ binds_treatment }}) %>%
      dplyr::group_by({{ protein_id }}) %>%
      dplyr::mutate({{ is_significant }} := ifelse(sum({{ is_significant }}, na.rm = TRUE) > 0,
        TRUE,
        FALSE
      )) %>%
      dplyr::ungroup() %>%
      dplyr::distinct()

    # Create contingency table
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
  }

  if (plot == FALSE) {
    return(cont_table)
  }

  # Add p-value to group name for plot
  # also group to correctly calculate the total number of proteins in the next step
  if (!group_missing) {
    cont_table <- cont_table %>%
      dplyr::mutate(group_pval = paste0(
        {{ group }}, " (p-value: ",
        ifelse(.data$pval < 0.01,
          formatC(.data$pval,
            format = "e", digits = 1
          ),
          round(.data$pval, digits = 2)
        ),
        ")"
      )) %>%
      dplyr::group_by({{ group }})
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
    dplyr::mutate(count = ifelse(.data$name == "All detected proteins",
      .data$total_interactor,
      .data$sig_interactor
    )) %>%
    {
      if (fill_by_group & !group_missing) {
        ggplot2::ggplot(., ggplot2::aes(.data$name, .data$value, fill = .data$group)) +
          ggplot2::geom_col(col = "black") +
          scale_fill_manual(values = fill_colours)
      } else {
        ggplot2::ggplot(., ggplot2::aes(.data$name, .data$value)) +
          ggplot2::geom_col(fill = fill_colours[1], col = "black")
      }
    } +
    {
      if (!group_missing) {
        ggplot2::facet_wrap(~ .data$group_pval, ncol = facet_n_col)
      }
    } +
    {
      if (!group_missing) {
        ggplot2::labs(
          title = paste0(
            "Proteins interacting with ",
            treatment_name
          ),
          x = "",
          y = paste("Interact with", treatment_name, "[%]")
        )
      } else {
        ggplot2::labs(
          title = paste0(
            "Proteins interacting with ",
            treatment_name,
            " (p-value: ",
            ifelse(cont_table$pval < 0.01,
              formatC(cont_table$pval,
                format = "e", digits = 1
              ),
              round(cont_table$pval, digits = 2)
            ),
            ")"
          ),
          x = "",
          y = paste("Interact with", treatment_name, "[%]")
        )
      }
    } +
    ggplot2::geom_text(aes(label = paste("n =", count)),
      position = position_stack(vjust = 0.5),
      size = 8
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 20),
      axis.text.x = ggplot2::element_text(size = 15),
      axis.text.y = ggplot2::element_text(size = 15),
      axis.title.y = ggplot2::element_text(size = 15),
      strip.text = ggplot2::element_text(size = 15),
      legend.title = ggplot2::element_text(size = 15),
      legend.text = ggplot2::element_text(size = 13),
      strip.background = element_blank()
    )
  if (plot == TRUE) {
    return(enrichment_plot)
  }
}
