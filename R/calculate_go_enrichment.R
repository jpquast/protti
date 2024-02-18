#' Perform gene ontology enrichment analysis
#'
#' `r lifecycle::badge('deprecated')`
#' This function was deprecated due to its name changing to `calculate_go_enrichment()`.
#'
#' @return A bar plot displaying negative log10 adjusted p-values for the top 10 enriched or
#' depleted gene ontology terms. Alternatively, plot cutoffs can be chosen individually with the
#' \code{plot_cutoff} argument. Bars are colored according to the direction of the enrichment
#' (enriched or deenriched). If \code{plot = FALSE}, a data frame is returned. P-values are
#' adjusted with Benjamini-Hochberg.
#' @keywords internal
#' @export
go_enrichment <- function(...) {
  # This function has been renamed and is therefore deprecated.
  lifecycle::deprecate_warn("0.2.0",
    "go_enrichment()",
    "calculate_go_enrichment()",
    details = "This function has been renamed."
  )

  calculate_go_enrichment(...)
}
#' Perform gene ontology enrichment analysis
#'
#' Analyses enrichment of gene ontology terms associated with proteins in the fraction of
#' significant proteins compared to all detected proteins. A two-sided Fisher's exact test is
#' performed to test significance of enrichment or depletion. GO annotations can be provided to
#' this function either through UniProt \code{go_annotations_uniprot}, through a table obtained
#' with \code{fetch_go} in the \code{go_data} argument or GO annotations are fetched automatically
#' by the function by providing \code{ontology_type} and \code{organism_id}.
#'
#' @param data a data frame that contains at least the input variables.
#' @param protein_id a character column in the \code{data} data frame that contains the protein
#' accession numbers.
#' @param is_significant a logical column in the \code{data} data frame that indicates if the
#' corresponding protein has a significantly changing peptide. The input data frame may contain
#' peptide level information with significance information. The function is able to extract
#' protein level information from this.
#' @param group optional, character column in the \code{data} data frame that contains information by
#' which the analysis should be grouped. The analysis will be performed separately for each of the
#' groups. This is most likely a column that labels separate comparisons of different conditions.
#' In protti the `assign_missingness()` function creates such a column automatically.
#' @param y_axis_free a logical value that specifies if the y-axis of the plot should be "free"
#' for each facet if a grouping variable is provided. Default is `TRUE`. If `FALSE` is selected
#' it is easier to compare GO categories directly with each other.
#' @param facet_n_col a numeric value that specifies the number of columns the faceted plot should have
#' if a column name is provided to group. The default is 2.
#' @param go_annotations_uniprot recommended, a character column in the \code{data} data frame
#' that contains gene ontology annotations obtained from UniProt using \code{fetch_uniprot}.
#' These annotations are already separated into the desired ontology type so the argument
#' \code{ontology_type} is not required.
#' @param ontology_type optional, character value specifying the type of ontology that should
#' be used. Possible values are molecular function (MF), biological process (BP), cellular component
#' (CC). This argument is not required if GO annotations are provided from UniProt in
#' \code{go_annotations_uniprot}. It is required if annotations are provided through \code{go_data} or
#' automatically fetched.
#' @param organism_id optional, character value specifying an NCBI taxonomy identifier of an
#' organism (TaxId). Possible inputs include only: "9606" (Human), "559292" (Yeast) and "83333"
#' (E. coli). Is only necessary if GO data is not provided either by \code{go_annotations_uniprot}
#' or in \code{go_data}.
#' @param go_data Optional, a data frame that can be obtained with `fetch_go()`. If you provide
#' data not obtained with `fetch_go()` make sure column names for protein ID (`db_id`) and GO ID
#' (`go_id`) are the same as for data obtained with `fetch_go()`.
#' @param plot a logical argument indicating whether the result should be plotted or returned as a table.
#' @param plot_title a character value that specifies the title of the plot. The default is "Gene ontology
#' enrichment of significant proteins".
#' @param label a logical argument indicating whether labels should be added to the plot.
#' Default is TRUE.
#' @param enrichment_type a character argument that is either "all", "enriched" or "deenriched". This
#' determines if the enrichment analysis should be performed in order to check for both enrichemnt and
#' deenrichemnt or only one of the two. This affects the statistics performed and therefore also the displayed
#' plot.
#' @param min_n_detected_proteins_in_process is a numeric argument that specifies the minimal number
#' of detected proteins a GO term needs to have to be displayed in the plot. The default is 1, meaning
#' no filtering of the plotted data is performed. This argument does not affect any computations or
#' the returned data if `plot = FALSE`. This argument is useful in order to remove terms that were only
#' detected in for example 1 protein. Even though these terms are sometimes significant, they are not
#' really relevant.
#' @param plot_cutoff a character value indicating if the plot should contain the top 10 most
#' significant proteins (p-value or adjusted p-value), or if a significance cutoff should be used
#' to determine the number of GO terms in the plot. This information should be provided with the
#' type first followed by the threshold separated by a space. Example are
#' `plot_cutoff = "adj_pval top10"`, `plot_cutoff = "pval 0.05"` or
#' `plot_cutoff = "adj_pval 0.01"`. The threshold can be chosen freely. The default value is
#' `"adj_pval top10"`.
#'
#' @return A bar plot displaying negative log10 adjusted p-values for the top 10 enriched or
#' depleted gene ontology terms. Alternatively, plot cutoffs can be chosen individually with the
#' `plot_cutoff` argument. Bars are colored according to the direction of the enrichment (enriched
#' or deenriched). If `plot = FALSE`, a data frame is returned. P-values are adjusted with
#' Benjamini-Hochberg.
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
#' \donttest{
#' # Load libraries
#' library(dplyr)
#' library(stringr)
#'
#' # Create example data
#' # Contains artificial de-enrichment for ribosomes.
#' uniprot_go_data <- fetch_uniprot_proteome(
#'   organism_id = 83333,
#'   columns = c(
#'     "accession",
#'     "go_f"
#'   )
#' )
#'
#' if (!is(uniprot_go_data, "character")) {
#'   data <- uniprot_go_data %>%
#'     mutate(significant = c(
#'       rep(TRUE, 1000),
#'       rep(FALSE, n() - 1000)
#'     )) %>%
#'     mutate(significant = ifelse(
#'       str_detect(
#'         go_f,
#'         pattern = "ribosome"
#'       ),
#'       FALSE,
#'       significant
#'     )) %>%
#'     mutate(group = c(rep("A", 500),
#'                      rep("B", 500),
#'                      rep("A", (n() - 1000)/2),
#'                      rep("B", round((n() - 1000)/2))))
#'
#'   # Plot gene ontology enrichment
#'   calculate_go_enrichment(
#'     data,
#'     protein_id = accession,
#'     go_annotations_uniprot = go_f,
#'     is_significant = significant,
#'     plot = TRUE,
#'     plot_cutoff = "pval 0.01"
#'   )
#'
#'  # Plot gene ontology enrichment with group
#'  calculate_go_enrichment(
#'    data,
#'    protein_id = accession,
#'    go_annotations_uniprot = go_f,
#'    is_significant = significant,
#'    group = group,
#'    facet_n_col = 1,
#'    plot = TRUE,
#'    plot_cutoff = "pval 0.01"
#'  )
#'
#'   # Calculate gene ontology enrichment
#'   go_enrichment <- calculate_go_enrichment(
#'     data,
#'     protein_id = accession,
#'     go_annotations_uniprot = go_f,
#'     is_significant = significant,
#'     plot = FALSE,
#'   )
#'
#'   head(go_enrichment, n = 10)
#' }
#' }
calculate_go_enrichment <- function(data,
                                    protein_id,
                                    is_significant,
                                    group = NULL,
                                    y_axis_free = TRUE,
                                    facet_n_col = 2,
                                    go_annotations_uniprot = NULL,
                                    ontology_type,
                                    organism_id = NULL,
                                    go_data = NULL,
                                    plot = TRUE,
                                    plot_title = "Gene ontology enrichment of significant proteins",
                                    label = TRUE,
                                    enrichment_type = "all",
                                    min_n_detected_proteins_in_process = 1,
                                    plot_cutoff = "adj_pval top10") {
  # to avoid note about no global variable binding. Usually this can be avoided with
  # .data$ but not in nesting in complete function.
  . <- NULL
  n_sig <- NULL
  group_missing <- missing(group)
  alternative_test <- switch(enrichment_type,
                             "all" = "two.sided",
                             "enriched" = "greater",
                             "deenriched" = "less",
                             "none")

  if (alternative_test == "none") stop("Please provide a valid enrichment_type! Can be 'all', 'enriched', 'deenriched'.")

  if (length(unique(dplyr::pull(data, {{ protein_id }}))) != nrow(data)) {
    data <- data %>%
      dplyr::ungroup() %>%
    {
      # conditional grouping
      if (!group_missing) dplyr::distinct(., {{ protein_id }}, {{ is_significant }}, {{ go_annotations_uniprot }}, {{ group }}) %>%
        dplyr::group_by({{ protein_id }}, {{ group }}) else
        dplyr::distinct(., {{ protein_id }}, {{ is_significant }}, {{ go_annotations_uniprot }}) %>%
        dplyr::group_by({{ protein_id }})
    } %>%
      dplyr::mutate({{ is_significant }} := ifelse(sum({{ is_significant }}, na.rm = TRUE) > 0, TRUE, FALSE)) %>%
      # do this to remove accidental double annotations
      dplyr::ungroup() %>%
      dplyr::distinct()
  }

  if (!missing(go_annotations_uniprot)) {
    if (!stringr::str_detect(
      dplyr::pull(
        data,
        {{ go_annotations_uniprot }}
      )[dplyr::pull(
        data,
        {{ go_annotations_uniprot }}
      ) != "" &
        !is.na(dplyr::pull(
          data,
          {{ go_annotations_uniprot }}
        ))][1],
      pattern = "\\[GO:"
    )) {
      stop(strwrap(paste("The column", rlang::as_name(rlang::enquo(go_annotations_uniprot)), "does not
contain the right GO annotation format. Please use the format provided by UniProt or provide data in
the go_data argument."), prefix = "\n", initial = ""))
    }
    input <- data %>%
      dplyr::mutate({{ go_annotations_uniprot }} := stringr::str_split({{ go_annotations_uniprot }}, ";")) %>%
      tidyr::unnest({{ go_annotations_uniprot }}) %>%
      dplyr::mutate({{ go_annotations_uniprot }} := stringr::str_trim({{ go_annotations_uniprot }})) %>%
      dplyr::mutate({{ go_annotations_uniprot }} := ifelse({{ go_annotations_uniprot }} == "",
        NA,
        {{ go_annotations_uniprot }}
      )) %>%
      dplyr::rename(go_id = {{ go_annotations_uniprot }}) %>%
      dplyr::rename(protein_id = {{ protein_id }})
  }

  if (missing(go_data) & missing(go_annotations_uniprot)) {
    if (missing(organism_id)) stop("Please provide the organism_id argument!")
    go_data <- fetch_go(organism_id)
    # if the provided protein_ids are not from uniprot but SGD, they are converted to uniprot
    if (!any(dplyr::pull(data, {{ protein_id }}) %in% go_data$db_id)) {
      if (organism_id != "559292") stop(strwrap("The provided protein IDs are not UniProt or SGD
IDs or the corresponding organism is wrong. Please provide IDs in the correct format and check
if you used the right organism ID.", prefix = "\n", initial = ""))
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
    ontology_type <- switch(ontology_type,
      "MF" = "F",
      "BP" = "P",
      "CC" = "C"
    )

    go_data <- go_data %>%
      dplyr::filter(.data$ontology == ontology_type) %>%
      # filter for right ontology
      dplyr::distinct(.data$go_id, .data$db_id) %>%
      dplyr::right_join(data, by = c("db_id" = rlang::as_name(rlang::enquo(protein_id)))) %>%
      dplyr::rename(protein_id = .data$db_id)
  }

  if (!missing(go_annotations_uniprot)) go_data <- input

  # group argument is not missing
  cont_table <- go_data %>%
    tidyr::drop_na(.data$go_id, {{ is_significant }}) %>%
    {   # group argument is missing
      if (group_missing) dplyr::group_by(., {{ is_significant }}) else
        # group argument is not missing
        dplyr::group_by(., {{ is_significant }}, {{ group }})
    } %>%
    dplyr::mutate(n_sig = dplyr::n_distinct(.data$protein_id)) %>%
    dplyr::group_by(.data$go_id, .add = TRUE) %>%
    dplyr::mutate(n_has_process = dplyr::n_distinct(.data$protein_id)) %>%
    # count number of proteins with process for sig and non-sig proteins
    {   # group argument is missing
      if (group_missing) dplyr::distinct(., .data$go_id, {{ is_significant }}, .data$n_sig, .data$n_has_process) %>%
        dplyr::ungroup() else
        # group argument is not missing
          dplyr::distinct(., .data$go_id, {{ is_significant }}, .data$n_sig, .data$n_has_process, {{ group }}) %>%
        dplyr::group_by({{ group }})
    } %>%
    tidyr::complete(.data$go_id, tidyr::nesting(!!rlang::ensym(is_significant), n_sig), fill = list(n_has_process = 0)) %>%
    dplyr::ungroup()


  if (group_missing) {
    fisher_test <- cont_table %>%
      split(dplyr::pull(., .data$go_id)) %>%
      purrr::map(.f = ~ dplyr::select(.x, -.data$go_id) %>%
        tibble::column_to_rownames(var = rlang::as_name(enquo(is_significant))) %>%
        as.matrix() %>%
        fisher.test(alternative = alternative_test)) %>%
      purrr::map2_df(
        .y = names(.),
        .f = ~ tibble::tibble(
          pval = .x$p.value,
          go_id = .y
        )
      )

  } else {
    fisher_test <- cont_table %>%
      split(dplyr::pull(., {{ group }})) %>%
      purrr::map2_dfr(
        .y = names(.),
        .f = ~ {
          .x %>%
            dplyr::select(-{{ group }}) %>%
            split(dplyr::pull(., .data$go_id)) %>%
            purrr::map(.f = ~ {
              dplyr::select(.x, -c(.data$go_id)) %>%
                tibble::column_to_rownames(var = rlang::as_name(enquo(is_significant))) %>%
                as.matrix() %>%
                fisher.test(alternative = alternative_test)
            }) %>%
            purrr::map2_dfr(
              .y = names(.),
              .f = ~ tibble::tibble(
                pval = .x$p.value,
                go_id = .y
              )
            ) %>%
            mutate({{ group }} := .y)
        }
      )
  }

  result_table <- cont_table %>%
    {   # group argument is missing
      if(group_missing) dplyr::left_join(., fisher_test, by = "go_id") else
        # group argument is not missing
        dplyr::left_join(., fisher_test, by = c("go_id", rlang::as_name(enquo(group)))) %>%
        dplyr::group_by({{ group }})
    } %>%
    dplyr::mutate(adj_pval = stats::p.adjust(.data$pval, method = "BH")) %>%
    dplyr::group_by(.data$go_id, .add = TRUE) %>%
    dplyr::mutate(
      n_detected_proteins = sum(.data$n_sig),
      n_detected_proteins_in_process = sum(.data$n_has_process),
      n_significant_proteins = ifelse({{ is_significant }} == TRUE, .data$n_sig, NA),
      n_significant_proteins_in_process = ifelse({{ is_significant }} == TRUE, .data$n_has_process, NA)
    ) %>%
    tidyr::drop_na() %>%
    dplyr::select(-c({{ is_significant }}, .data$n_sig, .data$n_has_process)) %>%
    dplyr::mutate(n_proteins_expected = round(
      .data$n_significant_proteins / .data$n_detected_proteins * .data$n_detected_proteins_in_process,
      digits = 2
    )) %>%
    dplyr::mutate(enrichment_type = ifelse(.data$n_proteins_expected < .data$n_significant_proteins_in_process, "Enriched", "Deenriched")) %>%
    dplyr::arrange(.data$pval)

  if (stringr::str_detect(result_table$go_id[1], pattern = "\\[GO:")) {
    result_table <- suppressWarnings(tidyr::separate(result_table, .data$go_id, c("term", "go_id"), "(?=\\[GO:)"))
  }

  if (plot == FALSE) {
    return(result_table)
  }

  filtered_result_table <- result_table %>%
    dplyr::filter(n_detected_proteins_in_process >= min_n_detected_proteins_in_process)

  if (!missing(group) & y_axis_free) {
    # arrange table by group and go term for plot
    # this ensures that the terms are in the right order for a facet plot with a free axis
    filtered_result_table <- filtered_result_table %>%
      dplyr::ungroup() %>%
      dplyr::mutate(term = paste(.data$term, {{ group }}, sep = "__"))
  }

  # add cutoff for plot
  if (stringr::str_detect(plot_cutoff, pattern = "top10")) {
    split_cutoff <- stringr::str_split(plot_cutoff, pattern = " ", simplify = TRUE)
    type <- split_cutoff[1]
    plot_input <- filtered_result_table %>%
      dplyr::ungroup() %>%
      dplyr::mutate(neg_log_sig = -log10(!!rlang::ensym(type))) %>%
      dplyr::slice(1:10)
  } else {
    split_cutoff <- stringr::str_split(plot_cutoff, pattern = " ", simplify = TRUE)
    type <- split_cutoff[1]
    threshold <- as.numeric(split_cutoff[2])
    plot_input <- filtered_result_table %>%
      dplyr::ungroup() %>%
      dplyr::mutate(neg_log_sig = -log10(!!rlang::ensym(type))) %>%
      dplyr::filter(!!rlang::ensym(type) <= threshold)
  }

  # plot
  enrichment_plot <-
    {
      if ("term" %in% colnames(plot_input)) {
        ggplot2::ggplot(
          plot_input,
          ggplot2::aes(stats::reorder(.data$term, .data$neg_log_sig),
            .data$neg_log_sig,
            fill = .data$enrichment_type
          )
        )
      } else {
        ggplot2::ggplot(
          plot_input,
          ggplot2::aes(stats::reorder(.data$go_id, .data$neg_log_sig),
            .data$neg_log_sig,
            fill = .data$enrichment_type
          )
        )
      }
    } +
    ggplot2::geom_col(col = "black", size = 1) +
    {
      if (label == TRUE & nrow(plot_input) > 0) {
        geom_text(
          data = plot_input,
          aes(
            label = paste0(
              .data$n_significant_proteins_in_process,
              "/",
              .data$n_detected_proteins_in_process,
              "(",
              round(.data$n_significant_proteins_in_process / .data$n_detected_proteins_in_process * 100, digits = 1),
              "%)"
            ),
            y = .data$neg_log_sig - 0.1,
            hjust = 1
          )
        )
      }
    } +
    ggplot2::scale_fill_manual(values = c(Deenriched = "#56B4E9", Enriched = "#E76145")) +
    ggplot2::scale_y_continuous(breaks = seq(0, 100, 1)) +
    ggplot2::coord_flip() +
    {
      if (!missing(group)) {
        if (y_axis_free) {
          ggplot2::facet_wrap(rlang::new_formula(NULL, rlang::enquo(group)), scales = "free_y", ncol = facet_n_col)
        } else {
          ggplot2::facet_wrap(rlang::new_formula(NULL, rlang::enquo(group)), ncol = facet_n_col)
        }
      }
    } +
    {
      if (!missing(group)) {
        # if axis were free then the special naming that ensures the right order needs to be removed again
        scale_x_discrete(labels = function(x) gsub("__.+$", "", x))
      }
    } +
    {
      if (type == "adj_pval") {
        ggplot2::labs(
          title = plot_title,
          y = "-log10 adjusted p-value"
        )
      }
    } +
    {
      if (type == "pval") {
        ggplot2::labs(
          title = plot_title,
          y = "-log10 p-value",
          fill = "Enrichment Type"
        )
      }
    } +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 20),
      axis.text.x = ggplot2::element_text(size = 15),
      axis.text.y = ggplot2::element_text(size = 15),
      axis.title.x = ggplot2::element_text(size = 15),
      axis.title.y = ggplot2::element_blank(),
      legend.title = ggplot2::element_text(size = 15),
      legend.text = ggplot2::element_text(size = 13),
      strip.text = ggplot2::element_text(size = 15),
      strip.background = element_blank()
    )

  return(enrichment_plot)
}
