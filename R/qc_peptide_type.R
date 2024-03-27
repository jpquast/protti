#' Check peptide type percentage share
#'
#' Calculates the percentage share of each peptide types (fully-tryptic, semi-tryptic,
#' non-tryptic) for each sample.
#'
#' @param data a data frame that contains at least the input columns.
#' @param sample a character or factor column in the \code{data} data frame that contains the sample names.
#' @param peptide a character column in the \code{data} data frame that contains the peptide
#' sequence.
#' @param pep_type a character column in the \code{data} data frame that contains the peptide
#' type. Can be obtained using the \code{find_peptide} and \code{assign_peptide_type} function
#' together.
#' @param intensity a numeric column in the \code{data} data frame that contains the corresponding
#' raw or normalised intensity values (not log2) for each peptide or precursor. Required when
#' "intensity" is chosen as the method.
#' @param remove_na_intensities a logical value that specifies if sample/peptide combinations with
#' intensities that are NA (not quantified IDs) should be dropped from the data frame for analysis
#' of peptide type distributions. Default is TRUE since we are usually interested in the peptide
#' type distribution of quantifiable IDs. This is only relevant for method = "count".
#' @param method a character value that indicates the method used for evaluation.
#' \code{method = "intensity"} calculates the peptide type percentage by intensity, whereas
#' \code{method = "count"} calculates the percentage by peptide ID count. Default is
#' \code{method = count}.
#' @param plot a logical value that indicates whether the result should be plotted.
#' @param interactive a logical value that indicates whether the plot should be interactive.
#'
#' @return A data frame that contains the calculated percentage shares of each peptide type per
#' sample. The \code{count} column contains the number of peptides with a specific type. The
#' \code{peptide_type_percent} column contains the percentage share of a specific peptide type.
#' @import dplyr
#' @import tidyr
#' @import ggplot2
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @importFrom plotly ggplotly
#' @importFrom tidyr drop_na
#' @importFrom utils data
#' @importFrom stringr str_sort
#' @export
#'
#' @examples
#' # Load libraries
#' library(dplyr)
#'
#' set.seed(123) # Makes example reproducible
#'
#' # Create example data
#' data <- create_synthetic_data(
#'   n_proteins = 100,
#'   frac_change = 0.05,
#'   n_replicates = 3,
#'   n_conditions = 2,
#'   method = "effect_random"
#' ) %>%
#'   mutate(intensity_non_log2 = 2^peptide_intensity_missing)
#'
#' # Determine peptide type percentages
#' qc_peptide_type(
#'   data = data,
#'   sample = sample,
#'   peptide = peptide,
#'   pep_type = pep_type,
#'   intensity = intensity_non_log2,
#'   method = "intensity",
#'   plot = FALSE
#' )
#'
#' # Plot peptide type
#' qc_peptide_type(
#'   data = data,
#'   sample = sample,
#'   peptide = peptide,
#'   pep_type = pep_type,
#'   intensity = intensity_non_log2,
#'   method = "intensity",
#'   plot = TRUE
#' )
qc_peptide_type <- function(data,
                            sample,
                            peptide,
                            pep_type,
                            intensity,
                            remove_na_intensities = TRUE,
                            method = "count",
                            plot = FALSE,
                            interactive = FALSE) {
  protti_colours <- "placeholder" # assign a placeholder to prevent a missing global variable warning
  utils::data("protti_colours", envir = environment()) # then overwrite it with real data
  if (remove_na_intensities == TRUE) {
    data <- data %>%
      tidyr::drop_na({{ intensity }})
  }
  if (method == "count") {
    result <- data %>%
      dplyr::distinct({{ sample }}, {{ peptide }}, {{ pep_type }}) %>%
      tidyr::drop_na({{ pep_type }}) %>%
      dplyr::count({{ sample }}, {{ pep_type }}, name = "count") %>%
      dplyr::group_by({{ sample }}) %>%
      dplyr::mutate(peptide_type_percent = .data$count / sum(.data$count) * 100) %>%
      dplyr::ungroup() %>%
      dplyr::distinct({{ sample }}, {{ pep_type }}, .data$peptide_type_percent) %>%
      dplyr::mutate(pep_type = factor({{ pep_type }},
        levels = c("fully-tryptic", "semi-tryptic", "non-tryptic")
      ))

    if (is(dplyr::pull(result, {{ sample }}), "character")) {
      result <- result %>%
        dplyr::mutate({{ sample }} := factor({{ sample }},
          levels = unique(stringr::str_sort({{ sample }}, numeric = TRUE))
        ))
    }

    label_positions <- result %>%
      dplyr::group_by({{ sample }}) %>%
      dplyr::arrange(desc(.data$pep_type)) %>%
      dplyr::mutate(label_y = cumsum(.data$peptide_type_percent)) %>%
      dplyr::filter(.data$peptide_type_percent > 5)

    if (plot == TRUE & interactive == FALSE) {
      plot <- result %>%
        ggplot2::ggplot(ggplot2::aes(
          x = {{ sample }},
          y = .data$peptide_type_percent,
          fill = .data$pep_type
        )) +
        ggplot2::geom_col(col = "black", size = 1) +
        ggplot2::geom_text(
          data = label_positions,
          aes(
            y = label_y,
            label = round(.data$peptide_type_percent, digits = 1)
          ),
          vjust = 1.5
        ) +
        ggplot2::labs(
          title = "Peptide types per .raw file",
          x = "",
          y = "Percentage of peptides",
          fill = "Type"
        ) +
        ggplot2::theme_bw() +
        ggplot2::theme(
          plot.title = ggplot2::element_text(size = 20),
          axis.title.x = ggplot2::element_text(size = 15),
          axis.text.y = ggplot2::element_text(size = 15),
          axis.text.x = ggplot2::element_text(size = 12, angle = 75, hjust = 1),
          axis.title.y = ggplot2::element_text(size = 15),
          legend.title = ggplot2::element_text(size = 15),
          legend.text = ggplot2::element_text(size = 15)
        ) +
        ggplot2::scale_fill_manual(values = protti_colours)

      return(plot)
    }
    if (plot == TRUE & interactive == TRUE) {
      plot <- result %>%
        ggplot2::ggplot(ggplot2::aes({{ sample }}, .data$peptide_type_percent, fill = .data$pep_type)) +
        ggplot2::geom_col(col = "black", size = 1) +
        ggplot2::labs(
          title = "Peptide types per .raw file",
          x = "Sample",
          y = "Percentage of peptides",
          fill = "Type"
        ) +
        ggplot2::theme_bw() +
        ggplot2::theme(
          plot.title = ggplot2::element_text(size = 20),
          axis.title.x = ggplot2::element_text(size = 15),
          axis.text.y = ggplot2::element_text(size = 15),
          axis.text.x = ggplot2::element_text(size = 12, angle = 75, hjust = 1),
          axis.title.y = ggplot2::element_text(size = 15),
          legend.title = ggplot2::element_text(size = 15),
          legend.text = ggplot2::element_text(size = 15)
        ) +
        ggplot2::scale_fill_manual(values = protti_colours)

      interactive_plot <- plotly::ggplotly(plot)

      return(interactive_plot)
    }
    if (plot == FALSE) {
      return(result)
    }
  }

  if (method == "intensity") {
    result <- data %>%
      tidyr::drop_na({{ intensity }}) %>%
      dplyr::distinct({{ sample }}, {{ peptide }}, {{ pep_type }}, {{ intensity }}) %>%
      tidyr::drop_na({{ pep_type }}) %>%
      dplyr::group_by({{ sample }}) %>%
      dplyr::mutate(total_int = sum({{ intensity }})) %>%
      dplyr::group_by({{ sample }}, {{ pep_type }}) %>%
      dplyr::mutate(pep_type_int = sum({{ intensity }})) %>%
      dplyr::group_by({{ sample }}) %>%
      dplyr::mutate(peptide_type_percent = (.data$pep_type_int / .data$total_int) * 100) %>%
      dplyr::ungroup() %>%
      dplyr::distinct({{ sample }}, {{ pep_type }}, .data$peptide_type_percent) %>%
      dplyr::mutate(pep_type = factor({{ pep_type }},
        levels = c("fully-tryptic", "semi-tryptic", "non-tryptic")
      ))

    if (is(dplyr::pull(result, {{ sample }}), "character")) {
      result <- result %>%
        dplyr::mutate({{ sample }} := factor({{ sample }},
          levels = unique(stringr::str_sort({{ sample }}, numeric = TRUE))
        ))
    }

    label_positions <- result %>%
      dplyr::group_by({{ sample }}) %>%
      dplyr::arrange(desc(.data$pep_type)) %>%
      dplyr::mutate(label_y = cumsum(.data$peptide_type_percent)) %>%
      dplyr::filter(.data$peptide_type_percent > 5)

    if (plot == TRUE & interactive == FALSE) {
      plot <- result %>%
        ggplot2::ggplot(ggplot2::aes(
          x = {{ sample }},
          y = .data$peptide_type_percent,
          fill = .data$pep_type
        )) +
        ggplot2::geom_col(col = "black", size = 1) +
        ggplot2::geom_text(
          data = label_positions,
          aes(
            y = label_y,
            label = round(.data$peptide_type_percent, digits = 1)
          ),
          vjust = 1.5
        ) +
        ggplot2::labs(
          title = "Peptide type intensity per .raw file",
          y = "Percentage of total peptide intensity",
          fill = "Type"
        ) +
        ggplot2::theme_bw() +
        ggplot2::theme(
          plot.title = ggplot2::element_text(size = 20),
          axis.title.x = ggplot2::element_blank(),
          axis.text.y = ggplot2::element_text(size = 15),
          axis.text.x = ggplot2::element_text(size = 12, angle = 75, hjust = 1),
          axis.title.y = ggplot2::element_text(size = 15),
          legend.title = ggplot2::element_text(size = 15),
          legend.text = ggplot2::element_text(size = 15)
        ) +
        ggplot2::scale_fill_manual(values = protti_colours)
      return(plot)
    }
    if (plot == TRUE & interactive == TRUE) {
      plot <- result %>%
        ggplot2::ggplot(ggplot2::aes(
          x = {{ sample }},
          .data$peptide_type_percent,
          fill = .data$pep_type
        )) +
        ggplot2::geom_col(col = "black", size = 1) +
        ggplot2::labs(
          title = "Peptide type intensity per .raw file",
          y = "Percentage of total peptide intensity",
          fill = "Type"
        ) +
        ggplot2::theme_bw() +
        ggplot2::theme(
          plot.title = ggplot2::element_text(size = 20),
          axis.title.x = ggplot2::element_blank(),
          axis.text.y = ggplot2::element_text(size = 15),
          axis.text.x = ggplot2::element_text(size = 12, angle = 75, hjust = 1),
          axis.title.y = ggplot2::element_text(size = 15),
          legend.title = ggplot2::element_text(size = 15),
          legend.text = ggplot2::element_text(size = 15)
        ) +
        ggplot2::scale_fill_manual(values = protti_colours)

      interactive_plot <- plotly::ggplotly(plot)

      return(interactive_plot)
    }
    if (plot == FALSE) {
      return(result)
    }
  }
}
