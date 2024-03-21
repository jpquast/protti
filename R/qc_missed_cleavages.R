#' Check missed cleavages
#'
#' Calculates the percentage of missed cleavages for each sample (by count or intensity). The
#' default settings remove grouping variables without quantitative information (intensity is NA).
#' These will not be used for the calculation of missed cleavage percentages.
#'
#' @param data a data frame containing at least sample names, peptide or precursor identifiers
#' and missed cleavage counts for each peptide or precursor.
#' @param sample a character or factor column in the \code{data} data frame that contains the sample name.
#' @param grouping a character column in the \code{data} data frame that contains either precursor or
#' peptide identifiers.
#' @param missed_cleavages a numeric column in the \code{data} data frame that contains the counts
#' of missed cleavages per peptide or precursor.
#' @param intensity a numeric column in the \code{data} data frame that contains the corresponding
#' raw or normalised intensity values (not log2) for each peptide or precursor. Required when
#' "intensity" is chosen as the method.
#' @param remove_na_intensities a logical value that specifies if sample/grouping combinations with
#' intensities that are NA (not quantified IDs) should be dropped from the data frame for analysis
#' of missed cleavages. Default is TRUE since we are usually interested in quantifiable peptides.
#' This is only relevant for method = "count".
#' @param plot a logical value that indicates whether the result should be plotted.
#' @param method a character value that indicates the method used for evaluation. "count"
#' calculates the percentage of missed cleavages based on counts of the corresponding peptide or
#' precursor, "intensity" calculates the percentage of missed cleavages by intensity of the
#' corresponding peptide or precursor.
#' @param interactive a logical value that specifies whether the plot should be interactive
#' (default is FALSE).
#'
#' @return A data frame that contains the calculated percentage made up by the sum of all peptides
#' or precursors containing the corresponding amount of missed cleavages.
#' @import dplyr
#' @import ggplot2
#' @importFrom magrittr %>%
#' @importFrom forcats fct_inorder
#' @importFrom rlang .data :=
#' @importFrom tidyr drop_na
#' @importFrom utils data
#' @importFrom stringr str_sort
#' @importFrom plotly ggplotly
#' @export
#'
#' @examples
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
#' # Calculate missed cleavage percentages
#' qc_missed_cleavages(
#'   data = data,
#'   sample = sample,
#'   grouping = peptide,
#'   missed_cleavages = n_missed_cleavage,
#'   intensity = intensity_non_log2,
#'   method = "intensity",
#'   plot = FALSE
#' )
#'
#' # Plot missed cleavages
#' qc_missed_cleavages(
#'   data = data,
#'   sample = sample,
#'   grouping = peptide,
#'   missed_cleavages = n_missed_cleavage,
#'   intensity = intensity_non_log2,
#'   method = "intensity",
#'   plot = TRUE
#' )
qc_missed_cleavages <-
  function(data,
           sample,
           grouping,
           missed_cleavages,
           intensity,
           remove_na_intensities = TRUE,
           method = "count",
           plot = FALSE,
           interactive = FALSE) {
    protti_colours <- "placeholder" # assign a placeholder to prevent a missing global variable warning
    utils::data("protti_colours", envir = environment()) # then overwrite it with real data
    if (remove_na_intensities == TRUE) {
      if (missing(intensity)) {
        stop(strwrap("Please provide a column containing
intensities or set remove_na_intensities to FALSE",
          prefix = "\n", initial = ""
        ))
      }

      data <- data %>%
        tidyr::drop_na({{ intensity }})
    }
    if (method == "count") {
      result <- data %>%
        dplyr::distinct({{ sample }}, {{ grouping }}, {{ missed_cleavages }}, {{ intensity }}) %>%
        dplyr::count({{ sample }}, {{ missed_cleavages }}) %>%
        dplyr::group_by({{ sample }}) %>%
        dplyr::mutate(total_peptide_count = sum(n)) %>%
        dplyr::group_by({{ sample }}, {{ missed_cleavages }}) %>%
        dplyr::summarise(mc_percent = n / .data$total_peptide_count * 100) %>%
        dplyr::ungroup() %>%
        dplyr::mutate({{ missed_cleavages }} := forcats::fct_inorder(factor({{ missed_cleavages }})))

      if (is(dplyr::pull(result, {{ sample }}), "character")) {
        result <- result %>%
          dplyr::mutate({{ sample }} := factor({{ sample }},
            levels = unique(stringr::str_sort({{ sample }}, numeric = TRUE))
          ))
      }

      label_positions <- result %>%
        dplyr::group_by({{ sample }}) %>%
        dplyr::arrange(desc({{ missed_cleavages }})) %>%
        dplyr::mutate(label_y = cumsum(.data$mc_percent)) %>%
        dplyr::filter(.data$mc_percent > 5)

      if (plot == FALSE) {
        return(result)
      } else {
        plot <- result %>%
          ggplot2::ggplot(aes(
            x = {{ sample }},
            y = .data$mc_percent,
            fill = {{ missed_cleavages }}
          )) +
          ggplot2::geom_col(col = "black", size = 1) +
          {
            if (interactive == FALSE) {
              ggplot2::geom_text(
                data = label_positions,
                aes(y = .data$label_y,
                  label = round(.data$mc_percent, digits = 1)),
                  vjust = 1.5
              )
            }
          } +
          ggplot2::labs(
            title = "Missed cleavages per .raw file",
            subtitle = "By percent of total peptide count",
            y = "% of total peptide count",
            fill = "Missed cleavages"
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
      }
    }

    if (method == "intensity") {
      result <- data %>%
        tidyr::drop_na({{ intensity }}) %>%
        dplyr::distinct(
          {{ sample }},
          {{ grouping }},
          {{ missed_cleavages }},
          {{ intensity }}
        ) %>%
        dplyr::group_by({{ sample }}) %>%
        dplyr::mutate(total_intensity = sum({{ intensity }})) %>%
        dplyr::group_by({{ sample }}, {{ missed_cleavages }}) %>%
        dplyr::mutate(sum_intensity_mc = sum({{ intensity }})) %>%
        dplyr::reframe(mc_percent = .data$sum_intensity_mc / .data$total_intensity * 100) %>%
        dplyr::mutate({{ missed_cleavages }} := forcats::fct_inorder(factor({{ missed_cleavages }}))) %>%
        dplyr::distinct()

      if (is(dplyr::pull(result, {{ sample }}), "character")) {
        result <- result %>%
          dplyr::mutate({{ sample }} := factor({{ sample }},
            levels = unique(stringr::str_sort({{ sample }}, numeric = TRUE))
          ))
      }

      label_positions <- result %>%
        dplyr::group_by({{ sample }}) %>%
        dplyr::arrange(desc({{ missed_cleavages }})) %>%
        dplyr::mutate(label_y = cumsum(.data$mc_percent)) %>%
        dplyr::filter(.data$mc_percent > 5)

      if (plot == FALSE) {
        return(result)
      } else {
        plot <- result %>%
          ggplot2::ggplot(aes(
            x = {{ sample }},
            y = .data$mc_percent,
            fill = {{ missed_cleavages }}
          )) +
          ggplot2::geom_col(col = "black", size = 1) +
          {
            if (interactive == FALSE) {
              ggplot2::geom_text(
                data = label_positions,
                aes(y = .data$label_y,
                    label = round(.data$mc_percent, digits = 1)),
                    vjust = 1.5
              )
            }
          } +
          ggplot2::labs(
            title = "Missed cleavages per .raw file",
            subtitle = "By percent of total intensity",
            y = "% of total intensity",
            fill = "Missed cleavages"
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
      }
    }
    if (interactive == TRUE) {
      return(plotly::ggplotly(plot))
    } else {
      return(plot)
    }
  }
