#' Check missed cleavages
#'
#' Calculates the percentage of missed cleavages for each sample (by count or intensity).The default settings remove grouping variables without quantitative information (intensity is NA).
#' These will not be used for the calculation of missed cleavage percentages.
#'
#' @param data A data frame containing at least sample names, peptide or precursor identifiers and missed cleavage counts for each peptide or precursor.
#' @param sample the name of the column in the data data frame containing the sample name.
#' @param grouping the name of the column in the data data frame containing either precursor or peptide identifiers.
#' @param missed_cleavages the name of the column in the data data frame containing the counts of missed cleavages per peptide or precursor.
#' @param intensity the name of the column containing the corresponding raw or normalised intensity values (not log2) for each peptide or precursor. Required when "intensity" is chosen as the method.
#' @param remove_na_intensities Logical specifying if sample/grouping combinations with intensities that are NA (not quantified IDs) should
#' be dropped from the data frame for analysis of missed cleavages. Default is TRUE since we are usually
#' interested in quantifiable peptides. This is only relevant for method = "count".
#' @param plot A logical indicating whether the result should be plotted.
#' @param method character vector indicating the method used for evaluation. "count" calculates the percentage of missed cleavages based on counts of the corresponding peptide or precursor, "intensity" calculates the percentage of missed cleavages by intensity of the corresponding peptide or precursor.
#' @param interactive Argument specifying whether the plot should be interactive (default is FALSE).
#'
#' @return A data frame that contains the calculated percentage made up by the sum of all peptides or precursors containing the corresponding amount of missed cleavages.
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
#' \dontrun{
#' qc_missed_cleavages(
#'   data,
#'   sample = r_file_name,
#'   grouping = pep_stripped_sequence,
#'   missed_cleavages = pep_nr_of_missed_cleavages,
#'   intensity = fg_quantity,
#'   method = "intensity",
#'   plot = TRUE
#' )
#' }
qc_missed_cleavages <-
  function(data, sample, grouping, missed_cleavages, intensity, remove_na_intensities = TRUE, method = "count", plot = FALSE, interactive = FALSE) {
    protti_colours <- "placeholder" # assign a placeholder to prevent a missing global variable warning
    utils::data("protti_colours", envir = environment()) # then overwrite it with real data
    if (remove_na_intensities == TRUE) {
      if (missing(intensity)) stop("Please provide a column containing intensities or set remove_na_intensities to FALSE")

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
        dplyr::mutate({{ missed_cleavages }} := forcats::fct_inorder(factor({{ missed_cleavages }}))) %>%
        dplyr::mutate({{ sample }} := factor({{ sample }}, levels = unique(stringr::str_sort({{ sample }}, numeric = TRUE))))

      if (plot == FALSE) {
        return(result)
      } else {
        plot <- result %>%
          ggplot2::ggplot(aes(x = {{ sample }}, y = .data$mc_percent, fill = {{ missed_cleavages }})) +
          geom_col(col = "black", size = 1) +
          {
            if (interactive == FALSE) {
              geom_text(
                data = result %>% dplyr::filter(.data$mc_percent > 5),
                aes(label = round(.data$mc_percent, digits = 1)),
                position = position_stack(vjust = 0.9)
              )
            }
          } +
          labs(
            title = "Missed cleavages per .raw file",
            subtitle = "By percent of total peptide count",
            x = "Sample",
            y = "% of total peptide count",
            fill = "Missed cleavages"
          ) +
          theme_bw() +
          theme(
            plot.title = ggplot2::element_text(size = 20),
            axis.title.x = ggplot2::element_text(size = 15),
            axis.text.y = ggplot2::element_text(size = 15),
            axis.text.x = ggplot2::element_text(size = 12, angle = 75, hjust = 1),
            axis.title.y = ggplot2::element_text(size = 15),
            legend.title = ggplot2::element_text(size = 15),
            legend.text = ggplot2::element_text(size = 15)
          ) +
          scale_fill_manual(values = protti_colours)
      }
    }

    if (method == "intensity") {
      result <- data %>%
        tidyr::drop_na({{ intensity }}) %>%
        dplyr::distinct({{ sample }}, {{ grouping }}, {{ missed_cleavages }}, {{ intensity }}) %>%
        dplyr::group_by({{ sample }}) %>%
        dplyr::mutate(total_intensity = sum({{ intensity }})) %>%
        dplyr::group_by({{ sample }}, {{ missed_cleavages }}) %>%
        dplyr::mutate(sum_intensity_mc = sum({{ intensity }})) %>%
        dplyr::summarise(mc_percent = .data$sum_intensity_mc / .data$total_intensity * 100) %>%
        dplyr::ungroup() %>%
        dplyr::mutate({{ missed_cleavages }} := forcats::fct_inorder(factor({{ missed_cleavages }}))) %>%
        dplyr::mutate({{ sample }} := factor({{ sample }}, levels = unique(stringr::str_sort({{ sample }}, numeric = TRUE)))) %>%
        dplyr::distinct()


      if (plot == FALSE) {
        return(result)
      } else {
        plot <- result %>%
          ggplot2::ggplot(aes(x = {{ sample }}, y = .data$mc_percent, fill = {{ missed_cleavages }})) +
          geom_col(col = "black", size = 1) +
          {
            if (interactive == FALSE) {
              geom_text(
                data = result %>% dplyr::filter(.data$mc_percent > 5),
                aes(label = round(.data$mc_percent, digits = 1)),
                position = position_stack(vjust = 0.9)
              )
            }
          } +
          labs(
            title = "Missed cleavages per .raw file",
            subtitle = "By percent of total intensity",
            x = "",
            y = "% of total intensity",
            fill = "Missed cleavages"
          ) +
          theme_bw() +
          theme(
            plot.title = ggplot2::element_text(size = 20),
            axis.title.x = ggplot2::element_text(size = 15),
            axis.text.y = ggplot2::element_text(size = 15),
            axis.text.x = ggplot2::element_text(size = 12, angle = 75, hjust = 1),
            axis.title.y = ggplot2::element_text(size = 15),
            legend.title = ggplot2::element_text(size = 15),
            legend.text = ggplot2::element_text(size = 15)
          ) +
          scale_fill_manual(values = protti_colours)
      }
    }
    if (interactive == TRUE) {
      return(plotly::ggplotly(plot))
    } else {
      return(plot)
    }
  }
