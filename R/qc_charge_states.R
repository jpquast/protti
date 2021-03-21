#' Check charge state distribution
#'
#' Calculates the charge state distribution for each sample (by count or intensity).
#'
#' @param data A data frame containing at least sample names, peptide or precursor identifiers and missed cleavage counts for each peptide or precursor.
#' @param sample the column in the data data frame containing the sample name.
#' @param grouping the column in the data data frame containing either precursor or peptide identifiers.
#' @param charge_states the column in the data data frame containing the different charge states assigned to the precursor or peptide.
#' @param intensity the name of the column containing the corresponding raw or normalised intensity values (not log2) for each peptide or precursor. Required when "intensity" is chosen as the method.
#' @param remove_na_intensities logical specifying if sample/grouping combinations with intensities that are NA (not quantified IDs) should
#' be dropped from the data frame for analysis of missed cleavages. Default is TRUE since we are usually
#' interested in quantifiable peptides. This is only relevant for method = "count".
#' @param plot logical indicating whether the result should be plotted.
#' @param method character vector indicating the method used for evaluation. "count" calculates the charge state distribution based on counts of the corresponding peptides or precursors in the charge state group, "intensity" calculates the percentage of precursors or peptides in each charge state group based on the corresponding intensity values.
#' @param interactive argument specifying whether the plot should be interactive (default is FALSE).
#'
#' @return A data frame that contains the calculated percentage made up by the sum of either all counts or intensities of peptides or precursors of the corresponding charge state (depending on which method is chosen).
#' @import dplyr
#' @import ggplot2
#' @importFrom magrittr %>%
#' @importFrom forcats fct_inorder
#' @importFrom plotly ggplotly
#' @importFrom rlang .data :=
#' @importFrom tidyr drop_na
#' @importFrom utils data
#' @importFrom stringr str_sort
#' @export
#'
#' @examples
#' \dontrun{
#' qc_charge_states(
#'   data,
#'   sample = r_file_name,
#'   grouping = pep_stripped_sequence,
#'   charge_states = fg_charge,
#'   method = "count",
#'   plot = TRUE
#' )
#' }
qc_charge_states <-
  function(data, sample, grouping, charge_states, intensity = NULL, remove_na_intensities = TRUE, method = "count", plot = FALSE, interactive = FALSE) {
    protti_colours <- "placeholder" # assign a placeholder to prevent a missing global variable warning
    utils::data("protti_colours", envir = environment()) # then overwrite it with real data
    if (remove_na_intensities == TRUE) {
      if (missing(intensity)) stop("Please provide a column containing intensities or set remove_na_intensities to FALSE")

      data <- data %>%
        tidyr::drop_na({{ intensity }})
    }
    if (method == "count") {
      result <- data %>%
        dplyr::distinct({{ sample }}, {{ grouping }}, {{ charge_states }}) %>%
        dplyr::count({{ sample }}, {{ charge_states }}) %>%
        dplyr::group_by({{ sample }}) %>%
        dplyr::mutate(total_peptides = sum(n)) %>%
        dplyr::group_by({{ sample }}, {{ charge_states }}) %>%
        dplyr::summarise(charge_per = n / .data$total_peptides * 100) %>%
        dplyr::ungroup() %>%
        dplyr::mutate({{ charge_states }} := forcats::fct_inorder(factor({{ charge_states }}))) %>%
        dplyr::mutate({{ sample }} := factor({{ sample }}, levels = unique(stringr::str_sort({{ sample }}, numeric = TRUE))))

      if (plot == FALSE) {
        return(result)
      } else {
        plot <- result %>%
          ggplot2::ggplot(aes(x = {{ sample }}, y = .data$charge_per, fill = {{ charge_states }})) +
          geom_col(col = "black", size = 1) +
          {
            if (interactive == FALSE) {
              geom_text(
                data = result %>% dplyr::filter(.data$charge_per > 5),
                aes(label = round(.data$charge_per, digits = 1)),
                position = position_stack(vjust = 0.9)
              )
            }
          } +
          labs(
            title = "Charge distribution per .raw file",
            subtitle = "By percent of total peptide count",
            x = "",
            y = "% of total peptide count",
            fill = "Charge"
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
        dplyr::distinct({{ sample }}, {{ grouping }}, {{ charge_states }}, {{ intensity }}) %>%
        dplyr::group_by({{ sample }}) %>%
        dplyr::mutate(total_intensity = sum({{ intensity }})) %>%
        dplyr::group_by({{ sample }}, {{ charge_states }}) %>%
        dplyr::mutate(sum_intensity_cs = sum({{ intensity }})) %>%
        dplyr::summarise(charge_per = .data$sum_intensity_cs / .data$total_intensity * 100) %>%
        dplyr::ungroup() %>%
        dplyr::mutate({{ charge_states }} := forcats::fct_inorder(factor({{ charge_states }}))) %>%
        dplyr::mutate({{ sample }} := factor({{ sample }}, levels = unique(stringr::str_sort({{ sample }}, numeric = TRUE)))) %>%
        dplyr::distinct()

      if (plot == FALSE) {
        return(result)
      } else {
        plot <- result %>%
          ggplot2::ggplot(aes(x = {{ sample }}, y = .data$charge_per, fill = {{ charge_states }})) +
          geom_col(col = "black", size = 1) +
          {
            if (interactive == FALSE) {
              geom_text(
                data = result %>% dplyr::filter(.data$charge_per > 5),
                aes(label = round(.data$charge_per, digits = 1)),
                position = position_stack(vjust = 0.9)
              )
            }
          } +
          labs(
            title = "Charge distribution per .raw file",
            subtitle = "By percent of total intensity",
            x = "Sample",
            y = "% of total intensity",
            fill = "Charge"
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
