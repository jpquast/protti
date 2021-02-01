#' Check charge state distribution
#'
#' Calculates the charge state distribution for each sample (by count or intensity).
#'
#' @param data A dataframe containing at least sample names, peptide or precursor identifiers and missed cleavage counts for each peptide or precursor.
#' @param sample The column in the data dataframe containing the sample name.
#' @param grouping The column in the data dataframe containing either precursor or peptide identifiers.
#' @param charge_states The column in the data dataframe containing the different charge states assigned to the precursor or peptide.
#' @param intensity Column containing the corresponding intensity values for each peptide or precursor. Required when "intensity" is chosen as the method.
#' @param plot A logical indicating whether the result should be plotted.
#' @param method Method used for evaluation. "count" calculates the charge state distribution based on counts of the corresponding peptides or precursors in the charge state group, "intensity" calculates the percentage of precursors or peptides in each charge state group based on the corresponding intensity values.
#'
#' @return A data frame that contains the calculated percentage made up by the sum of either all counts or intensities of peptides or precursors of the corresponding charge state (depending on which method is chosen).
#' @import dplyr
#' @import ggplot2
#' @importFrom magrittr %>%
#' @importFrom forcats fct_inorder
#' @importFrom rlang .data :=
#' @importFrom tidyr drop_na
#' @importFrom utils data
#' @importFrom stringr str_sort
#' @export
#'
#' @examples
#' \dontrun{
#' qc_charge_states(
#' data,
#' sample = r_file_name,
#' grouping = pep_stripped_sequence,
#' charge_states = fg_charge,
#' intensity = NULL,
#' method = "count",
#' plot = TRUE)
#' }
qc_charge_states <-
  function(data, sample, grouping, charge_states, intensity = NULL, method, plot = FALSE)
  {
    protti_colors <- "placeholder" # assign a placeholder to prevent a missing global variable warning
    utils::data("protti_colors", envir=environment()) # then overwrite it with real data
    if (method == "count")
    {
      result <- data %>%
          dplyr::distinct({{sample}}, {{grouping}}, {{charge_states}}) %>%
          dplyr::count({{sample}}, {{charge_states}}) %>%
          dplyr::group_by({{sample}}) %>%
          dplyr::mutate(total_peptides = sum(n)) %>%
          dplyr::group_by({{sample}}, {{charge_states}}) %>%
          dplyr::summarise(charge_per = n / .data$total_peptides * 100) %>%
          dplyr::ungroup() %>%
          dplyr::mutate({{charge_states}} := forcats::fct_inorder(factor({{charge_states}}))) %>% 
          dplyr::mutate({{sample}} := factor({{sample}}, levels = unique(stringr::str_sort({{sample}}, numeric = TRUE))))

      if(plot == FALSE) {
        return(result)
      } else {
        plot <- result %>%
          ggplot2::ggplot(aes(x = {{sample}}, y = .data$charge_per, fill = {{charge_states}})) +
          geom_col(col = "black", size = 1) +
          geom_text(
            data = result %>% dplyr::filter(.data$charge_per > 5),
            aes(label = round(.data$charge_per, digits = 1)),
            position = position_stack(vjust = 0.9)
          ) +
          labs(title = "Charge distribution per .raw file",
               subtitle = "By percent of total peptide count",
               x = "",
               y = "% of total peptide count",
               fill = "Charge") +
          theme_bw() +
          theme(plot.title = ggplot2::element_text(size = 20),
                axis.title.x = ggplot2::element_text(size = 15),
                axis.text.y = ggplot2::element_text(size = 15),
                axis.text.x = ggplot2::element_text(size = 12, angle = 75, hjust =1),
                axis.title.y = ggplot2::element_text(size = 15),
                legend.title = ggplot2::element_text(size = 15),
                legend.text = ggplot2::element_text(size = 15)) +
          scale_fill_manual(values = protti_colors)
        return(plot)
      }
    }

    if (method == "intensity")
    {
      result <- data %>%
        tidyr::drop_na({{intensity}}) %>% 
        dplyr::distinct({{sample}}, {{grouping}}, {{charge_states}}, {{intensity}}) %>%
        dplyr::group_by({{sample}}) %>%
        dplyr::mutate(total_intensity = sum({{intensity}})) %>%
        dplyr::group_by({{sample}}, {{charge_states}}) %>%
        dplyr::mutate(sum_intensity_cs = sum({{intensity}})) %>%
        dplyr::summarise(charge_per = .data$sum_intensity_cs / .data$total_intensity * 100) %>%
        dplyr::ungroup() %>%
        dplyr::mutate({{charge_states}} := forcats::fct_inorder(factor({{charge_states}}))) %>%
        dplyr::mutate({{sample}} := factor({{sample}}, levels = unique(stringr::str_sort({{sample}}, numeric = TRUE)))) %>% 
        dplyr::distinct()

      if(plot == FALSE)
      {return(result)
      } else {
        plot <- result %>%
          ggplot2::ggplot(aes(x = {{sample}}, y = .data$charge_per, fill = {{charge_states}})) +
          geom_col(col = "black", size = 1) +
          geom_text(
            data = result %>% dplyr::filter(.data$charge_per > 5),
            aes(label = round(.data$charge_per, digits = 1)),
            position = position_stack(vjust = 0.9)
          ) +
          labs(title = "Charge distribution per .raw file",
               subtitle = "By percent of total intensity",
               x = "Sample",
               y = "% of total intensity",
               fill = "Charge") +
          theme_bw() +
          theme(plot.title = ggplot2::element_text(size = 20),
                axis.title.x = ggplot2::element_text(size = 15),
                axis.text.y = ggplot2::element_text(size = 15),
                axis.text.x = ggplot2::element_text(size = 12, angle = 75, hjust =1),
                axis.title.y = ggplot2::element_text(size = 15),
                legend.title = ggplot2::element_text(size = 15),
                legend.text = ggplot2::element_text(size = 15)) +
          scale_fill_manual(values = protti_colors)
        return(plot)
      }
    }
  }
