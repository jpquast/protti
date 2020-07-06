#' Check charge state distribution
#'
#' Calculates the charge state distribution for each sample (by count or intensity).
#'
#' @param data A dataframe containing at least sample names, peptide or precursor identifiers and missed cleavage counts for each peptide or precursor.
#' @param sample The column in the data dataframe containing the sample name.
#' @param grouping The column in the data dataframe containing either precursor or peptide identifiers.
#' @param charge_states The column in the data dataframe containing the different charge states assigned to the precursor or peptide.
#' @param intensity Optional column containing the corresponding intensity values to each peptide or precursor.
#' @param plot A logical indicating whether the result should be plotted.
#' @param method Method used for evaluation. "count" calculates the charge state distribution based on counts of the corresponding peptides or precursors in the charge state group, "intensity" calculates the percentage of precursors or peptides in each charge state group.
#'
#' @return A data frame that contains the calculated percentage made up by the sum of all peptides or precursors of the corresponding of a specific charge state.
#' @import dplyr
#' @import ggplot2
#' @importFrom magrittr %>%
#' @importFrom forcats fct_inorder
#' @importFrom rlang .data
#' @importFrom rlang :=
#' @export
#'
#' @examples
#' \dontrun{
#' qc_charge_states(data, sample, grouping, charge_startes, intensity = NULL, plot = TRUE)
#' }

qc_charge_states <-
  function(data, sample, grouping, charge_states, intensity = NULL, method, plot = FALSE)
  {
    if (method == "count")
    {
      result <- data %>%
          dplyr::count({{sample}}, {{charge_states}}) %>%
          dplyr::group_by({{sample}}) %>%
          dplyr::mutate(total_peptides = sum(n)) %>%
          dplyr::group_by({{sample}}, {{charge_states}}) %>%
          dplyr::summarise(charge_per = n / .data$total_peptides * 100) %>%
          dplyr::ungroup() %>%
          dplyr::mutate({{charge_states}} := forcats::fct_inorder(factor({{charge_states}})))

      if(plot == FALSE)
      {return(result)
      } else {
        plot <- result %>%
          ggplot2::ggplot(aes(x = {{sample}}, y = .data$charge_per, fill = {{charge_states}})) +
          geom_col(col = "black") +
          geom_text(
            data = result %>% dplyr::filter(.data$charge_per > 5),
            aes(label = round(.data$charge_per, digits = 1)),
            position = position_stack(vjust = 0.5)
          ) +
          labs(title = "Charge distribution per .raw file",
               subtitle = "By percent of total peptide count",
               x = "Sample",
               y = "% of total peptide count",
               fill = "Charge") +
          theme_bw() +
          theme(axis.text.x = element_text(angle = 45, hjust = 1))
        return(plot)
      }
    }

    if (method == "intensity")
    {
      result <- data %>%
        dplyr::distinct({{sample}}, {{grouping}}, {{charge_states}}, {{intensity}}) %>%
        dplyr::group_by({{sample}}) %>%
        dplyr::mutate(total_intensity = sum({{intensity}})) %>%
        dplyr::group_by({{sample}}, {{charge_states}}) %>%
        dplyr::mutate(sum_intensity_cs = sum({{intensity}})) %>%
        dplyr::summarise(charge_per = .data$sum_intensity_cs / .data$total_intensity * 100) %>%
        dplyr::ungroup() %>%
        dplyr::mutate({{charge_states}} := forcats::fct_inorder(factor({{charge_states}}))) %>%
        dplyr::distinct()

      if(plot == FALSE)
      {return(result)
      } else {
        plot <- result %>%
          ggplot2::ggplot(aes(x = {{sample}}, y = .data$charge_per, fill = {{charge_states}})) +
          geom_col(col = "black") +
          geom_text(
            data = result %>% dplyr::filter(.data$charge_per > 5),
            aes(label = round(.data$charge_per, digits = 1)),
            position = position_stack(vjust = 0.5)
          ) +
          labs(title = "Charge distribution per .raw file",
               subtitle = "By percent of total intensity",
               x = "Sample",
               y = "% of total intensity",
               fill = "Charge") +
          theme_bw() +
          theme(axis.text.x = element_text(angle = 45, hjust = 1))
        return(plot)
      }
    }
  }
