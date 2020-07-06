#' Check missed cleavages
#'
#' Calculates the percentage of missed cleavages for each sample (by count or intensity).
#'
#' @param data A dataframe containing at least sample names, peptide or precursor identifiers and missed cleavage counts for each peptide or precursor.
#' @param sample The column in the data dataframe containing the sample name.
#' @param grouping The column in the data dataframe containing either precursor or peptide identifiers.
#' @param missed_cleavages The column in the data dataframe containing the counts of missed cleavages per peptide or precursor.
#' @param intensity Optional column containing the corresponding intensity values to each peptide or precursor.
#' @param plot A logical indicating whether the result should be plotted.
#' @param method Method used for evaluation. "Count" calculates the percentage of missed cleavages based on counts or the corresponding peptide or precursor, "intensity" calculates the percentage of missed cleavages by intensity of the corresponding peptide or precursor.
#'
#' @return A data frame that contains the calculated percentage made up by the sum of all peptides or precursors containing the corresponding amount of missed cleavages.
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
#' qc_missed_cleavages(data, sample, grouping, missed_cleavages, intensity = NULL, plot = TRUE)
#' }
qc_missed_cleavages <-
  function(data, sample, grouping, missed_cleavages, intensity = NULL, method, plot = FALSE)
  {
    if (method == "count")
    {
    result <- data %>%
      dplyr::distinct({{sample}}, {{grouping}}, {{missed_cleavages}}, {{intensity}}) %>%
      dplyr::count({{sample}}, {{missed_cleavages}}) %>%
      dplyr::group_by({{sample}}) %>%
      dplyr::mutate(total_peptide_count = sum(n)) %>%
      dplyr::group_by({{sample}}, {{missed_cleavages}}) %>%
      dplyr::summarise(mc_percent = n / .data$total_peptide_count * 100) %>%
      dplyr::ungroup() %>%
      dplyr::mutate({{missed_cleavages}} := forcats::fct_inorder(factor({{missed_cleavages}})))

    if(plot == FALSE)
    {return(result)
    } else {
        result %>%
        ggplot2::ggplot(aes(x = {{sample}}, y = .data$mc_percent, fill = {{missed_cleavages}})) +
        geom_col(col = "black") +
        geom_text(
          data = subset(result, .data$mc_percent > 5),
          aes(label = round(.data$mc_percent, digits = 1)),
          position = position_stack(vjust = 0.5)
        ) +
        labs(title = "Missed cleavages per .raw file",
             subtitle = "By percent of total peptide count",
             x = "Sample",
             y = "% of total peptide count",
             fill = "Missed cleavages") +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
    }
   }

    if (method == "intensity")
    {
      result <- data %>%
        dplyr::distinct({{sample}}, {{grouping}}, {{missed_cleavages}}, {{intensity}}) %>%
        dplyr::group_by({{sample}}) %>%
        dplyr::mutate(total_intensity = sum({{intensity}})) %>%
        dplyr::group_by({{sample}}, {{missed_cleavages}}) %>%
        dplyr::mutate(sum_intensity_mc = sum({{intensity}})) %>%
        dplyr::summarise(mc_percent = .data$sum_intensity_mc / .data$total_intensity * 100) %>%
        dplyr::ungroup() %>%
        dplyr::mutate({{missed_cleavages}} := forcats::fct_inorder(factor({{missed_cleavages}}))) %>%
        dplyr::distinct()


        if(plot == FALSE)
        {return(result)
        } else {
          result %>%
            ggplot2::ggplot(aes(x = {{sample}}, y = .data$mc_percent, fill = {{missed_cleavages}})) +
            geom_col(col = "black") +
            geom_text(
              data = subset(result, .data$mc_percent > 5),
              aes(label = round(.data$mc_percent, digits = 1)),
              position = position_stack(vjust = 0.5)
            ) +
            labs(title = "Missed cleavages per .raw file",
                 subtitle = "By percent of total intensity",
                 x = "Sample",
                 y = "% of total intensity",
                 fill = "Missed cleavages") +
            theme_bw() +
            theme(axis.text.x = element_text(angle = 45, hjust = 1))
        }
      }
    }


