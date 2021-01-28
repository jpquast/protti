#' Check missed cleavages
#'
#' Calculates the percentage of missed cleavages for each sample (by count or intensity).
#'
#' @param data A dataframe containing at least sample names, peptide or precursor identifiers and missed cleavage counts for each peptide or precursor.
#' @param sample The column in the data dataframe containing the sample name.
#' @param grouping The column in the data dataframe containing either precursor or peptide identifiers.
#' @param missed_cleavages The column in the data dataframe containing the counts of missed cleavages per peptide or precursor.
#' @param intensity Column containing the corresponding intensity values to each peptide or precursor (not log2 transformed).
#' @param remove_na_intensities Logical specifying if sample/grouping combinations with intensities that are NA (not quantified IDs) should 
#' be dropped from the data frame for analysis of missed cleavages. Default is TRUE since we are usually 
#' interested in quantifiable peptides.
#' @param plot A logical indicating whether the result should be plotted.
#' @param method Method used for evaluation. "count" calculates the percentage of missed cleavages based on counts of the corresponding peptide or precursor, "intensity" calculates the percentage of missed cleavages by intensity of the corresponding peptide or precursor.
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
#' @export
#'
#' @examples
#' \dontrun{
#' qc_missed_cleavages(
#' data,
#' sample = r_file_name,
#' grouping = pep_stripped_sequence,
#' missed_cleavages = pep_nr_of_missed_clavages,
#' intensity = fg_quantity,
#' method = "intensity",
#' plot = TRUE)
#' }
qc_missed_cleavages <-
  function(data, sample, grouping, missed_cleavages, intensity, remove_na_intensities = TRUE, method, plot = FALSE) {
    protti_colors <- "placeholder" # assign a placeholder to prevent a missing global variable warning
    utils::data("protti_colors", envir=environment()) # then overwrite it with real data
    if(remove_na_intensities == TRUE){
      data <- data %>% 
        tidyr::drop_na({{intensity}})
    }
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
      dplyr::mutate({{missed_cleavages}} := forcats::fct_inorder(factor({{missed_cleavages}}))) %>% 
      dplyr::mutate({{sample}} := factor({{sample}}, levels = unique(stringr::str_sort({{sample}}, numeric = TRUE))))

    if(plot == FALSE){
      return(result)
    } else {
        plot <- result %>%
        ggplot2::ggplot(aes(x = {{sample}}, y = .data$mc_percent, fill = {{missed_cleavages}})) +
        geom_col(col = "black", size = 1) +
        geom_text(
          data = result %>% dplyr::filter(.data$mc_percent > 5),
          aes(label = round(.data$mc_percent, digits = 1)),
          position = position_stack(vjust = 0.9)
        ) +
        labs(title = "Missed cleavages per .raw file",
             subtitle = "By percent of total peptide count",
             x = "Sample",
             y = "% of total peptide count",
             fill = "Missed cleavages") +
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
        tidyr::drop_na() %>% 
        dplyr::distinct({{sample}}, {{grouping}}, {{missed_cleavages}}, {{intensity}}) %>%
        dplyr::group_by({{sample}}) %>%
        dplyr::mutate(total_intensity = sum({{intensity}})) %>%
        dplyr::group_by({{sample}}, {{missed_cleavages}}) %>%
        dplyr::mutate(sum_intensity_mc = sum({{intensity}})) %>%
        dplyr::summarise(mc_percent = .data$sum_intensity_mc / .data$total_intensity * 100) %>%
        dplyr::ungroup() %>%
        dplyr::mutate({{missed_cleavages}} := forcats::fct_inorder(factor({{missed_cleavages}}))) %>%
        dplyr::mutate({{sample}} := factor({{sample}}, levels = unique(stringr::str_sort({{sample}}, numeric = TRUE)))) %>% 
        dplyr::distinct()


        if(plot == FALSE){
          return(result)
        } else {
          plot <- result %>%
            ggplot2::ggplot(aes(x = {{sample}}, y = .data$mc_percent, fill = {{missed_cleavages}})) +
            geom_col(col = "black", size = 1) +
            geom_text(
              data = result %>% dplyr::filter(.data$mc_percent > 5),
              aes(label = round(.data$mc_percent, digits = 1)),
              position = position_stack(vjust = 0.9)
            ) +
            labs(title = "Missed cleavages per .raw file",
                 subtitle = "By percent of total intensity",
                 x = "Sample",
                 y = "% of total intensity",
                 fill = "Missed cleavages") +
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


