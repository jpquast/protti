#' Ranked intensities
#'
#' Plots ranked intensities
#'
#' @param data a data frame that contains at least the input variables.
#' @param sample a character column in the \code{data} data frame that contains the sample names.
#' @param grouping a character column in the \code{data} data frame that contains protein, precursor,
#' or peptide identifiers.
#' @param log2_intensity a numeric column in the \code{data} data frame that contains the log2
#' transfromed intensities of the selected grouping variable.
#' @param facet a logical value that specifies whether the plot should be faceted by samples
#' (default is FALSE).If \code{facet = FALSE} the median of each protein intensity will be shown in the plot.
#' @param interactive a logical value that specifies whether the plot should be interactive
#' (default is FALSE).
#'
#' @return A plot showing the ranked intensities is returned.
#' @import dplyr
#' @import ggplot2
#' @importFrom magrittr %>%
#' @importFrom tidyr drop_na
#' @importFrom plotly ggplotly
#' @importFrom utils data
#' @export
#'
#' @examples
#' set.seed(123) # Makes example reproducible
#'
#' # Create synthetic data
#' data <- create_synthetic_data(
#'   n_proteins = 10,
#'   frac_change = 0.5,
#'   n_replicates = 4,
#'   n_conditions = 3,
#'   method = "effect_random",
#'   additional_metadata = FALSE
#' )
#'
#' # Calculate protein abundance
#' protein_intensities <- calculate_protein_abundance(
#'                            data = data,
#'                            sample = sample, 
#'                            protein_id = protein, 
#'                            precursor = peptide,
#'                            peptide = peptide, 
#'                            intensity_log2 = peptide_intensity)
#'                           
#' # Plot ranked intensities for all samples combined
#' qc_ranked_intensities(data = protein_intensities, 
#'                      sample = sample, 
#'                      grouping = protein, 
#'                      log2_intensity = peptide_intensity, 
#'                      interactive = TRUE)
#'
#' # Plot ranked intensities for each sample separately
#' qc_ranked_intensities(data = protein_intensities, 
#'                      sample = sample, 
#'                      grouping = protein, 
#'                      log2_intensity = peptide_intensity, 
#'                      interactive = TRUE, 
#'                      facet = TRUE)
#'
qc_ranked_intensities <- function(data, 
                                sample,
                                grouping, 
                                log2_intensity,
                                facet = FALSE, 
                                interactive = FALSE) {
  if(facet == FALSE){
    input <- data %>%
      dplyr::distinct({{ sample }}, {{ grouping }}, {{ log2_intensity }}) %>%
      tidyr::drop_na({{ log2_intensity }}) %>%
      dplyr::mutate(intensity = 2^{{ log2_intensity}}) %>%
      dplyr::group_by({{ grouping }}) %>%
      dplyr::summarize(median_intensity = stats::median(log10(.data$intensity), na.rm = TRUE), .groups = "drop") %>%
      dplyr::distinct({{ grouping }}, .data$median_intensity) %>%
      dplyr::arrange(desc(.data$median_intensity)) %>%
      dplyr::mutate(rank = dplyr::row_number())
    
    plot <- input %>%
      ggplot2::ggplot(ggplot2::aes(.data$rank, .data$median_intensity, text = paste("grouping:" = {{ grouping }}))) +
      ggplot2::geom_point(size = 3,
                          colour = "#5680C1") +
      ggplot2::labs(title = "Ranked intensitites", x = "Rank", y = "Median log10 intensity") +
      ggplot2::theme_bw() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(size = 20),
        axis.title.x = ggplot2::element_text(size = 15),
        axis.text.y = ggplot2::element_text(size = 15),
        axis.text.x = ggplot2::element_text(size = 12),
        axis.title.y = ggplot2::element_text(size = 15),
        strip.text = ggplot2::element_text(size = 15),
        panel.border = ggplot2::element_rect(fill = NA),
        strip.background = element_rect(fill = "white")
      )
    
    if (interactive == FALSE) {
      return(plot)
    }
    suppressWarnings(plotly::ggplotly(plot))
  }else{
    input <- data %>%
      dplyr::distinct({{ sample }}, {{ grouping }}, {{ log2_intensity }}) %>%
      tidyr::drop_na({{ log2_intensity }}) %>%
      dplyr::mutate(intensity = 2^{{ log2_intensity}}) %>%
      dplyr::group_by({{ grouping }}, {{ sample }}) %>%
      dplyr::mutate(log10_intensity = log10(.data$intensity)) %>%
      dplyr::distinct({{ sample }}, {{ grouping }}, .data$log10_intensity) %>%
      dplyr::group_by({{ sample }}) %>%
      dplyr::arrange(desc(.data$log10_intensity)) %>%
      dplyr::mutate(rank = dplyr::row_number())
    
    plot <- input %>%
      ggplot2::ggplot(ggplot2::aes(.data$rank, .data$log10_intensity, text = paste("grouping:" = {{ grouping }}))) +
      ggplot2::geom_point(colour = "#5680C1") +
      ggplot2::labs(title = "Ranked intensitites", x = "Rank", y = "Median log10 intensity") +
      ggplot2::facet_wrap(rlang::new_formula(NULL, rlang::enquo(sample)), scales = "free", ncol = 4) +
      ggplot2::theme_bw() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(size = 20),
        axis.title.x = ggplot2::element_text(size = 15),
        axis.text.y = ggplot2::element_text(size = 15),
        axis.text.x = ggplot2::element_text(size = 12),
        axis.title.y = ggplot2::element_text(size = 15),
        strip.text = ggplot2::element_text(size = 15),
        panel.border = ggplot2::element_rect(fill = NA),
        strip.background = element_rect(fill = "white")
      )
    
    if (interactive == FALSE) {
      return(plot)
    }
    suppressWarnings(plotly::ggplotly(plot))
  }
}
