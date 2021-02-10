#' Peptide abundance profile plot
#'
#' Creates a plot of peptide abundances across samples. This is helpful to investigate effects of peptide and protein abundance changes
#' in different samples and conditions.
#'
#' @param data Data frame containing at least the input variables.
#' @param sample Column in the data frame containing sample names.
#' @param peptide Column in the data frame containing peptide or precursor names.
#' @param intensity Column in the data frame containing log2 transformed intensities.
#' @param grouping Column in the data frame containing groups by which the data should be split. This can be for example protein IDs.
#' @param targets Character vector specifying elements of the grouping column for which a plot should be returned.
#' @param split_all Logical, specifying if a plot for each element of the grouping column should be returned. The default is FALSE. Be careful
#' with this arguments as the number of plots can be quite high.
#' @param protein_abundance_plot Logical, if the input for this plot comes directly from \code{calculate_protein_abundance} this argument
#' can be set to \code{TRUE}. This displays all peptides in gray, while the protein abundance is displayed in green.
#'
#' @return A list of peptide profile plots.
#' @import ggplot2
#' @importFrom magrittr %>%
#' @importFrom dplyr distinct pull filter
#' @importFrom tidyr drop_na
#' @importFrom rlang !! ensym
#' @importFrom purrr map2
#' @importFrom utils data
#' @export
#'
#' @examples
#' \dontrun{
#' plot_peptide_abundance(
#'   data,
#'   sample = r_file_name,
#'   peptide = eg_precursor_id,
#'   intensity = log2_intensity,
#'   grouping = pg_protein_accessions,
#'   targets = c("P03421")
#' )
#' }
plot_peptide_profiles <- function(data, sample, peptide, intensity, grouping, targets, split_all = FALSE, protein_abundance_plot = FALSE) {
  . <- NULL
  protti_colours <- "placeholder" # assign a placeholder to prevent a missing global variable warning
  utils::data("protti_colours", envir = environment()) # then overwrite it with real data
  if (split_all == FALSE) {
    if (missing(targets)) stop("Please provide at least one target to plot!")
    input <- data %>%
      dplyr::distinct({{ sample }}, {{ peptide }}, {{ intensity }}, {{ grouping }}) %>%
      tidyr::drop_na({{ intensity }}) %>% 
      dplyr::filter({{ grouping }} %in% targets) %>%
      split(dplyr::pull(., !!ensym(grouping)))
  }
  if (split_all == TRUE) {
    groups <- length(unique(dplyr::pull(data, {{ grouping }})))
    message("Splitting into ", groups, " groups and returning ", groups, " plots.")
    input <- data %>%
      dplyr::distinct({{ sample }}, {{ peptide }}, {{ intensity }}, {{ grouping }}) %>%
      tidyr::drop_na({{ intensity }}) %>% 
      split(dplyr::pull(., !!ensym(grouping)))
  }

  plot <- purrr::map2(
    .x = input,
    .y = names(input),
    .f = ~ {
      ggplot2::ggplot(.x, ggplot2::aes({{ sample }}, {{ intensity }}, group = {{ peptide }}, col = {{ peptide }})) +
        ggplot2::geom_point() +
        ggplot2::geom_line(size = 1) +
        ggplot2::labs(title = paste("Peptide profiles:", .y), x = "Sample", y = "Intensity [log2]", col = "Peptides") +
        ggplot2::theme_bw() +
        {if(length(unique(dplyr::pull(.x, {{ peptide }}))) > 25) ggplot2::theme(legend.position = "none")} +
        ggplot2::theme(
          plot.title = ggplot2::element_text(size = 20),
          axis.title.x = ggplot2::element_text(size = 15),
          axis.text.y = ggplot2::element_text(size = 15),
          axis.text.x = ggplot2::element_text(size = 12, angle = 75, hjust = 1),
          axis.title.y = ggplot2::element_text(size = 15),
          legend.title = ggplot2::element_text(size = 15),
          legend.text = ggplot2::element_text(size = 15),
          strip.text.x = ggplot2::element_text(size = 15),
          strip.text = ggplot2::element_text(size = 15),
          strip.background = ggplot2::element_blank()
        ) +
        {
          if (protein_abundance_plot == FALSE) ggplot2::scale_color_manual(values = rep(protti_colours, 10)) # repeated colours to have enough even for proteins with many peptides
        } +
        {
          if (protein_abundance_plot == TRUE) ggplot2::scale_color_manual(values = c(rep("gray", length(unique(dplyr::pull(.x, {{ peptide }}))) - 1), "green"))
        }
    }
  )

  plot
}