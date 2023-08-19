#' Check ranked intensities
#'
#' Calculates and plots ranked intensities for proteins, peptides or precursors.
#'
#' @param data a data frame that contains at least sample names, grouping identifiers (precursor,
#' peptide or protein) and log2 transformed intensities for each grouping identifier.
#' @param sample a character column in the `data` data frame that contains the sample names.
#' @param grouping a character column in the `data` data frame that contains protein, precursor,
#' or peptide identifiers.
#' @param intensity_log2 a numeric column in the `data` data frame that contains the log2
#' transformed intensities of the selected grouping variable.
#' @param facet a logical value that specifies whether the calculation should be done group wise by
#' sample and if the resulting plot should be faceted by sample. (default is `FALSE`).
#' If `facet = FALSE` the median of each protein intensity will be returned.
#' @param plot a logical value that specifies whether the result should be plotted (default is `FALSE`).
#' @param y_axis_transformation a character value that determines that y-axis transformation. The
#' value is either "log2" or "log10" (default is "log10").
#' @param interactive a logical value that specifies whether the plot should be interactive
#' (default is `FALSE`).
#'
#' @return A data frame containing the ranked intensities is returned. If `plot = TRUE` a plot
#' is returned. The intensities are log10 transformed for the plot.
#' @import dplyr
#' @import ggplot2
#' @importFrom magrittr %>%
#' @importFrom tidyr drop_na
#' @importFrom plotly ggplotly
#' @importFrom utils data
#' @importFrom stringr str_sort
#' @importFrom ggrepel geom_text_repel
#' @export
#'
#' @examples
#' set.seed(123) # Makes example reproducible
#'
#' # Create synthetic data
#' data <- create_synthetic_data(
#'   n_proteins = 50,
#'   frac_change = 0.05,
#'   n_replicates = 4,
#'   n_conditions = 3,
#'   method = "effect_random",
#'   additional_metadata = FALSE
#' )
#'
#' # Plot ranked intensities for all samples combined
#' qc_ranked_intensities(
#'   data = data,
#'   sample = sample,
#'   grouping = peptide,
#'   intensity_log2 = peptide_intensity,
#'   plot = TRUE,
#' )
#'
#' # Plot ranked intensities for each sample separately
#' qc_ranked_intensities(
#'   data = data,
#'   sample = sample,
#'   grouping = peptide,
#'   intensity_log2 = peptide_intensity,
#'   plot = TRUE,
#'   facet = TRUE
#' )
#'
qc_ranked_intensities <- function(data,
                                  sample,
                                  grouping,
                                  intensity_log2,
                                  facet = FALSE,
                                  plot = FALSE,
                                  y_axis_transformation = "log10",
                                  interactive = FALSE) {
  protti_colours <- "placeholder" # assign a placeholder to prevent a missing global variable warning
  utils::data("protti_colours", envir = environment()) # then overwrite it with real data

  data_prep <- data %>%
    dplyr::distinct({{ sample }}, {{ grouping }}, {{ intensity_log2 }}) %>%
    tidyr::drop_na({{ intensity_log2 }}) %>%
    dplyr::mutate(intensity = 2^{{ intensity_log2 }})

  if (is(dplyr::pull(data, {{ sample }}), "character")) {
    # reorder sample names if they are not already a factor
    data_prep <- data_prep %>%
      dplyr::mutate({{ sample }} := factor({{ sample }},
        levels = unique(stringr::str_sort({{ sample }}, numeric = TRUE))
      ))
  }

  if (facet) {
    input <- data_prep %>%
      dplyr::mutate(log10_intensity = log10(.data$intensity)) %>% # log10 transform intensities
      dplyr::distinct({{ sample }}, {{ grouping }}, {{ intensity_log2 }}, .data$log10_intensity) %>%
      dplyr::arrange(dplyr::desc(.data$log10_intensity)) %>% # sort intensities
      dplyr::group_by({{ sample }}) %>%
      dplyr::mutate(rank = 1:dplyr::n()) %>% # give intensities a rank
      dplyr::ungroup()

    if (y_axis_transformation == "log2") {
      input <- input %>%
        dplyr::mutate(intensity_plot = {{ intensity_log2 }})
    }

    if (y_axis_transformation == "log10") {
      input <- input %>%
        dplyr::mutate(intensity_plot = .data$log10_intensity)
    }

    input_top <- input %>%
      dplyr::group_by({{ sample }}) %>%
      filter(rank <= 10 | rank >= dplyr::n() - 9)
  } else {
    input <- data_prep %>%
      dplyr::group_by({{ grouping }}) %>%
      dplyr::summarize(median_intensity = stats::median(.data$intensity, na.rm = TRUE), .groups = "drop") %>%
      # calculate median intensities for each group
      dplyr::mutate(
        log2_median_intensity = log2(.data$median_intensity), # log2 transform intensities
        log10_median_intensity = log10(.data$median_intensity)
      ) %>% # log10 transform intensities
      dplyr::distinct({{ grouping }}, .data$log2_median_intensity, .data$log10_median_intensity) %>%
      dplyr::arrange(dplyr::desc(.data$log2_median_intensity)) %>% # sort intensities
      dplyr::mutate(rank = 1:dplyr::n()) # give intensities a rank

    if (y_axis_transformation == "log2") {
      input <- input %>%
        dplyr::mutate(intensity_plot = .data$log2_median_intensity)
    }

    if (y_axis_transformation == "log10") {
      input <- input %>%
        dplyr::mutate(intensity_plot = .data$log10_median_intensity)
    }

    input_top <- input %>%
      filter(rank <= 10 | rank >= dplyr::n() - 9)
  }

  if (plot == FALSE) {
    output <- input %>%
      dplyr::select(-"intensity_plot")

    # return data frame
    return(output)
  }

  plot <- input %>%
    ggplot2::ggplot(ggplot2::aes(.data$rank, .data$intensity_plot)) +
    ggplot2::geom_point(
      size = 3,
      colour = protti_colours[1]
    ) +
    {
      if (facet) {
        ggplot2::labs(title = "Ranked Intensitites", x = "Rank", y = paste0(y_axis_transformation, " Intensity"))
      } else {
        ggplot2::labs(title = "Ranked Intensitites", x = "Rank", y = paste0("Median ", y_axis_transformation, " Intensity"))
      }
    } +
    {
      if (facet) {
        ggplot2::facet_wrap(rlang::new_formula(NULL, rlang::enquo(sample)), scales = "free", ncol = 4)
      }
    } +
    {
      if (!interactive) {
        ggrepel::geom_text_repel(
          data = input_top,
          aes(label = {{ grouping }}),
          size = 4
        )
      }
    } +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 20),
      axis.title.x = ggplot2::element_text(size = 15),
      axis.text.y = ggplot2::element_text(size = 15),
      axis.text.x = ggplot2::element_text(size = 12),
      axis.title.y = ggplot2::element_text(size = 15),
      strip.text = ggplot2::element_text(size = 15),
      strip.background = element_blank()
    )

  if (interactive == FALSE) {
    return(plot)
  } else {
    suppressWarnings(plotly::ggplotly(plot))
  }
}
