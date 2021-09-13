#' Peptide abundance profile plot
#'
#' `r lifecycle::badge('deprecated')`
#' This function was deprecated due to its name changing to `peptide_profile_plot()`.
#'
#' @keywords internal
#' @export
plot_peptide_profiles <- function(...) {
  # This function has been renamed and is therefore deprecated.
  lifecycle::deprecate_warn("0.2.0",
    "plot_peptide_profiles()",
    "peptide_profile_plot()",
    details = "This function has been renamed."
  )

  peptide_profile_plot(...)
}
#' Peptide abundance profile plot
#'
#' Creates a plot of peptide abundances across samples. This is helpful to investigate effects of
#' peptide and protein abundance changes in different samples and conditions.
#'
#' @param data a data frame that contains at least the input variables.
#' @param sample a character column in the \code{data} data frame that contains sample names.
#' @param peptide a character column in the \code{data} data frame that contains peptide or
#' precursor names.
#' @param intensity_log2 a numeric column in the \code{data} data frame that contains log2
#' transformed intensities.
#' @param grouping a character column in the \code{data} data frame that contains groups by which
#' the data should be split. This can be for example protein IDs.
#' @param targets a character vector that specifies elements of the grouping column which should
#' be plotted. This can also be \code{"all"} if plots for all groups should be created. Depending
#' on the number of elements in your grouping column this can be many plots.
#' @param protein_abundance_plot a logical value. If the input for this plot comes directly from
#' \code{calculate_protein_abundance} this argument can be set to \code{TRUE}. This displays all
#' peptides in gray, while the protein abundance is displayed in green.
#' @param interactive a logical value that indicates whether the plot should be interactive
#' (default is FALSE). If this is TRUE only one target can be supplied to the function. Interactive
#' plots cannot be exported either.
#' @param export a logical value that indicates if plots should be exported as PDF. The output
#' directory will be the current working directory. The name of the file can be chosen using the
#' \code{export_name} argument.
#' @param export_name A character vector that provides the name of the exported file if
#' \code{export = TRUE}.
#'
#' @return A list of peptide profile plots.
#' @import ggplot2
#' @import progress
#' @importFrom magrittr %>%
#' @importFrom dplyr distinct pull filter
#' @importFrom tidyr drop_na
#' @importFrom rlang !! ensym
#' @importFrom plotly ggplotly
#' @importFrom purrr map2
#' @importFrom utils data
#' @export
#'
#' @examples
#' # Create example data
#' data <- data.frame(
#'   sample = c(
#'     rep("S1", 6),
#'     rep("S2", 6),
#'     rep("S1", 2),
#'     rep("S2", 2)
#'   ),
#'   protein_id = c(
#'     rep("P1", 12),
#'     rep("P2", 4)
#'   ),
#'   precursor = c(
#'     rep(c("A1", "A2", "B1", "B2", "C1", "D1"), 2),
#'     rep(c("E1", "F1"), 2)
#'   ),
#'   peptide = c(
#'     rep(c("A", "A", "B", "B", "C", "D"), 2),
#'     rep(c("E", "F"), 2)
#'   ),
#'   intensity = c(
#'     rnorm(n = 6, mean = 15, sd = 2),
#'     rnorm(n = 6, mean = 21, sd = 1),
#'     rnorm(n = 2, mean = 15, sd = 1),
#'     rnorm(n = 2, mean = 15, sd = 2)
#'   )
#' )
#'
#' # Calculate protein abundances and retain precursor
#' # abundances that can be used in a peptide profile plot
#' complete_abundances <- calculate_protein_abundance(
#'   data,
#'   sample = sample,
#'   protein_id = protein_id,
#'   precursor = precursor,
#'   peptide = peptide,
#'   intensity_log2 = intensity,
#'   method = "iq",
#'   for_plot = TRUE
#' )
#'
#' # Plot protein abundance profile
#' # protein_abundance_plot can be set to
#' # FALSE to to also colour precursors
#' peptide_profile_plot(
#'   data = complete_abundances,
#'   sample = sample,
#'   peptide = precursor,
#'   intensity_log2 = intensity,
#'   grouping = protein_id,
#'   targets = c("P1"),
#'   protein_abundance_plot = TRUE
#' )
peptide_profile_plot <- function(data,
                                 sample,
                                 peptide,
                                 intensity_log2,
                                 grouping,
                                 targets,
                                 protein_abundance_plot = FALSE,
                                 interactive = FALSE,
                                 export = FALSE,
                                 export_name = "peptide_profile_plots") {
  . <- NULL
  n_samples <- length(unique(dplyr::pull(data, {{ sample }})))
  protti_colours <- "placeholder" # assign a placeholder to prevent a missing global variable warning
  utils::data("protti_colours", envir = environment()) # then overwrite it with real data
  if (missing(targets)) stop("Please provide at least one target to plot!")
  if (!("all" %in% targets)) {
    input <- data %>%
      dplyr::distinct({{ sample }}, {{ peptide }}, {{ intensity_log2 }}, {{ grouping }}) %>%
      tidyr::drop_na({{ intensity_log2 }}) %>%
      dplyr::filter({{ grouping }} %in% targets) %>%
      split(dplyr::pull(., !!ensym(grouping)))
  }
  if ("all" %in% targets) {
    groups <- length(unique(dplyr::pull(data, {{ grouping }})))
    message("Splitting into ", groups, " groups and returning ", groups, " plots.")
    input <- data %>%
      dplyr::distinct({{ sample }}, {{ peptide }}, {{ intensity_log2 }}, {{ grouping }}) %>%
      tidyr::drop_na({{ intensity_log2 }}) %>%
      split(dplyr::pull(., !!ensym(grouping)))
  }
  pb <- progress::progress_bar$new(
    total = length(input),
    format = " Creating plots [:bar] :current/:total (:percent) :eta"
  )
  plot <- purrr::map2(
    .x = input,
    .y = names(input),
    .f = ~ {
      pb$tick()
      ggplot2::ggplot(.x, ggplot2::aes({{ sample }}, {{ intensity_log2 }}, group = {{ peptide }}, color = {{ peptide }})) +
        ggplot2::geom_point() +
        ggplot2::geom_line(size = 1) +
        ggplot2::labs(
          title = paste("Peptide profiles:", .y),
          x = "Sample",
          y = "Intensity [log2]",
          col = "Peptides"
        ) +
        ggplot2::theme_bw() +
        {
          if (length(unique(dplyr::pull(.x, {{ peptide }}))) > 25) ggplot2::theme(legend.position = "none")
        } +
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
          # repeated colours to have enough even for proteins with many peptides
          if (protein_abundance_plot == FALSE) ggplot2::scale_color_manual(values = rep(protti_colours, 10))
        } +
        {
          if (protein_abundance_plot == TRUE) {
            ggplot2::scale_color_manual(values = c(rep("gray", length(unique(dplyr::pull(.x, {{ peptide }}))) - 1), "green"))
          }
        }
    }
  )
  if (interactive == FALSE) {
    if (export == TRUE) {
      grDevices::pdf(
        file = paste0(export_name, ".pdf"),
        width = 10 * ceiling(n_samples / 10), # more samples need wider plot
        height = 6
      )
      suppressWarnings(print(plot))
      grDevices::dev.off()
    } else {
      return(plot)
    }
  }
  if (interactive == TRUE) {
    if (length(targets) > 1) {
      stop(strwrap("Please only provide one target for interactive plots!",
        prefix = "\n", initial = ""
      ))
    }
    if (export == TRUE) {
      stop(strwrap("Interactive plots cannot be exported! Please decide if you
want an interactive plot or if you want to export your plots.",
        prefix = "\n", initial = ""
      ))
    }
    plotly::ggplotly(plot[[1]], tooltip = c("x", "y", "group"))
  }
}
